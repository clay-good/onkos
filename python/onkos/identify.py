"""Practical identifiability — could a realistic trial even estimate this parameter?

Onkos surfaces that the kill/resistance terms carry ~90% IIV CV, and the spec's
stated reason is that *resistance is poorly identifiable from short trials* (spec
§1, §4). This module measures that claim instead of asserting it: given a model and
a realistic clinical observation schedule, it computes — from the **Fisher
information of the design** — the best precision (Cramér-Rao lower bound) with which
each structural parameter could be estimated, and flags the parameters a trial of
that shape cannot pin down.

The headline pairing is the stored ``iiv_cv_percent`` next to the predicted
``rse_percent``: a parameter that is both high-CV *and* high-RSE is one whose large
reported variability is, at least in part, a flat-likelihood artifact of the
originating design — not a clean estimate of biological spread. That is exactly the
parameter whose out-of-context transport (the project's load-bearing risk) is least
defensible.

This is the design-level companion to the parameter-level sensitivity triage
(:mod:`onkos.sensitivity`) and the model-level model-selection triage
(:mod:`onkos.combine`): sensitivity asks *which* uncertainty drives the forecast;
identifiability asks which uncertainty is *reducible by better data* versus
*structurally stuck under this design*.

Population / design level only. NOT a per-patient quantity, NOT a prognosis, NOT a
therapy ranking. "Unidentifiable" always means "under this schedule and residual-
error model", and the result says so. A singular design reports ``inf`` (honest),
never a fabricated finite bound; identifiability never moves a tier.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .export.registry import get_kernel, kernel_values
from .load import Dataset
from .simulate import simulate

__all__ = [
    "Identifiability",
    "ParamIdentifiability",
    "fisher_information",
    "crlb_rse",
    "collinearity_index",
    "identifiability",
]

# A realistic RECIST-style scan cadence: baseline then ~q6w for a year (weeks).
DEFAULT_SCHEDULE = (0.0, 6.0, 12.0, 18.0, 24.0, 36.0, 48.0)
# Declared ceilings (conventions, printed with every verdict — not fitted).
RSE_CEILING = 0.50  # a parameter with predicted RSE > 50% is practically unidentifiable
COLLINEARITY_CEILING = 15.0  # γ_K above this ⇒ a confounded (non-identifiable) combination
_CV_ARTIFACT_THRESHOLD = 50.0  # an IIV CV at/above this, paired with a failing RSE, is flagged


# --------------------------------------------------------------------------- #
# Pure information algebra — the landmark-tested estimator core.              #
# The input ``scaled_sens`` is the residual-weighted sensitivity matrix       #
# S̃_ij = (∂f(t_i)/∂θ_j) / σ_i   (rows = observations, cols = parameters),     #
# so the Fisher information is simply S̃ᵀ S̃.                                  #
# --------------------------------------------------------------------------- #


def fisher_information(scaled_sens) -> np.ndarray:
    """Design Fisher-information matrix ``M = S̃ᵀ S̃`` from the residual-weighted
    sensitivity matrix (rows = observations, cols = parameters). Information adds
    over independent observations, which is exactly the row-sum of outer products."""
    s = np.atleast_2d(np.asarray(scaled_sens, dtype=float))
    return s.T @ s


def crlb_rse(scaled_sens, theta) -> tuple:
    """Cramér-Rao precision bound from the residual-weighted sensitivity matrix.

    Returns ``(cov, se, rse)`` where ``cov = M⁻¹`` (the lower-bound parameter
    covariance), ``se = sqrt(diag(cov))``, and ``rse = se / |theta|`` (a *fraction*).
    A non-invertible (structurally non-identifiable) design yields ``inf`` standard
    errors and RSEs — reported honestly, never silently regularized.
    """
    theta = np.abs(np.asarray(theta, dtype=float))
    m = fisher_information(scaled_sens)
    p = m.shape[0]
    # Singular / ill-conditioned FIM -> the design does not identify the parameters.
    if not np.all(np.isfinite(m)) or np.linalg.matrix_rank(m) < p:
        inf = np.full(p, np.inf)
        return np.full((p, p), np.inf), inf, inf
    cov = np.linalg.inv(m)
    var = np.diag(cov)
    se = np.where(var >= 0, np.sqrt(np.abs(var)), np.inf)
    with np.errstate(divide="ignore", invalid="ignore"):
        rse = np.where(theta > 0, se / theta, np.inf)
    return cov, se, rse


def collinearity_index(scaled_sens) -> float:
    """Brun-Reichert-Künsch collinearity index ``γ_K = 1/sqrt(λ_min(S̃ᵀS̃))`` over
    the *column-normalized* weighted sensitivity matrix.

    ``γ_K = 1`` ⇔ orthogonal (perfectly separable parameter directions); ``γ_K → ∞``
    ⇔ collinear (a non-identifiable parameter combination). Scale-free: rescaling a
    parameter's column leaves ``γ_K`` unchanged. A zero or parallel column yields
    ``inf``."""
    s = np.atleast_2d(np.asarray(scaled_sens, dtype=float))
    norms = np.linalg.norm(s, axis=0)
    if np.any(norms == 0) or not np.all(np.isfinite(s)):
        return float("inf")
    sn = s / norms  # unit-length columns -> Gram matrix has unit diagonal
    eigs = np.linalg.eigvalsh(sn.T @ sn)
    lam_min = float(eigs[0])
    if lam_min <= 0:
        return float("inf")
    return float(1.0 / np.sqrt(lam_min))


# --------------------------------------------------------------------------- #
# The result objects.                                                         #
# --------------------------------------------------------------------------- #


@dataclass
class ParamIdentifiability:
    symbol: str  # record-facing parameter symbol (e.g. "lambda")
    central: float  # central parameter value the analysis is local to
    rse_percent: float  # predicted relative standard error (Cramér-Rao), in %
    iiv_cv_percent: float | None  # the stored IIV CV (for the CV-vs-RSE pairing)
    identifiable: bool  # rse_percent < rse_ceiling


@dataclass
class Identifiability:
    record_id: str
    schedule: np.ndarray
    sigma_prop: float
    sigma_add: float
    tier: str
    collinearity_index: float
    rse_ceiling_percent: float
    collinearity_ceiling: float
    params: list = field(default_factory=list)  # sorted desc by rse_percent (worst first)
    warnings: list = field(default_factory=list)
    clinical_use: str = CLINICAL_USE

    @property
    def practically_identifiable(self) -> bool:
        """Every parameter precise enough AND no near-non-identifiable combination,
        under this design. The conjunction so neither failure mode hides the other."""
        per_param = all(p.identifiable for p in self.params)
        return bool(per_param and self.collinearity_index < self.collinearity_ceiling)

    @property
    def worst(self) -> ParamIdentifiability | None:
        """The least-identifiable parameter — the curation-triage target."""
        return self.params[0] if self.params else None

    def to_dict(self) -> dict:
        """JSON-serializable result. Carries the clinical-use prohibition and names
        its design (schedule + residual error) — "unidentifiable" is always
        relative to a stated design, never an in-principle claim."""
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "record_id": self.record_id,
            "design": {
                "schedule_weeks": np.asarray(self.schedule, dtype=float).tolist(),
                "sigma_prop": self.sigma_prop,
                "sigma_add": self.sigma_add,
            },
            "tier": self.tier,
            "practically_identifiable": self.practically_identifiable,
            "collinearity_index": self.collinearity_index,
            "rse_ceiling_percent": self.rse_ceiling_percent,
            "collinearity_ceiling": self.collinearity_ceiling,
            "parameters": [
                {
                    "symbol": p.symbol,
                    "central": p.central,
                    "rse_percent": p.rse_percent,
                    "iiv_cv_percent": p.iiv_cv_percent,
                    "identifiable": p.identifiable,
                }
                for p in self.params
            ],
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


# --------------------------------------------------------------------------- #
# Binding the information algebra to a record + a clinical design.            #
# --------------------------------------------------------------------------- #


def _tumor_at(ds, record_id, *, context, drug_effect, exposure, exposure_response, t, overrides):
    return simulate(
        ds, record_id, context=context, drug_effect=drug_effect, exposure=exposure,
        exposure_response=exposure_response, t=t, param_overrides=overrides,
    ).tumor_size


def identifiability(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    schedule=DEFAULT_SCHEDULE,
    drug_effect: float | None = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    sigma_prop: float = 0.2,
    sigma_add: float = 0.0,
    h: float = 1e-4,
    rse_ceiling: float = RSE_CEILING,
    collinearity_ceiling: float = COLLINEARITY_CEILING,
) -> Identifiability:
    """Predict, from the design Fisher information, how precisely each structural
    parameter of ``record_id`` could be estimated under ``schedule`` (the scan times,
    in the kernel's time unit) and a combined residual-error model
    ``σ_i = sqrt(sigma_add² + (sigma_prop·y_i)²)``.

    Returns an :class:`Identifiability` with per-parameter predicted RSE paired with
    the stored IIV CV, the collinearity index, and the ``practically_identifiable``
    verdict. Restricted to dynamic (ODE) TGI/growth records — identifiability of a
    *trajectory* is what the question asks; survival and transform kernels are not
    trajectories in this sense.
    """
    record = ds[record_id]
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError(
            f"identifiability is defined for dynamic (ODE) TGI/growth records; "
            f"'{record_id}' binds a '{spec.kind}' kernel"
        )
    t = np.asarray(schedule, dtype=float)
    if t.size < 1:
        raise ValueError("schedule must contain at least one observation time")

    base = kernel_values(record)
    sim_kw = dict(context=context, drug_effect=drug_effect, exposure=exposure,
                  exposure_response=exposure_response, t=t)
    central_tr = simulate(ds, record_id, **sim_kw)
    y = central_tr.tumor_size

    # Combined residual-error standard deviation per observation (declared model).
    sigma = np.sqrt(sigma_add**2 + (sigma_prop * y) ** 2)
    sigma = np.where(sigma > 0, sigma, 1e-12)  # guard the (rare) y_i == 0 point

    kparams = list(spec.params)  # kernel-internal names (perturbation targets)
    symbols = list(spec.record_symbols)  # record-facing symbols (display + CV lookup)
    cv_by_symbol = {p.symbol: p.iiv_cv_percent for p in record.parameters}

    # Local sensitivity matrix by central finite difference over the kernel.
    sens = np.zeros((t.size, len(kparams)))
    theta = np.zeros(len(kparams))
    for j, kname in enumerate(kparams):
        theta0 = float(base[kname])
        theta[j] = theta0
        step = h * abs(theta0) if theta0 != 0 else h
        up = _tumor_at(ds, record_id, overrides={kname: theta0 + step}, **sim_kw)
        dn = _tumor_at(ds, record_id, overrides={kname: theta0 - step}, **sim_kw)
        sens[:, j] = (up - dn) / (2.0 * step)

    scaled = sens / sigma[:, None]
    _, _, rse = crlb_rse(scaled, theta)
    gamma = collinearity_index(scaled)

    params = [
        ParamIdentifiability(
            symbol=symbols[j],
            central=float(theta[j]),
            rse_percent=float(rse[j] * 100.0),
            iiv_cv_percent=cv_by_symbol.get(symbols[j]),
            identifiable=bool(rse[j] < rse_ceiling),
        )
        for j in range(len(kparams))
    ]
    # Worst (least identifiable) first — the triage order.
    params.sort(key=lambda p: (p.rse_percent if np.isfinite(p.rse_percent) else np.inf),
                reverse=True)

    warnings: list[str] = []
    if any(not np.isfinite(p.rse_percent) for p in params) or not np.isfinite(gamma):
        warnings.append(
            "singular_fim: the design does not identify the parameters "
            "(non-invertible Fisher information); RSE/collinearity reported as inf"
        )
    if gamma >= collinearity_ceiling and np.isfinite(gamma):
        warnings.append(
            f"collinear_parameters: collinearity index {gamma:.1f} >= {collinearity_ceiling:.0f} "
            "— a parameter combination is poorly separable under this design"
        )
    for p in params:
        cv = p.iiv_cv_percent
        if cv is not None and cv >= _CV_ARTIFACT_THRESHOLD and not p.identifiable:
            warnings.append(
                f"cv_is_identifiability_artifact: '{p.symbol}' carries IIV CV {cv:.0f}% and is "
                f"practically unidentifiable (predicted RSE {p.rse_percent:.0f}%) under this "
                "design — its large reported variability is at least partly a flat-likelihood "
                "artifact, not a clean estimate of biological spread"
            )

    return Identifiability(
        record_id=record_id,
        schedule=t,
        sigma_prop=sigma_prop,
        sigma_add=sigma_add,
        tier=central_tr.tier,  # the record's propagated tier — analysis cannot move it
        collinearity_index=gamma,
        rse_ceiling_percent=rse_ceiling * 100.0,
        collinearity_ceiling=collinearity_ceiling,
        params=params,
        warnings=warnings,
    )
