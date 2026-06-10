"""D-optimal trial design — the best sampling schedule a fixed budget allows.

:mod:`onkos.identify` (v0.22) evaluates whether a *given* observation schedule can
estimate a model's parameters, from the design Fisher information ``M = S̃ᵀS̃``. This
module answers the question that left open: a pharmacometrician *chooses* the sampling
times, so what is the **best** schedule under a fixed budget of ``N`` measurements, and
does even it cross the identifiability line?

The information is **additive over timepoints** — ``M = Σᵢ s̃ᵢ s̃ᵢᵀ`` — so the
sensitivity rows are computed once on a dense candidate grid and schedule optimization is
pure linear algebra over row subsets (no re-simulation per candidate). The shipped
criterion is **D-optimality**: maximize ``det M`` (minimize the volume of the joint
confidence ellipsoid), the textbook optimal-design objective. The combinatorial subset
search uses greedy forward selection (baseline-anchored), and the reported optimal is the
better of the greedy and uniform designs, so ``D-efficiency ≥ 1`` by construction.

The honest payload: for the Claret resistance model the optimal design lowers every RSE
yet still cannot identify ``kL`` / ``λ`` — their flatness is *structural*, not a design
failure — while for a 2-parameter biexponential the optimal design does identify both. So
optimal design is the control that separates "badly designed trial" from "structurally
unidentifiable parameter."

Design / population level only. NOT a per-patient schedule, NOT a prognosis, NOT a dosing
or therapy choice. Greedy selection is a heuristic (guaranteed only ``≥ uniform``); the
result says so. A design analysis never moves a tier.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .export.registry import get_kernel, kernel_values
from .identify import (
    RSE_CEILING,
    collinearity_index,
    crlb_rse,
    fisher_information,
)
from .load import Dataset
from .simulate import simulate

__all__ = [
    "d_optimal_rows",
    "CandidateDesign",
    "OptimalDesign",
    "optimal_schedule",
]

_RIDGE = 1e-9  # tiny regularizer so log-det is finite during the under-determined greedy phase


# --------------------------------------------------------------------------- #
# Pure selection core — D-optimal row selection over a scaled-sensitivity      #
# matrix. Landmark-tested in isolation (no record needed).                     #
# --------------------------------------------------------------------------- #


def _logdet_fim(scaled_rows: np.ndarray, p: int) -> float:
    """Log-determinant of the (ridge-stabilized) Fisher information of the given rows."""
    m = fisher_information(scaled_rows) + _RIDGE * np.eye(p)
    sign, ld = np.linalg.slogdet(m)
    return float(ld) if sign > 0 else float("-inf")


def d_optimal_rows(scaled_sens, n: int, seed_rows=()) -> list:
    """Greedy forward D-optimal selection of ``n`` rows from ``scaled_sens`` (the
    residual-weighted sensitivity matrix, one row per candidate timepoint): repeatedly add
    the row that most increases ``log det(Σ rows rowᵀ)``, starting from ``seed_rows`` (e.g.
    the mandatory baseline). Returns the chosen row indices in selection order.

    Greedy is a standard D-optimal heuristic, not a global optimum — the caller guarantees
    only ``≥ uniform`` by comparison."""
    ss = np.asarray(scaled_sens, dtype=float)
    if ss.ndim != 2:
        raise ValueError("scaled_sens must be a 2-D (timepoints x parameters) matrix")
    m, p = ss.shape
    if n < 1:
        raise ValueError("n must be at least 1")
    if n > m:
        raise ValueError(f"n={n} exceeds the {m} candidate rows")
    chosen = list(dict.fromkeys(int(r) for r in seed_rows))  # dedup, preserve order
    if len(chosen) > n:
        raise ValueError("more seed_rows than the budget n")
    while len(chosen) < n:
        remaining = [c for c in range(m) if c not in chosen]
        best = max(remaining, key=lambda c: _logdet_fim(ss[chosen + [c]], p))
        chosen.append(best)
    return chosen


# --------------------------------------------------------------------------- #
# Binding the selection core to a record + a measurement budget.              #
# --------------------------------------------------------------------------- #


@dataclass
class CandidateDesign:
    label: str  # "uniform" | "D-optimal"
    schedule: list  # the sampling times (in the kernel's time unit)
    rse_percent: dict  # symbol -> predicted Cramér-Rao RSE %
    max_rse_percent: float
    collinearity_index: float
    log_det_fim: float
    identifiable: list  # symbols with RSE below the ceiling

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "schedule": [float(t) for t in self.schedule],
            "rse_percent": self.rse_percent,
            "max_rse_percent": self.max_rse_percent,
            "collinearity_index": self.collinearity_index,
            "log_det_fim": self.log_det_fim,
            "identifiable": list(self.identifiable),
        }


@dataclass
class OptimalDesign:
    record_id: str
    context: dict
    tier: str
    n_samples: int
    horizon: float
    criterion: str  # "D"
    rse_ceiling_percent: float
    uniform: CandidateDesign
    optimal: CandidateDesign
    d_efficiency: float  # (det M_opt / det M_uni)^(1/p) ≥ 1 — how much more informative
    rescues_any: bool  # did the optimal design identify a parameter uniform could not?
    structurally_flat: list  # parameters unidentifiable even under the optimal design
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "record_id": self.record_id,
            "context": self.context,
            "tier": self.tier,
            "n_samples": self.n_samples,
            "horizon": self.horizon,
            "criterion": self.criterion,
            "rse_ceiling_percent": self.rse_ceiling_percent,
            "uniform": self.uniform.to_dict(),
            "optimal": self.optimal.to_dict(),
            "d_efficiency": self.d_efficiency,
            "rescues_any": self.rescues_any,
            "structurally_flat": list(self.structurally_flat),
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def _scaled_sensitivities(
    ds, record_id, *, context, drug_effect, exposure, exposure_response, grid, sigma_prop,
    sigma_add, h,
):
    """Residual-weighted parameter-sensitivity matrix on ``grid`` (one row per time), plus
    the parameter central values and record-facing symbols. Same finite-difference and
    error model as :func:`onkos.identify.identifiability`."""
    spec = get_kernel(ds[record_id])
    base = kernel_values(ds[record_id])
    kparams = list(spec.params)
    symbols = list(spec.record_symbols)
    sim_kw = dict(context=context, drug_effect=drug_effect, exposure=exposure,
                  exposure_response=exposure_response, t=grid)

    y = simulate(ds, record_id, **sim_kw).tumor_size
    sigma = np.sqrt(sigma_add**2 + (sigma_prop * y) ** 2)
    sigma = np.where(sigma > 0, sigma, 1e-12)

    sens = np.zeros((grid.size, len(kparams)))
    theta = np.zeros(len(kparams))
    for j, kname in enumerate(kparams):
        t0 = float(base[kname])
        theta[j] = t0
        step = h * abs(t0) if t0 != 0 else h
        up = simulate(ds, record_id, param_overrides={kname: t0 + step}, **sim_kw).tumor_size
        dn = simulate(ds, record_id, param_overrides={kname: t0 - step}, **sim_kw).tumor_size
        sens[:, j] = (up - dn) / (2.0 * step)

    return sens / sigma[:, None], theta, symbols


def _score(ss_rows, theta, symbols, label, schedule, rse_ceiling):
    """Build a CandidateDesign from the chosen scaled-sensitivity rows."""
    _, _, rse = crlb_rse(ss_rows, theta)
    gamma = collinearity_index(ss_rows)
    rse_pct = {symbols[j]: float(rse[j] * 100.0) for j in range(len(symbols))}
    identifiable = [symbols[j] for j in range(len(symbols)) if rse[j] < rse_ceiling]
    finite = [v for v in rse_pct.values() if np.isfinite(v)]
    return CandidateDesign(
        label=label,
        schedule=list(schedule),
        rse_percent=rse_pct,
        max_rse_percent=float(max(rse_pct.values())) if finite and len(finite) == len(rse_pct)
        else float("inf"),
        collinearity_index=float(gamma),
        log_det_fim=_logdet_fim(ss_rows, len(symbols)),
        identifiable=identifiable,
    )


def optimal_schedule(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    n_samples: int = 7,
    horizon: float = 48.0,
    grid_points: int = 49,
    drug_effect: float | None = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    sigma_prop: float = 0.2,
    sigma_add: float = 0.0,
    h: float = 1e-4,
    rse_ceiling: float = RSE_CEILING,
    include_baseline: bool = True,
) -> OptimalDesign:
    """The D-optimal sampling schedule of ``n_samples`` measurements over ``[0, horizon]``
    for the dynamic (ODE) TGI record ``record_id``, scored against a uniform schedule of the
    same budget. Reuses the :mod:`onkos.identify` Fisher-information core; reports the
    per-parameter Cramér-Rao RSE, the collinearity index, the D-efficiency over uniform, and
    whether the optimal design rescues any parameter uniform could not.

    Design / population level only; restricted to ODE records (a survival/transform record is
    not a trajectory to design for)."""
    spec = get_kernel(ds[record_id])
    if spec.kind != "ode":
        raise ValueError(
            f"optimal_schedule is defined for dynamic (ODE) TGI/growth records; "
            f"'{record_id}' binds a '{spec.kind}' kernel"
        )
    p = len(list(spec.params))
    if n_samples < p:
        raise ValueError(
            f"n_samples={n_samples} is below the {p} parameters; the design is rank-deficient"
        )
    if grid_points < n_samples:
        raise ValueError("grid_points must be at least n_samples")

    grid = np.linspace(0.0, float(horizon), int(grid_points))
    central = simulate(ds, record_id, context=context, drug_effect=drug_effect,
                       exposure=exposure, exposure_response=exposure_response, t=grid)
    ss, theta, symbols = _scaled_sensitivities(
        ds, record_id, context=context, drug_effect=drug_effect, exposure=exposure,
        exposure_response=exposure_response, grid=grid, sigma_prop=sigma_prop,
        sigma_add=sigma_add, h=h,
    )

    # Uniform comparator: n_samples evenly spaced grid indices.
    uniform_rows = sorted(set(np.linspace(0, grid.size - 1, n_samples).round().astype(int)))
    # D-optimal: greedy, baseline-anchored, then take the better of greedy / uniform.
    seed = (0,) if include_baseline else ()
    greedy_rows = sorted(d_optimal_rows(ss, n_samples, seed_rows=seed))
    optimal_rows = (greedy_rows
                    if _logdet_fim(ss[greedy_rows], p) >= _logdet_fim(ss[uniform_rows], p)
                    else uniform_rows)

    uniform = _score(ss[uniform_rows], theta, symbols, "uniform", grid[uniform_rows], rse_ceiling)
    optimal = _score(ss[optimal_rows], theta, symbols, "D-optimal", grid[optimal_rows], rse_ceiling)

    d_eff = float(np.exp((optimal.log_det_fim - uniform.log_det_fim) / p))
    rescues_any = bool(set(optimal.identifiable) - set(uniform.identifiable))
    structurally_flat = [s for s in symbols if s not in optimal.identifiable]

    warnings: list[str] = []
    warnings.append(
        "greedy_heuristic: the optimal schedule is a greedy D-optimal design, guaranteed "
        "only no worse than uniform (D-efficiency ≥ 1), not a proven global optimum"
    )
    if structurally_flat:
        verb = "remains" if len(structurally_flat) == 1 else "remain"
        warnings.append(
            f"structurally_flat: {', '.join(structurally_flat)} {verb} unidentifiable "
            f"(RSE ≥ {rse_ceiling * 100:.0f}%) even under the D-optimal design — a structural "
            "flat, not a design failure"
        )

    return OptimalDesign(
        record_id=record_id,
        context=context or {},
        tier=central.tier,
        n_samples=n_samples,
        horizon=float(horizon),
        criterion="D",
        rse_ceiling_percent=rse_ceiling * 100.0,
        uniform=uniform,
        optimal=optimal,
        d_efficiency=d_eff,
        rescues_any=rescues_any,
        structurally_flat=structurally_flat,
        warnings=warnings,
    )
