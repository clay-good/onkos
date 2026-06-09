"""Variance-based parameter sensitivity — *which* uncertainty matters most.

Uncertainty propagation (:mod:`onkos.uncertainty`) gives prediction *bands*; this
module attributes that variance to individual parameters. Because the dataset
samples each IIV-bearing parameter independently (lognormal), a first-order
decomposition is exact-in-form: the standardized regression coefficient of a
parameter equals its correlation with the target, and the squared coefficients
partition the explained variance.

The point is curation triage. The spec's highest-leverage contribution is
verifying records against the source PDF (§9); a sensitivity ranking says which
parameter's uncertainty actually moves the survival prediction — i.e. where that
verification effort pays off first.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .export.registry import get_kernel, kernel_values
from .load import Dataset
from .simulate import simulate
from .uncertainty import _iiv_by_kernel_name

__all__ = ["ParamSensitivity", "SensitivityResult", "sensitivity"]


@dataclass
class ParamSensitivity:
    symbol: str  # record-facing parameter symbol (e.g. "lambda")
    iiv_cv_percent: float
    src: float  # standardized regression coefficient (signed; direction of effect)
    contribution: float  # SRC^2 normalized to the explained variance (0..1)


@dataclass
class SensitivityResult:
    record_id: str
    target: str
    n: int
    n_used: int
    r_squared: float  # first-order variance of the target explained by the inputs
    indices: list = field(default_factory=list)  # sorted desc by contribution

    @property
    def dominant(self) -> ParamSensitivity | None:
        return self.indices[0] if self.indices else None


def _target_value(traj, target: str) -> float:
    if target in ("median_os_weeks", "median_os"):
        m = traj.median_os
        return float("nan") if m is None else float(m)
    if target in ("median_pfs_weeks", "median_pfs"):
        m = traj.median_pfs
        return float("nan") if m is None else float(m)
    return float(traj.metrics.get(target, float("nan")))


def sensitivity(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    target: str = "median_os_weeks",
    drug_effect: float | None = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    t: np.ndarray | None = None,
    n: int = 400,
    seed: int = 0,
) -> SensitivityResult:
    """Rank the IIV-bearing parameters of ``record_id`` by how much their
    variability drives ``target`` (a metric key or ``"median_os_weeks"``)."""
    if t is None:
        t = np.linspace(0.0, 156.0, 313)
    t = np.asarray(t, dtype=float)
    record = ds[record_id]
    spec = get_kernel(record)
    cvs = _iiv_by_kernel_name(record)  # kernel-internal name -> CV (fraction)
    if not cvs:
        raise ValueError(f"record '{record_id}' has no IIV parameters to analyze")

    kernel_to_symbol = dict(zip(spec.params, spec.record_symbols))
    cv_by_symbol = {p.symbol: p.iiv_cv_percent for p in record.parameters if p.iiv_cv_percent}
    base = kernel_values(record)
    names = list(cvs)  # kernel-internal names, fixed order
    rng = np.random.default_rng(seed)

    draws = np.zeros((n, len(names)))  # centered log-deviations (the sampled z)
    y = np.zeros(n)
    for i in range(n):
        overrides = {}
        for j, kname in enumerate(names):
            sd = np.sqrt(np.log(1.0 + cvs[kname] ** 2))
            z = float(rng.normal(0.0, sd))
            draws[i, j] = z
            overrides[kname] = base[kname] * np.exp(z)
        traj = simulate(
            ds, record_id, context=context, drug_effect=drug_effect, exposure=exposure,
            exposure_response=exposure_response, t=t, param_overrides=overrides,
        )
        y[i] = _target_value(traj, target)

    mask = np.isfinite(y)
    X, yv = draws[mask], y[mask]
    n_used = int(mask.sum())

    src = np.zeros(len(names))
    if n_used >= 3 and np.std(yv) > 0:
        ys = (yv - yv.mean()) / yv.std()
        for j in range(len(names)):
            col = X[:, j]
            if np.std(col) > 0:
                src[j] = float(np.corrcoef(col, ys)[0, 1])  # = SRC for independent inputs

    ss = float(np.sum(src**2))
    r_squared = min(ss, 1.0)
    indices = [
        ParamSensitivity(
            symbol=kernel_to_symbol[kname],
            iiv_cv_percent=cv_by_symbol.get(kernel_to_symbol[kname], cvs[kname] * 100),
            src=float(src[j]),
            contribution=float(src[j] ** 2 / ss) if ss > 0 else 0.0,
        )
        for j, kname in enumerate(names)
    ]
    indices.sort(key=lambda p: p.contribution, reverse=True)
    return SensitivityResult(
        record_id=record_id, target=target, n=n, n_used=n_used,
        r_squared=r_squared, indices=indices,
    )
