"""Monte-Carlo parameter-uncertainty propagation.

The dataset stores inter-individual variability (``iiv_cv_percent``) on its
high-uncertainty kill/resistance terms precisely so it *cannot* masquerade as a
point estimate (spec §4). This module makes that stored uncertainty flow into the
simulation: parameters with an IIV CV are sampled lognormally (the standard
pharmacometric convention) and the resulting tumor-size, TGI-metric, and
population-OS distributions are summarized as bands.

NOT a per-patient prognosis. The bands describe how the *population/trial-level*
prediction moves under the parameters' reported variability — another axis of the
project's "honest about uncertainty by default" stance, complementing the
model-selection uncertainty of the virtual-trial divergence view.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .export.registry import get_kernel, kernel_values
from .load import Dataset
from .simulate import median_survival, simulate

__all__ = ["Ensemble", "simulate_ensemble"]


@dataclass
class Band:
    median: np.ndarray
    lo: np.ndarray
    hi: np.ndarray


@dataclass
class Ensemble:
    record_id: str
    t: np.ndarray
    n: int
    ci: tuple
    tier: str
    tumor_size: Band
    metrics: dict  # name -> {"median","lo","hi"}
    warnings: list = field(default_factory=list)
    os_curve: Band | None = None

    @property
    def median_os(self):
        return self.metrics.get("median_os_weeks")


def _iiv_by_kernel_name(record) -> dict:
    """Map kernel-internal parameter name -> IIV CV (fraction) for sampled params."""
    spec = get_kernel(record)
    sym_to_kernel = dict(zip(spec.record_symbols, spec.params))
    out = {}
    for p in record.parameters:
        if p.iiv_cv_percent and p.symbol in sym_to_kernel:
            out[sym_to_kernel[p.symbol]] = p.iiv_cv_percent / 100.0
    return out


def _band(samples: np.ndarray, ci: tuple) -> Band:
    lo, hi = np.percentile(samples, [ci[0], ci[1]], axis=0)
    return Band(median=np.median(samples, axis=0), lo=lo, hi=hi)


def simulate_ensemble(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    drug_effect: float | None = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
    n: int = 200,
    seed: int = 0,
    ci: tuple = (5.0, 95.0),
) -> Ensemble:
    """Run ``n`` simulations with parameters sampled from their lognormal IIV and
    return median + ``ci`` percentile bands on tumor size, OS, and the metrics.

    Parameters without an IIV CV are held at their central value, so a record with
    no reported variability yields a degenerate (zero-width) band.
    """
    if t is None:
        t = np.linspace(0.0, 104.0, 209)
    t = np.asarray(t, dtype=float)
    record = ds[record_id]
    base = kernel_values(record)
    cvs = _iiv_by_kernel_name(record)
    rng = np.random.default_rng(seed)

    def run(overrides):
        return simulate(
            ds, record_id, context=context, drug_effect=drug_effect, exposure=exposure,
            exposure_response=exposure_response, t=t, survival_link=survival_link,
            param_overrides=overrides,
        )

    central = run(None)  # tier + warnings are parameter-independent

    tumor_samples = np.empty((n, t.size))
    has_os = central.os_curve is not None
    os_samples = np.empty((n, t.size)) if has_os else None
    metric_samples: dict = {k: np.empty(n) for k in central.metrics}
    median_os_samples = np.empty(n)

    for i in range(n):
        overrides = {}
        for kname, cv in cvs.items():
            sd = np.sqrt(np.log(1.0 + cv * cv))  # lognormal: preserves the median
            overrides[kname] = base[kname] * float(np.exp(rng.normal(0.0, sd)))
        tr = run(overrides)
        tumor_samples[i] = tr.tumor_size
        if has_os:
            os_samples[i] = tr.os_curve
        for k, v in tr.metrics.items():
            metric_samples[k][i] = v
        m = median_survival(t, tr.os_curve) if has_os else np.nan
        median_os_samples[i] = np.nan if m is None else m

    def _summary(arr: np.ndarray) -> dict:
        # Metrics like k_g / duration-of-response are nan when they do not apply
        # (no regrowth, no response) — summarize over the samples where they do.
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return {"median": float("nan"), "lo": float("nan"), "hi": float("nan")}
        lo, md, hi = np.percentile(finite, [ci[0], 50, ci[1]])
        return {"median": float(md), "lo": float(lo), "hi": float(hi)}

    metrics = {k: _summary(arr) for k, arr in metric_samples.items()}
    if has_os:
        metrics["median_os_weeks"] = _summary(median_os_samples)

    return Ensemble(
        record_id=record_id,
        t=t,
        n=n,
        ci=ci,
        tier=central.tier,
        warnings=central.warnings,
        tumor_size=_band(tumor_samples, ci),
        os_curve=_band(os_samples, ci) if has_os else None,
        metrics=metrics,
    )
