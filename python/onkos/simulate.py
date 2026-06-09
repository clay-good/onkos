"""Population-level forward simulation.

NOT a prognostic engine and NOT a treatment optimizer. ``simulate`` produces
tumor-size and population overall-survival *trajectories* for research, model
comparison, and export validation. It never returns an individual prognosis or
ranks therapies (see spec §6, §10).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .export.reference import effect as er_effect
from .export.reference import integrate_observable
from .export.registry import get_kernel, kernel_values
from .load import Dataset
from .models import Record
from .tiers import propagate

__all__ = ["CLINICAL_USE", "Trajectory", "simulate", "median_survival"]


@dataclass
class Trajectory:
    record_id: str
    t: np.ndarray
    tumor_size: np.ndarray
    tier: str
    warnings: list[str] = field(default_factory=list)
    os_curve: np.ndarray | None = None
    metrics: dict[str, float] = field(default_factory=dict)
    clinical_use: str = CLINICAL_USE

    @property
    def median_os(self) -> float | None:
        return median_survival(self.t, self.os_curve) if self.os_curve is not None else None


def median_survival(t: np.ndarray, s: np.ndarray) -> float | None:
    """First time the survival fraction crosses 0.5 (linear interpolation)."""
    s = np.asarray(s)
    below = np.where(s <= 0.5)[0]
    if len(below) == 0:
        return None
    i = below[0]
    if i == 0:
        return float(t[0])
    t0, t1, s0, s1 = t[i - 1], t[i], s[i - 1], s[i]
    if s1 == s0:
        return float(t1)
    return float(t0 + (0.5 - s0) * (t1 - t0) / (s1 - s0))


def _baseline_y0(ds: Dataset, tumor_type: str | None, line: str | None) -> float:
    for r in ds:
        if r.kind != "context_baseline":
            continue
        dc = r.derivation_context
        if dc and dc.tumor_type == tumor_type and (line is None or dc.line_of_therapy == line):
            if "baseline_tumor_size" in r:
                return float(r["baseline_tumor_size"].central)
    return 100.0


def _find_survival_link(ds: Dataset, tumor_type: str | None) -> Record | None:
    """Return the survival link whose context matches this tumor type.

    No fallback: an unmatched context (including any preclinical context) gets no
    OS curve rather than a survival model from an unrelated tumor type."""
    for r in ds:
        if r.purpose != "survival_link":
            continue
        dc = r.derivation_context
        if dc and dc.tumor_type == tumor_type:
            return r
    return None


def _tumor_metrics(t: np.ndarray, y: np.ndarray, y0: float) -> dict[str, float]:
    week8 = float(np.interp(8.0, t, y))
    nadir = float(np.min(y))
    nadir_t = float(t[int(np.argmin(y))])
    return {
        "week8_tumor_size": week8,
        "week8_relative_change": (week8 - y0) / y0,
        "nadir_tumor_size": nadir,
        "time_to_nadir_weeks": nadir_t,
        "depth_of_response": (y0 - nadir) / y0,
    }


def _resolve_effect(
    ds: Dataset,
    *,
    drug_effect: float | None,
    exposure,
    exposure_response: str | None,
    contributing: list[Record],
):
    """Determine the drug-effect magnitude E driving the kill term.

    Resolution order:
    * ``exposure`` + ``exposure_response`` -> E is the ER transform of the
      exposure (scalar or time series); the ER record joins tier propagation.
    * ``exposure`` alone -> the exposure drives E directly (identity), e.g. a
      concentration profile feeding a Simeoni-style ``k2*E*x1`` kill term.
    * otherwise -> the scalar ``drug_effect`` (default 1.0).
    """
    if exposure is not None:
        if exposure_response is not None:
            er = ds[exposure_response]
            contributing.append(er)
            return er_effect(get_kernel(er), exposure, kernel_values(er))
        return exposure
    return float(drug_effect if drug_effect is not None else 1.0)


def simulate(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    drug_effect: float | None = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> Trajectory:
    """Forward-simulate a TGI (or growth) record and, where a survival link is
    available, the resulting population OS curve.

    The drug effect E may be given directly (``drug_effect``) or derived from a
    PK exposure metric through an exposure-response record (``exposure`` +
    ``exposure_response``). A time-varying ``exposure`` (array aligned to ``t``,
    e.g. a Hypnos PK profile) yields a time-varying E(t) and the tumor ODE is
    integrated numerically; a scalar exposure uses the fast closed form.
    """
    context = context or {}
    tumor_type = context.get("tumor_type")
    line = context.get("line") or context.get("line_of_therapy")
    if t is None:
        t = np.linspace(0.0, 104.0, 209)  # two years, weekly-ish
    t = np.asarray(t, dtype=float)

    record = ds[record_id]
    spec = get_kernel(record)
    y0 = float(context.get("y0", _baseline_y0(ds, tumor_type, line)))

    contributing: list[Record] = [record]
    e_value = _resolve_effect(
        ds,
        drug_effect=drug_effect,
        exposure=exposure,
        exposure_response=exposure_response,
        contributing=contributing,
    )
    e_arr = np.atleast_1d(np.asarray(e_value, dtype=float))
    time_varying = e_arr.size == t.size and e_arr.size > 1

    vals = kernel_values(record)
    for inp in spec.inputs:
        if inp in ("V0", "y0", "w0"):
            vals[inp] = y0
        elif inp == "E":
            vals[inp] = float(e_arr[0])

    # Closed form for the single-state, constant-effect case; numerical
    # integration for multi-state systems (Simeoni) or a time-varying effect.
    if spec.analytic is not None and spec.n_states == 1 and not time_varying:
        tumor = np.asarray(spec.analytic(t, vals), dtype=float)
    else:
        tumor = integrate_observable(spec, t, vals, e_series=e_arr if time_varying else None)
    metrics = _tumor_metrics(t, tumor, y0)

    baseline = _baseline_record(ds, tumor_type, line)
    if baseline is not None:
        contributing.append(baseline)

    os_curve = None
    link = None
    if record.purpose in ("tgi", "metric"):
        link = ds[survival_link] if survival_link else _find_survival_link(ds, tumor_type)
    if link is not None:
        link_spec = get_kernel(link)
        link_vals = kernel_values(link)
        link_vals["x"] = metrics["week8_relative_change"]
        os_curve = np.asarray(link_spec.analytic(t, link_vals), dtype=float)
        contributing.append(link)

    prop = propagate(contributing, tumor_type=tumor_type, line=line)
    return Trajectory(
        record_id=record_id,
        t=t,
        tumor_size=tumor,
        tier=prop.tier,
        warnings=prop.warnings,
        os_curve=os_curve,
        metrics=metrics,
    )


def _baseline_record(ds: Dataset, tumor_type, line) -> Record | None:
    for r in ds:
        if r.kind != "context_baseline":
            continue
        dc = r.derivation_context
        if dc and dc.tumor_type == tumor_type and (line is None or dc.line_of_therapy == line):
            return r
    return None
