"""Population-level forward simulation.

NOT a prognostic engine and NOT a treatment optimizer. ``simulate`` produces
tumor-size and population overall-survival *trajectories* for research, model
comparison, and export validation. It never returns an individual prognosis or
ranks therapies (see spec §6, §10).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ._const import CLINICAL_USE
from .export.reference import effect as er_effect
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
    candidates = [r for r in ds if r.purpose == "survival_link"]
    for r in candidates:
        dc = r.derivation_context
        if dc and dc.tumor_type == tumor_type:
            return r
    return candidates[0] if candidates else None


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

    If an exposure-response record and an exposure metric are supplied, E is the
    ER transform of the exposure (scalar or time series), and the ER record joins
    the tier-propagation set. Otherwise E is the scalar ``drug_effect``.
    """
    if exposure_response is not None and exposure is not None:
        er = ds[exposure_response]
        er_spec = get_kernel(er)
        e = er_effect(er_spec, exposure, kernel_values(er))
        contributing.append(er)
        return e
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
        if inp in ("V0", "y0"):
            vals[inp] = y0
        elif inp == "E":
            vals[inp] = float(e_arr[0])

    if time_varying:
        tumor = _integrate_timevarying(spec, t, vals, y0, e_arr)
    else:
        tumor = np.asarray(spec.analytic(t, vals), dtype=float)
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


def _integrate_timevarying(spec, t: np.ndarray, vals: dict, y0: float, e_arr: np.ndarray):
    """Integrate a TGI ODE with a time-varying drug effect E(t) (PK-driven)."""
    if spec.rhs is None:
        raise ValueError(f"kernel '{spec.name}' has no ODE rhs for time-varying exposure")

    def rhs_t(tt, yy):
        v = dict(vals)
        v["E"] = float(np.interp(tt, t, e_arr))
        return spec.rhs(tt, yy, v)

    sol = solve_ivp(
        rhs_t, (float(t[0]), float(t[-1])), [float(y0)], t_eval=t, rtol=1e-8, atol=1e-10,
        method="LSODA",
    )
    return sol.y[0]


def _baseline_record(ds: Dataset, tumor_type, line) -> Record | None:
    for r in ds:
        if r.kind != "context_baseline":
            continue
        dc = r.derivation_context
        if dc and dc.tumor_type == tumor_type and (line is None or dc.line_of_therapy == line):
            return r
    return None
