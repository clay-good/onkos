"""Model-agnostic TGI-metric extraction (the Stein/Bruno panel).

Given any simulated tumor-size trajectory, extract the derived metrics oncology
pharmacometrics reports (spec §3, §6):

- **depth of response (DpR)** — fractional shrinkage from baseline to nadir;
- **nadir** and **time to nadir**;
- **tumor shrinkage-rate constant k_s** — log-linear decline rate over the
  on-treatment shrinkage phase;
- **tumor growth-rate constant k_g** — log-linear regrowth rate over the late
  (post-nadir) phase; the strongly prognostic Stein/Bruno quantity;
- **time to growth (TTG)** — nadir time when genuine regrowth follows;
- **duration of response (DoR)** — RECIST-style time from partial response
  (≥30% shrinkage from baseline) to progression (≥20% growth from nadir);
- **week-8 relative change** — the covariate the survival link consumes.

The extractor is deliberately model-agnostic: it estimates k_g / k_s from the
trajectory the same way the Stein method estimates them from RECIST data, so the
metrics are comparable across the Claret, biexponential, Simeoni, and any future
kernels. Values that do not apply (no shrinkage, no regrowth, no response) are
returned as ``nan`` rather than fabricated.
"""

from __future__ import annotations

import numpy as np

# RECIST 1.1-style thresholds on the sum of longest diameters.
_PR_SHRINK = 0.30  # partial response: ≥30% decrease from baseline
_PD_GROWTH = 0.20  # progressive disease: ≥20% increase from the nadir

# Lower limit on relative tumor size for the integrated-burden metric: a tumor at
# 0.1% of baseline is a complete response, clinically indistinguishable from zero.
# Flooring here keeps log(v/y0) finite under eradication (v → 0) without changing the
# ranking, so the integral is a stable summary rather than a −∞-dominated one.
_BURDEN_FLOOR = 1e-3


def _loglinear_slope(t: np.ndarray, v: np.ndarray) -> float:
    """d ln(v)/dt by least squares; nan if under-determined or non-positive v."""
    mask = v > 0
    if mask.sum() < 2 or np.ptp(t[mask]) <= 0:
        return float("nan")
    return float(np.polyfit(t[mask], np.log(v[mask]), 1)[0])


def _duration_of_response(t, v, y0, nadir, nadir_i) -> float:
    pr_level = (1.0 - _PR_SHRINK) * y0
    onset = next((float(t[i]) for i in range(len(t)) if v[i] <= pr_level), None)
    if onset is None:
        return float("nan")
    pd_level = (1.0 + _PD_GROWTH) * nadir
    prog = next((float(t[i]) for i in range(nadir_i, len(t)) if v[i] >= pd_level), None)
    if prog is None:
        return float("nan")
    return prog - onset


def extract_tgi_metrics(t: np.ndarray, v: np.ndarray, y0: float) -> dict:
    """Extract the TGI-metric panel from a tumor-size trajectory ``v`` over ``t``
    with baseline ``y0``."""
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    n = len(t)
    nadir_i = int(np.argmin(v))
    nadir = float(v[nadir_i])
    t_nadir = float(t[nadir_i])
    week8 = float(np.interp(8.0, t, v))

    m = {
        "week8_tumor_size": week8,
        "week8_relative_change": (week8 - y0) / y0,
        "nadir_tumor_size": nadir,
        "time_to_nadir_weeks": t_nadir,
        "depth_of_response": (y0 - nadir) / y0,
    }

    # k_g — growth-rate constant over the late (post-nadir) regrowth phase. Using
    # the tail of the horizon isolates the growth regime from the transition; this
    # recovers the generating growth rate for both the Claret (kL) and the
    # biexponential (kg) kernels.
    kg = float("nan")
    tail_start = max(nadir_i, int(n * 0.75))
    if n - tail_start >= 2:
        s = _loglinear_slope(t[tail_start:], v[tail_start:])
        kg = s if (np.isfinite(s) and s > 0) else float("nan")
    m["tumor_growth_rate_kg"] = kg

    # k_s — shrinkage-rate constant. For a shrink-then-grow trajectory the initial
    # instantaneous log-slope s0 equals (k_g − k_s), so k_s = k_g − s0 recovers the
    # generating shrink rate (biexp ks; Claret kD·E). With no growth phase, the
    # initial decline rate itself is k_s.
    head_end = max(2, round(0.03 * n))
    s0 = _loglinear_slope(t[:head_end], v[:head_end])
    ks = float("nan")
    if np.isfinite(s0):
        ks = (kg - s0) if np.isfinite(kg) else (-s0)
        ks = ks if ks > 1e-9 else float("nan")
    m["tumor_shrinkage_rate_ks"] = ks

    regrew = nadir_i < n - 1 and v[-1] > nadir + 1e-9 and np.isfinite(kg)
    m["time_to_growth_weeks"] = t_nadir if regrew else float("nan")
    m["duration_of_response_weeks"] = _duration_of_response(t, v, y0, nadir, nadir_i)

    # Integrated tumor burden — the time-averaged log relative tumor size over the
    # observation horizon (the AUC of the log-size curve, i.e. the log geometric-mean
    # relative burden). Unlike week-8 (one early point, blind to the tail) and unlike
    # k_g (the terminal regrowth slope, blind to depth), this single number integrates
    # *both* the depth of response and the regrowth tail. It is therefore a candidate
    # survival bridge metric that re-ranks models a third way (spec: burden-AUC). The
    # value is horizon-dependent by construction (it is a cumulative-burden summary).
    if n >= 2 and np.ptp(t) > 0:
        rel = np.maximum(v / y0, _BURDEN_FLOOR)
        logb = np.log(rel)
        # Trapezoidal integral, version-agnostic (np.trapz is deprecated in numpy 2.x,
        # np.trapezoid is absent below 2.0); divide by the span to time-average.
        area = float(np.sum((logb[1:] + logb[:-1]) * 0.5 * np.diff(t)))
        m["log_burden_auc"] = area / float(np.ptp(t))
    else:
        m["log_burden_auc"] = float("nan")
    return m
