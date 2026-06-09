"""PK bridge — turn a dose/regimen (or an external Hypnos PK record) into the
exposure metric the exposure-response kernels expect.

Onkos consumes exposure; it does not model PK (that is Hypnos's job). This module
is a deliberately small, *illustrative* bridge so the spec's headline composability
claim — a full **PK → exposure → tumor-dynamics → survival** chain — can be run
self-contained:

- :func:`steady_state_metrics` returns the standard exposure metrics (C_avg, C_max,
  C_min, AUC over a dosing interval) for a one-compartment oral regimen. The
  cornerstone relation ``C_avg = F·Dose/(CL·tau)`` is what links dose to the
  exposure that drives an Emax/power ER record.
- :func:`concentration_profile` builds the multiple-dose concentration-time curve
  (Bateman superposition) for plotting or as a time-varying exposure drive.
- :func:`from_profile` ingests an external concentration-time profile (e.g. emitted
  by a Hypnos PK model) and resamples it onto the simulation grid.

For real PK, fit/simulate with Hypnos and feed the resulting profile or C_avg in;
the generators here are illustrative and clearly labelled as such.
"""

from __future__ import annotations

import numpy as np

__all__ = ["steady_state_metrics", "concentration_profile", "from_profile"]


def _single_dose(t: np.ndarray, dose: float, ka: float, ke: float, v: float, f: float) -> np.ndarray:
    """One-compartment oral single-dose concentration (Bateman function), t>=0."""
    t = np.asarray(t, dtype=float)
    c = np.zeros_like(t)
    pos = t >= 0
    if abs(ka - ke) < 1e-12:  # flip-flop / equal-rate limit
        c[pos] = (f * dose * ka / v) * t[pos] * np.exp(-ke * t[pos])
    else:
        c[pos] = (f * dose * ka) / (v * (ka - ke)) * (
            np.exp(-ke * t[pos]) - np.exp(-ka * t[pos])
        )
    return c


def concentration_profile(
    dose: float,
    tau: float,
    n_doses: int,
    *,
    ka: float,
    ke: float,
    v: float = 1.0,
    f: float = 1.0,
    t: np.ndarray,
) -> np.ndarray:
    """Multiple-dose oral concentration over ``t`` (superposition of Bateman curves).

    ``tau`` is the dosing interval and ``n_doses`` the number of doses given at
    0, tau, 2·tau, …"""
    t = np.asarray(t, dtype=float)
    total = np.zeros_like(t)
    for i in range(n_doses):
        total += _single_dose(t - i * tau, dose, ka, ke, v, f)
    return total


def steady_state_metrics(
    dose: float, tau: float, *, ka: float, ke: float, v: float = 1.0, f: float = 1.0
) -> dict:
    """Steady-state exposure metrics for a one-compartment oral regimen.

    Returns ``c_avg``, ``c_max``, ``c_min``, ``auc_tau`` and the ``accumulation``
    ratio. ``c_avg`` and ``auc_tau`` are closed-form (``C_avg = F·Dose/(CL·tau)``,
    ``CL = ke·V``); ``c_max``/``c_min`` are read from a converged superposition."""
    cl = ke * v
    c_avg = f * dose / (cl * tau)
    auc_tau = f * dose / cl  # AUC over one interval at steady state == single-dose AUC_inf

    # Converge the superposition, then read the last interval.
    n = 40
    grid = np.linspace((n - 1) * tau, n * tau, 400)
    ss = concentration_profile(dose, tau, n, ka=ka, ke=ke, v=v, f=f, t=grid)
    c_max, c_min = float(np.max(ss)), float(np.min(ss))
    single_cmax = float(np.max(_single_dose(np.linspace(0, tau, 400), dose, ka, ke, v, f)))
    accumulation = c_max / single_cmax if single_cmax > 0 else float("nan")
    return {
        "c_avg": float(c_avg),
        "c_max": c_max,
        "c_min": c_min,
        "auc_tau": float(auc_tau),
        "accumulation": accumulation,
    }


def from_profile(times, concentrations, t: np.ndarray) -> np.ndarray:
    """Resample an external concentration-time profile (e.g. a Hypnos PK record)
    onto the simulation grid ``t`` by linear interpolation.

    Outside the supplied range the profile is held at its endpoints."""
    times = np.asarray(times, dtype=float)
    concentrations = np.asarray(concentrations, dtype=float)
    order = np.argsort(times)
    return np.interp(np.asarray(t, dtype=float), times[order], concentrations[order])
