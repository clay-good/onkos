"""Scientific landmark validation — a second, independent validation axis.

The export round-trip tests (``test_reference.py``) prove the *exported* form
agrees with the reference kernel. They do **not** prove the kernel reproduces
the *scientifically meaningful behaviour* of the published model it claims to
implement: a kernel can be internally self-consistent yet wrong.

These tests close that gap. Each reference kernel is checked against the
characteristic, analytically-derivable **landmark** of its published model — a
quantitative property a correct implementation must exhibit (an inflection
point, a static-tumor concentration, a stationary growth condition, a survival
median). The landmarks and their sources are catalogued in
``docs/validation-landmarks.md``.

This is the honest form of the spec §9 directive to compare kernel output
against "published example simulations": the landmark *is* the published
property, derived from the model's own equations — no digitized data is
fabricated.
"""

import numpy as np
from onkos.export.reference import KERNELS


def _slope(fn, t, h=1e-6):
    """Central finite difference of fn at t."""
    return (fn(t + h) - fn(t - h)) / (2 * h)


# --- growth laws -----------------------------------------------------------


def test_exponential_doubling_time():
    """Exponential growth doubles every ln2/kg (definitional doubling time)."""
    s = KERNELS["growth_exponential"]
    v = {"kg": 0.05, "V0": 30.0}
    td = np.log(2) / v["kg"]
    assert np.isclose(s.analytic(td, v), 2 * v["V0"], rtol=1e-9)


def test_logistic_inflection_at_half_carrying_capacity():
    """Logistic growth rate is maximal at V = Vmax/2 (the inflection point)."""
    s = KERNELS["growth_logistic"]
    v = {"kg": 0.1, "Vmax": 200.0}

    def rate(V):
        return s.rhs(0.0, [V], v)[0]

    peak = v["Vmax"] / 2
    assert rate(peak) > rate(peak * 0.8)
    assert rate(peak) > rate(peak * 1.2)


def test_gompertz_inflection_at_vmax_over_e():
    """Gompertz growth rate is maximal at V = Vmax/e (the published inflection)."""
    s = KERNELS["growth_gompertz"]
    v = {"kg": 0.08, "Vmax": 200.0}

    def rate(V):
        return s.rhs(0.0, [V], v)[0]

    peak = v["Vmax"] / np.e
    grid = np.linspace(1.0, v["Vmax"] - 1.0, 4000)
    argmax = grid[int(np.argmax([rate(V) for V in grid]))]
    assert np.isclose(argmax, peak, rtol=2e-3)


def test_von_bertalanffy_carrying_capacity_and_inflection():
    """Von Bertalanffy growth dV/dt = a*V^(2/3) - b*V is stationary at the carrying
    capacity V_inf = (a/b)^3, and (being surface-limited) its growth rate peaks at the
    published inflection V_inf*(2/3)^3 = (2a/3b)^3 — strictly below Vmax/2, unlike logistic."""
    s = KERNELS["growth_von_bertalanffy"]
    a, b = 0.2924, 0.05
    v = {"a": a, "b": b}

    def rate(V):
        return s.rhs(0.0, [V], v)[0]

    v_inf = (a / b) ** 3
    assert np.isclose(v_inf, 200.0, rtol=1e-3)  # the record's illustrative carrying capacity
    assert np.isclose(rate(v_inf), 0.0, atol=1e-9)  # stationary at carrying capacity (exact)

    # the absolute growth rate peaks at (2a/3b)^3 (where d/dV[a V^2/3 - bV] = 0)
    peak = (2.0 * a / (3.0 * b)) ** 3
    grid = np.linspace(1.0, v_inf - 1.0, 4000)
    argmax = grid[int(np.argmax([rate(V) for V in grid]))]
    assert np.isclose(argmax, peak, rtol=3e-3)
    assert peak < v_inf / 2  # surface-limited: inflection below the logistic half-capacity


def test_von_bertalanffy_is_sub_exponential():
    """The specific growth rate (1/V)dV/dt = a*V^(-1/3) - b falls monotonically with size —
    growth is sub-exponential from the first cell, the defining surface-limited behavior."""
    s = KERNELS["growth_von_bertalanffy"]
    v = {"a": 0.2924, "b": 0.05}

    def specific_rate(V):
        return s.rhs(0.0, [V], v)[0] / V

    sizes = [5.0, 25.0, 75.0, 150.0]
    rates = [specific_rate(V) for V in sizes]
    assert all(rates[i] > rates[i + 1] for i in range(len(rates) - 1))


# --- drug-effect / TGI -----------------------------------------------------


def test_claret_initial_log_slope():
    """Claret: the instantaneous log-growth at t=0 is kL - kD*E (perturbed slope)."""
    s = KERNELS["claret_tgi"]
    v = {"kL": 0.05, "kD": 0.12, "lam": 0.20, "E": 1.0, "y0": 50.0}

    def log_y(t):
        return np.log(s.analytic(t, v))

    assert np.isclose(_slope(log_y, 1e-3), v["kL"] - v["kD"] * v["E"], atol=1e-3)


def test_norton_simon_stationary_at_e_equals_g_over_k():
    """Norton-Simon is stationary (dV/dt = 0 at every V) when E = g/k."""
    s = KERNELS["norton_simon"]
    g, k, Vmax = 0.10, 0.08, 200.0
    v = {"g": g, "k": k, "Vmax": Vmax, "E": g / k}
    for V in (20.0, 80.0, 160.0):
        assert abs(s.rhs(0.0, [V], v)[0]) < 1e-9


def test_biexp_nadir_time():
    """Bi-exponential nadir is at t* = ln(ks*E/kg)/(kg+ks*E) — a minimum, slope 0."""
    s = KERNELS["biexp_tgi"]
    v = {"kg": 0.01, "ks": 0.20, "E": 1.0, "y0": 60.0}
    tstar = np.log(v["ks"] * v["E"] / v["kg"]) / (v["kg"] + v["ks"] * v["E"])
    assert abs(_slope(lambda t: s.analytic(t, v), tstar)) < 1e-6
    assert s.analytic(tstar, v) < s.analytic(tstar * 0.5, v)
    assert s.analytic(tstar, v) < s.analytic(tstar * 2.0, v)


def test_simeoni_exp_and_linear_phase_slopes():
    """Simeoni unperturbed growth: fractional rate -> lam0 when small, absolute
    rate -> lam1 when large (the exp->linear transition the model is named for)."""
    s = KERNELS["simeoni_exp_linear"]
    v = {"lam0": 0.20, "lam1": 0.70, "psi": 20.0}
    small = 1e-3  # w << lam1/lam0 -> exponential phase
    assert np.isclose(s.rhs(0.0, [small], v)[0] / small, v["lam0"], rtol=1e-3)
    large = 1e3  # w >> lam1/lam0 -> linear phase
    assert np.isclose(s.rhs(0.0, [large], v)[0], v["lam1"], rtol=1e-3)


def test_simeoni_tgi_tumor_static_concentration():
    """Simeoni TGI: the proliferating compartment is stationary at the
    tumor-static concentration c* = lam0/k2 (growth exactly balances kill)."""
    s = KERNELS["simeoni_tgi"]
    lam0, k2 = 0.20, 0.012
    v = {"lam0": lam0, "lam1": 0.70, "psi": 20.0, "k2": k2, "k1": 0.50, "E": lam0 / k2}
    w0 = 0.05  # small tumor -> exponential phase, growth ~ lam0*x1
    dx1 = s.rhs(0.0, [w0, 0.0, 0.0, 0.0], v)[0]
    assert abs(dx1) < 1e-4


# --- survival --------------------------------------------------------------


def test_weibull_median_survival():
    """Weibull-PH median (x=0) is at t = scale*(ln2)^(1/shape) where S = 0.5."""
    s = KERNELS["survival_weibull_ph"]
    v = {"weibull_shape": 1.3, "weibull_scale": 80.0, "beta": 0.6, "x": 0.0}
    t_med = v["weibull_scale"] * np.log(2) ** (1.0 / v["weibull_shape"])
    assert np.isclose(s.analytic(t_med, v), 0.5, rtol=1e-9)


def test_weibull_proportional_hazards():
    """A unit covariate scales the cumulative hazard by exp(beta): S_x = S_0^exp(beta)."""
    s = KERNELS["survival_weibull_ph"]
    base = {"weibull_shape": 1.3, "weibull_scale": 80.0, "beta": 0.6}
    t = 50.0
    s0 = s.analytic(t, {**base, "x": 0.0})
    s1 = s.analytic(t, {**base, "x": 1.0})
    assert np.isclose(s1, s0 ** np.exp(base["beta"]), rtol=1e-9)


def test_cox_median_tracks_baseline():
    """Cox-PH at x=0 reproduces the tabulated baseline; its median is the
    baseline median (S0 = 0.5)."""
    s = KERNELS["survival_cox_ph"]
    times = [0.0, 10.0, 20.0, 30.0, 40.0]
    surv = [1.0, 0.8, 0.5, 0.3, 0.1]
    v = {"baseline_times": times, "baseline_survival": surv, "beta": 0.5, "x": 0.0}
    assert np.isclose(s.analytic(20.0, v), 0.5, rtol=1e-9)


# --- exposure-response -----------------------------------------------------


def test_emax_half_maximal_at_ec50():
    """Emax: effect is exactly Emax/2 at C = EC50 (definition of EC50)."""
    s = KERNELS["er_emax"]
    v = {"Emax": 0.9, "EC50": 25.0}
    assert np.isclose(s.analytic(v["EC50"], v), v["Emax"] / 2, rtol=1e-12)


def test_sigmoid_emax_half_maximal_at_ec50_any_hill():
    """Sigmoid Emax: still Emax/2 at C = EC50 regardless of the Hill coefficient."""
    s = KERNELS["er_sigmoid_emax"]
    for gamma in (0.7, 1.0, 2.5):
        v = {"Emax": 0.9, "EC50": 25.0, "gamma": gamma}
        assert np.isclose(s.analytic(v["EC50"], v), v["Emax"] / 2, rtol=1e-12)


def test_power_models_are_scale_free():
    """Power laws (er_power, ivive_power): doubling C scales effect by 2^exponent."""
    for name, exp_key in (("er_power", "theta"), ("ivive_power", "power")):
        s = KERNELS[name]
        coeff = "slope" if name == "er_power" else "scale"
        v = {coeff: 0.4, exp_key: 0.75}
        c = 10.0
        assert np.isclose(s.analytic(2 * c, v) / s.analytic(c, v), 2 ** v[exp_key], rtol=1e-12)


# --- immuno-oncology (hypothesis tier) -------------------------------------


def test_io_immune_homeostasis_without_tumor():
    """IO tumor-immune: with no tumor, effectors relax to s/delta (homeostasis)."""
    s = KERNELS["io_tumor_immune"]
    v = {"alpha": 0.1, "beta": 0.01, "s": 0.2, "rho": 0.3, "eta": 1.0,
         "mu": 0.05, "delta": 0.1, "E": 0.0}
    eff_star = v["s"] / v["delta"]
    d_eff = s.rhs(0.0, [0.0, eff_star], v)[1]
    assert abs(d_eff) < 1e-12
