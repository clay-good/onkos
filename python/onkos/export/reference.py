"""Pure-NumPy/SciPy reference kernels — the single computational ground truth.

Every model binds to one kernel here. Each ODE kernel provides three
*independent* expressions of the same dynamics:

1. ``analytic(t, vals)``      — closed-form solution (hand-written);
2. ``rhs(t, y, vals)``        — the ODE right-hand side (hand-written);
3. ``rhs_infix[state]``       — a symbolic string used to generate SBML/NONMEM.

The validation discipline (see ``tests/test_roundtrip.py``):

* analytic vs. SciPy integration of ``rhs``  -> agreement to ~1e-4 (ODE);
* exported SBML/NONMEM re-parsed and the rate law evaluated against ``rhs``
  -> agreement to ~1e-6 (algebraic). An export bug cannot ship silently.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

# Namespace allowed inside reference infix expressions.
_FUNCS = {"exp": np.exp, "log": np.log, "ln": np.log, "sqrt": np.sqrt, "pow": pow}


def eval_infix(expr: str, env: dict[str, float]):
    """Evaluate an internal infix expression in a restricted namespace.

    Only used on strings defined in this module — never on external input."""
    return eval(expr, {"__builtins__": {}}, {**_FUNCS, **env})  # noqa: S307


@dataclass
class KernelSpec:
    name: str
    kind: str  # "ode" | "survival" | "exposure_response"
    states: list[str]
    params: list[str]  # kernel-internal parameter names (infix-safe)
    record_symbols: list[str]  # corresponding record symbols (same order)
    inputs: list[str]  # simulation-supplied inputs (initial conditions, effect, covariate)
    analytic: Callable | None = None  # closed form; absent for multi-state systems
    rhs: Callable | None = None
    rhs_infix: dict[str, str] = field(default_factory=dict)
    # Infix mapping the state vector to the observed quantity (e.g. total tumor
    # weight = sum of compartments). None -> the first state is observed.
    observable: str | None = None
    # Which input seeds the first state's initial condition (V0/y0/w0).
    init_input: str | None = None
    # For multi-state systems whose non-first states start nonzero: the value
    # name (input or parameter) seeding each state, by position.
    init_inputs: list[str] = field(default_factory=list)
    # Survival kernels with a nonparametric (tabulated) baseline hazard (Cox):
    # simulate injects the record's ``structure.baseline_survival`` table.
    uses_baseline: bool = False

    def map_values(self, record_values: dict[str, float]) -> dict[str, float]:
        """Translate record symbols to kernel-internal names (e.g. lambda->lam)."""
        out = {}
        for kname, rsym in zip(self.params, self.record_symbols):
            if rsym in record_values:
                out[kname] = record_values[rsym]
        return out

    @property
    def n_states(self) -> int:
        return len(self.states)


def init_vector(spec: KernelSpec, vals: dict[str, float]) -> np.ndarray:
    """Initial state vector.

    With ``init_inputs`` (multi-state), each state is seeded by the named value
    (input or parameter). Otherwise the single seed input fills the first state
    and the rest start at zero."""
    y0 = np.zeros(spec.n_states, dtype=float)
    if spec.init_inputs:
        for i, name in enumerate(spec.init_inputs):
            y0[i] = float(vals.get(name, 0.0))
        return y0
    seed = spec.init_input or next((i for i in spec.inputs if i in ("V0", "y0", "w0")), None)
    if seed is not None and seed in vals:
        y0[0] = float(vals[seed])
    return y0


def integrate(spec: KernelSpec, t: np.ndarray, vals: dict[str, float], y0) -> np.ndarray:
    """Integrate a single-state ODE rhs over ``t`` (for ODE-vs-analytic checks)."""
    y0 = np.atleast_1d(np.asarray(y0, dtype=float))
    sol = solve_ivp(
        lambda tt, yy: spec.rhs(tt, yy, vals),
        (float(t[0]), float(t[-1])),
        y0,
        t_eval=t,
        rtol=1e-8,
        atol=1e-10,
        method="LSODA",
    )
    return sol.y[0]


def observe(spec: KernelSpec, states: np.ndarray) -> np.ndarray:
    """Map the integrated state array (n_states x len(t)) to the observable."""
    if spec.observable is None:
        return states[0]
    env = {name: states[i] for i, name in enumerate(spec.states)}
    return np.asarray(eval_infix(spec.observable, env), dtype=float)


def integrate_observable(
    spec: KernelSpec, t: np.ndarray, vals: dict[str, float], e_series=None
) -> np.ndarray:
    """Integrate a (possibly multi-state) ODE system and return the observable.

    ``e_series`` may be None (E taken from ``vals``), a scalar, or an array
    aligned to ``t`` (a time-varying drug effect / PK profile)."""
    y0 = init_vector(spec, vals)
    e_arr = None if e_series is None else np.atleast_1d(np.asarray(e_series, dtype=float))
    time_varying = e_arr is not None and e_arr.size == t.size and e_arr.size > 1

    def rhs(tt, yy):
        v = vals
        if e_arr is not None:
            v = dict(vals)
            v["E"] = float(np.interp(tt, t, e_arr)) if time_varying else float(e_arr.reshape(-1)[0])
        return spec.rhs(tt, yy, v)

    sol = solve_ivp(
        rhs, (float(t[0]), float(t[-1])), y0, t_eval=t, rtol=1e-8, atol=1e-10, method="LSODA"
    )
    return observe(spec, sol.y)


# ----------------------------------------------------------------------------
# Growth laws
# ----------------------------------------------------------------------------
def _exp_analytic(t, v):
    return v["V0"] * np.exp(v["kg"] * t)


def _exp_rhs(t, y, v):
    return [v["kg"] * y[0]]


def _logistic_analytic(t, v):
    V0, Vmax, kg = v["V0"], v["Vmax"], v["kg"]
    return Vmax / (1.0 + ((Vmax - V0) / V0) * np.exp(-kg * t))


def _logistic_rhs(t, y, v):
    return [v["kg"] * y[0] * (1.0 - y[0] / v["Vmax"])]


def _gompertz_analytic(t, v):
    V0, Vmax, kg = v["V0"], v["Vmax"], v["kg"]
    return Vmax * np.exp(np.log(V0 / Vmax) * np.exp(-kg * t))


def _gompertz_rhs(t, y, v):
    return [v["kg"] * y[0] * np.log(v["Vmax"] / y[0])]


# ----------------------------------------------------------------------------
# Von Bertalanffy (surface-area-limited) growth: dV/dt = a*V^(2/3) - b*V.
# Proliferation scales with the tumor SURFACE (V^2/3, nutrient/oxygen access at
# the rim) while loss scales with VOLUME (V) — the classic ontogenetic-growth
# law. Carrying capacity V_inf = (a/b)^3 (where dV/dt = 0). The substitution
# u = V^(1/3) linearizes it (du/dt = (a - b*u)/3), giving the closed form below.
# ----------------------------------------------------------------------------
def _von_bertalanffy_analytic(t, v):
    V0, a, b = v["V0"], v["a"], v["b"]
    c = a / b  # = V_inf^(1/3)
    return (c + (V0 ** (1.0 / 3.0) - c) * np.exp(-b * t / 3.0)) ** 3


def _von_bertalanffy_rhs(t, y, v):
    return [v["a"] * y[0] ** (2.0 / 3.0) - v["b"] * y[0]]


# ----------------------------------------------------------------------------
# Power-law (sub-exponential) growth: dV/dt = a*V^p, p < 1. Benzekry et al.
# (2014) found this the best-fitting unperturbed law across many tumor datasets:
# the specific growth rate a*V^(p-1) falls with size, so growth is sub-exponential
# (slower than ANY exponential) but unbounded — there is no carrying capacity.
# Separable: V^(1-p) is linear in t, giving the closed form below (p != 1; the
# p -> 1 limit is the dedicated exponential kernel).
# ----------------------------------------------------------------------------
def _power_law_analytic(t, v):
    V0, a, p = v["V0"], v["a"], v["p"]
    return (V0 ** (1.0 - p) + a * (1.0 - p) * t) ** (1.0 / (1.0 - p))


def _power_law_rhs(t, y, v):
    return [v["a"] * y[0] ** v["p"]]


# ----------------------------------------------------------------------------
# Claret 2009 clinical TGI: dy/dt = kL*y - kD*E*exp(-lam*t)*y
# ----------------------------------------------------------------------------
def _claret_analytic(t, v):
    kL, kD, lam, E, y0 = v["kL"], v["kD"], v["lam"], v["E"], v["y0"]
    return y0 * np.exp(kL * t - (kD * E / lam) * (1.0 - np.exp(-lam * t)))


def _claret_rhs(t, y, v):
    return [(v["kL"] - v["kD"] * v["E"] * np.exp(-v["lam"] * t)) * y[0]]


# ----------------------------------------------------------------------------
# Norton-Simon kill on a Gompertz growth law: dV/dt = (g - k*E)*V*ln(Vmax/V).
# The kill term is proportional to the GROWTH rate (not to tumor size), so a
# small (fast-growing) tumor is more chemo-sensitive than a large one near
# carrying capacity — the Norton-Simon hypothesis. Distinct from the log-kill
# mechanism (kill proportional to V) used by the Claret model.
# ----------------------------------------------------------------------------
def _norton_simon_analytic(t, v):
    g, Vmax, k, E, V0 = v["g"], v["Vmax"], v["k"], v["E"], v["V0"]
    return Vmax * np.exp(np.log(V0 / Vmax) * np.exp(-(g - k * E) * t))


def _norton_simon_rhs(t, y, v):
    vc = max(float(y[0]), 1e-12)  # clip away from the ln singularity at V=0
    return [(v["g"] - v["k"] * v["E"]) * vc * np.log(v["Vmax"] / vc)]


# ----------------------------------------------------------------------------
# Two-population (Goldie-Coldman) resistance: a tumor of a drug-SENSITIVE clone S
# and a pre-existing drug-RESISTANT clone R, observed together as V = S + R.
#   dS/dt = (kg - kd*E)*S      sensitive: net growth kg, killed at potency kd by E
#   dR/dt =  kgr*R             resistant: grows at kgr, NOT killed by the drug
# S starts at the measured baseline V0; R starts at a small pre-existing burden R0.
# The kill potency kd matches the Claret parameterization (kd*E is the per-time kill
# on sensitive cells), so a comparison against the phenomenological decay-of-effect
# (Claret) model isolates the RESISTANCE MECHANISM rather than an effect-scale
# difference. Distinct from Claret's exponential-decay-of-kill: here regrowth is the
# resistant clone outgrowing, and R0 is a biologically interpretable parameter.
# ----------------------------------------------------------------------------
def _two_pop_rhs(t, y, v):
    s, r = y
    return [(v["kg"] - v["kd"] * v["E"]) * s, v["kgr"] * r]


# ----------------------------------------------------------------------------
# Acquired (drug-induced) resistance: same sensitive/resistant two-clone system,
# but resistance is GENERATED during treatment rather than pre-existing. Under
# drug pressure E, sensitive cells convert to resistant at switching rate alpha:
#   dS/dt = (kg - kd*E)*S - alpha*E*S      sensitive: grow, killed, AND switching out
#   dR/dt =  kgr*R        + alpha*E*S      resistant: grow + the acquired influx
# The resistant clone typically starts at R0 = 0 (none pre-exists); the regrowth
# tail is the treatment-generated resistant pool. Contrast with two_population
# (R0 > 0, alpha = 0): same kg/kd/kgr, so the difference is purely the resistance
# ORIGIN — pre-existing burden vs acquired switching — a model-selection axis.
# ----------------------------------------------------------------------------
def _acquired_rhs(t, y, v):
    s, r = y
    flux = v["alpha"] * v["E"] * s
    return [(v["kg"] - v["kd"] * v["E"]) * s - flux, v["kgr"] * r + flux]


# ----------------------------------------------------------------------------
# Biexponential clinical TGI: y = y0*(exp(-ks*E*t) + exp(kg*t) - 1)
# ----------------------------------------------------------------------------
def _biexp_analytic(t, v):
    kg, ks, E, y0 = v["kg"], v["ks"], v["E"], v["y0"]
    return y0 * (np.exp(-ks * E * t) + np.exp(kg * t) - 1.0)


def _biexp_rhs(t, y, v):
    kg, ks, E, y0 = v["kg"], v["ks"], v["E"], v["y0"]
    return [y0 * (kg * np.exp(kg * t) - ks * E * np.exp(-ks * E * t))]


# ----------------------------------------------------------------------------
# Survival: Weibull proportional hazards
# S(t) = exp(-(t/scale)**shape * exp(beta*x))
# ----------------------------------------------------------------------------
def _weibull_ph_analytic(t, v):
    shape, scale, beta, x = v["weibull_shape"], v["weibull_scale"], v["beta"], v["x"]
    return np.exp(-((t / scale) ** shape) * np.exp(beta * x))


# Cox proportional hazards with a NONPARAMETRIC (tabulated) baseline survival
# S0(t): S(t | x) = S0(t) ** exp(beta * x). The baseline is interpolated from the
# record's structure.baseline_survival table; this is what distinguishes a Cox
# model from the parametric Weibull form (the baseline comes from data, not a
# closed-form distribution).
def _cox_ph_analytic(t, v):
    bt, b_s0 = v["baseline_times"], v["baseline_survival"]
    s0 = np.interp(t, bt, b_s0, left=1.0, right=float(b_s0[-1]))
    s0 = np.clip(s0, 1e-12, 1.0)
    return s0 ** np.exp(v["beta"] * v["x"])


# ----------------------------------------------------------------------------
# Exposure-response transforms: PK exposure metric C -> drug-effect magnitude E.
# These are algebraic (kind="exposure_response"); ``analytic(C, vals)`` returns E.
# They let a PK exposure (optionally piped from a Hypnos record) drive the kill
# term of a TGI model, completing the PK -> exposure -> tumor-dynamics chain.
# ----------------------------------------------------------------------------
def _er_emax(c, v):
    return v["Emax"] * c / (v["EC50"] + c)


def _er_sigmoid_emax(c, v):
    g = v["gamma"]
    return v["Emax"] * c**g / (v["EC50"] ** g + c**g)


def _er_power(c, v):
    return v["slope"] * c ** v["theta"]


def effect(spec: KernelSpec, exposure, vals: dict[str, float]):
    """Evaluate an exposure-response (or IVIVE) transform kernel."""
    if spec.kind != "exposure_response":
        raise ValueError(f"kernel '{spec.name}' is not a transform kernel")
    return spec.analytic(np.asarray(exposure, dtype=float), vals)


# ----------------------------------------------------------------------------
# Simeoni 2004 preclinical xenograft model (multi-state; no closed form).
#
# Unperturbed growth is exponential then linear:
#     g(w) = lam0 * x1 / (1 + (lam0 * w / lam1)^psi)^(1/psi)
# (psi large -> hard switch: rate ~ lam0*w while small, ~ lam1 once large).
# Drug (concentration E) damages proliferating cells, which then traverse a
# signal-distribution transit chain x2->x3->x4 before dying (delayed cell death):
#     dx1/dt = g(w) - k2*E*x1
#     dx2/dt = k2*E*x1 - k1*x2
#     dx3/dt = k1*x2   - k1*x3
#     dx4/dt = k1*x3   - k1*x4
# Observed tumor weight w = x1 + x2 + x3 + x4.
# ----------------------------------------------------------------------------
def _simeoni_growth(x1, w, v):
    return v["lam0"] * x1 / (1.0 + (v["lam0"] * w / v["lam1"]) ** v["psi"]) ** (1.0 / v["psi"])


def _simeoni_exp_linear_rhs(t, y, v):
    w = y[0]
    return [_simeoni_growth(w, w, v)]


def _simeoni_tgi_rhs(t, y, v):
    x1, x2, x3, x4 = y
    w = x1 + x2 + x3 + x4
    g = _simeoni_growth(x1, w, v)
    kill = v["k2"] * v["E"] * x1
    return [g - kill, kill - v["k1"] * x2, v["k1"] * x2 - v["k1"] * x3, v["k1"] * x3 - v["k1"] * x4]


# In-vitro -> in-vivo potency translation (power-law): in-vivo potency from an
# in-vitro metric (e.g. IC50). power=1 is linear scaling.
def _ivive_power(c, v):
    return v["scale"] * c ** v["power"]


# ----------------------------------------------------------------------------
# Immuno-oncology tumor-immune QSP (Kuznetsov 1994, nondimensional).
#
# HYPOTHESIS-TIER, NOT FOR PREDICTION. Effector cells and tumor cells interact
# predator-prey-style; the model reproduces immune control / dormancy / escape
# QUALITATIVELY. An immunotherapy effect E (e.g. checkpoint blockade) augments
# the immune-mediated kill (1+E). Effectors start at 0 (immune-naive).
#     d tumor/dt    = alpha*tumor*(1 - beta*tumor) - (1+E)*effector*tumor
#     d effector/dt = s + rho*effector*tumor/(eta+tumor) - mu*effector*tumor - delta*effector
# ----------------------------------------------------------------------------
def _io_tumor_immune_rhs(t, y, v):
    tumor, eff = y
    d_tumor = v["alpha"] * tumor * (1 - v["beta"] * tumor) - (1 + v["E"]) * eff * tumor
    d_eff = (
        v["s"]
        + v["rho"] * eff * tumor / (v["eta"] + tumor)
        - v["mu"] * eff * tumor
        - v["delta"] * eff
    )
    return [d_tumor, d_eff]


KERNELS: dict[str, KernelSpec] = {
    "growth_exponential": KernelSpec(
        name="growth_exponential",
        kind="ode",
        states=["tumor_size"],
        params=["kg"],
        record_symbols=["kg"],
        inputs=["V0"],
        analytic=_exp_analytic,
        rhs=_exp_rhs,
        rhs_infix={"tumor_size": "kg * tumor_size"},
    ),
    "growth_logistic": KernelSpec(
        name="growth_logistic",
        kind="ode",
        states=["tumor_size"],
        params=["kg", "Vmax"],
        record_symbols=["kg", "Vmax"],
        inputs=["V0"],
        analytic=_logistic_analytic,
        rhs=_logistic_rhs,
        rhs_infix={"tumor_size": "kg * tumor_size * (1 - tumor_size / Vmax)"},
    ),
    "growth_gompertz": KernelSpec(
        name="growth_gompertz",
        kind="ode",
        states=["tumor_size"],
        params=["kg", "Vmax"],
        record_symbols=["kg", "Vmax"],
        inputs=["V0"],
        analytic=_gompertz_analytic,
        rhs=_gompertz_rhs,
        rhs_infix={"tumor_size": "kg * tumor_size * ln(Vmax / tumor_size)"},
    ),
    "growth_von_bertalanffy": KernelSpec(
        name="growth_von_bertalanffy",
        kind="ode",
        states=["tumor_size"],
        params=["a", "b"],
        record_symbols=["a", "b"],
        inputs=["V0"],
        analytic=_von_bertalanffy_analytic,
        rhs=_von_bertalanffy_rhs,
        rhs_infix={"tumor_size": "a * tumor_size ** (2 / 3) - b * tumor_size"},
    ),
    "growth_power_law": KernelSpec(
        name="growth_power_law",
        kind="ode",
        states=["tumor_size"],
        params=["a", "p"],
        record_symbols=["a", "p"],
        inputs=["V0"],
        analytic=_power_law_analytic,
        rhs=_power_law_rhs,
        rhs_infix={"tumor_size": "a * tumor_size ** p"},
    ),
    "claret_tgi": KernelSpec(
        name="claret_tgi",
        kind="ode",
        states=["tumor_size"],
        params=["kL", "kD", "lam"],
        record_symbols=["kL", "kD", "lambda"],
        inputs=["y0", "E"],
        analytic=_claret_analytic,
        rhs=_claret_rhs,
        rhs_infix={"tumor_size": "kL * tumor_size - kD * E * exp(-lam * t) * tumor_size"},
    ),
    "norton_simon": KernelSpec(
        name="norton_simon",
        kind="ode",
        states=["tumor_size"],
        params=["g", "Vmax", "k"],
        record_symbols=["g", "Vmax", "k"],
        inputs=["V0", "E"],
        analytic=_norton_simon_analytic,
        rhs=_norton_simon_rhs,
        rhs_infix={"tumor_size": "(g - k * E) * tumor_size * ln(Vmax / tumor_size)"},
        init_input="V0",
    ),
    "two_population_resistance": KernelSpec(
        name="two_population_resistance",
        kind="ode",
        states=["sensitive", "resistant"],
        params=["kg", "kd", "kgr", "R0"],
        record_symbols=["kg", "kd", "kg_resistant", "R0"],
        inputs=["V0", "E"],
        rhs=_two_pop_rhs,
        rhs_infix={
            "sensitive": "(kg - kd * E) * sensitive",
            "resistant": "kgr * resistant",
        },
        observable="sensitive + resistant",
        # Sensitive clone seeded by the measured baseline V0; resistant clone by the
        # pre-existing burden R0 (a parameter) — the multi-state init pattern.
        init_inputs=["V0", "R0"],
    ),
    "acquired_resistance": KernelSpec(
        name="acquired_resistance",
        kind="ode",
        states=["sensitive", "resistant"],
        params=["kg", "kd", "kgr", "alpha", "R0"],
        record_symbols=["kg", "kd", "kg_resistant", "alpha", "R0"],
        inputs=["V0", "E"],
        rhs=_acquired_rhs,
        rhs_infix={
            "sensitive": "(kg - kd * E) * sensitive - alpha * E * sensitive",
            "resistant": "kgr * resistant + alpha * E * sensitive",
        },
        observable="sensitive + resistant",
        # Sensitive clone seeded by V0; resistant clone by R0 (≈ 0 — resistance is
        # acquired, not pre-existing — the regrowth pool is generated by switching).
        init_inputs=["V0", "R0"],
    ),
    "biexp_tgi": KernelSpec(
        name="biexp_tgi",
        kind="ode",
        states=["tumor_size"],
        params=["kg", "ks"],
        record_symbols=["kg", "ks"],
        inputs=["y0", "E"],
        analytic=_biexp_analytic,
        rhs=_biexp_rhs,
        rhs_infix={"tumor_size": "y0 * (kg * exp(kg * t) - ks * E * exp(-ks * E * t))"},
    ),
    "survival_weibull_ph": KernelSpec(
        name="survival_weibull_ph",
        kind="survival",
        states=["survival_fraction"],
        params=["weibull_shape", "weibull_scale", "beta"],
        record_symbols=["weibull_shape", "weibull_scale", "beta"],
        inputs=["x"],
        analytic=_weibull_ph_analytic,
    ),
    "survival_cox_ph": KernelSpec(
        name="survival_cox_ph",
        kind="survival",
        states=["survival_fraction"],
        params=["beta"],
        record_symbols=["beta"],
        inputs=["x"],
        analytic=_cox_ph_analytic,
        uses_baseline=True,
    ),
    "er_emax": KernelSpec(
        name="er_emax",
        kind="exposure_response",
        states=["effect"],
        params=["Emax", "EC50"],
        record_symbols=["Emax", "EC50"],
        inputs=["C"],
        analytic=_er_emax,
    ),
    "er_sigmoid_emax": KernelSpec(
        name="er_sigmoid_emax",
        kind="exposure_response",
        states=["effect"],
        params=["Emax", "EC50", "gamma"],
        record_symbols=["Emax", "EC50", "gamma"],
        inputs=["C"],
        analytic=_er_sigmoid_emax,
    ),
    "er_power": KernelSpec(
        name="er_power",
        kind="exposure_response",
        states=["effect"],
        params=["slope", "theta"],
        record_symbols=["slope", "theta"],
        inputs=["C"],
        analytic=_er_power,
    ),
    "simeoni_exp_linear": KernelSpec(
        name="simeoni_exp_linear",
        kind="ode",
        states=["tumor_size"],
        params=["lam0", "lam1", "psi"],
        record_symbols=["lambda0", "lambda1", "psi"],
        inputs=["w0"],
        rhs=_simeoni_exp_linear_rhs,
        rhs_infix={
            "tumor_size": "lam0 * tumor_size / (1 + (lam0 * tumor_size / lam1)**psi)**(1/psi)"
        },
        init_input="w0",
    ),
    "simeoni_tgi": KernelSpec(
        name="simeoni_tgi",
        kind="ode",
        states=["x1", "x2", "x3", "x4"],
        params=["lam0", "lam1", "psi", "k2", "k1"],
        record_symbols=["lambda0", "lambda1", "psi", "k2", "k1"],
        inputs=["w0", "E"],
        rhs=_simeoni_tgi_rhs,
        rhs_infix={
            "x1": "lam0 * x1 / (1 + (lam0 * (x1 + x2 + x3 + x4) / lam1)**psi)**(1/psi) - k2 * E * x1",
            "x2": "k2 * E * x1 - k1 * x2",
            "x3": "k1 * x2 - k1 * x3",
            "x4": "k1 * x3 - k1 * x4",
        },
        observable="x1 + x2 + x3 + x4",
        init_input="w0",
    ),
    "ivive_power": KernelSpec(
        name="ivive_power",
        kind="exposure_response",
        states=["in_vivo_potency"],
        params=["scale", "power"],
        record_symbols=["scale", "power"],
        inputs=["C"],
        analytic=_ivive_power,
    ),
    "io_tumor_immune": KernelSpec(
        name="io_tumor_immune",
        kind="ode",
        states=["tumor", "effector"],
        params=["alpha", "beta", "s", "rho", "eta", "mu", "delta", "eff0"],
        record_symbols=["alpha", "beta", "s", "rho", "eta", "mu", "delta", "eff0"],
        inputs=["T0", "E"],
        rhs=_io_tumor_immune_rhs,
        rhs_infix={
            "tumor": "alpha * tumor * (1 - beta * tumor) - (1 + E) * effector * tumor",
            "effector": "s + rho * effector * tumor / (eta + tumor) "
            "- mu * effector * tumor - delta * effector",
        },
        init_inputs=["T0", "eff0"],
    ),
}
