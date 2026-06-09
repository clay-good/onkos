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
    """Initial state vector: the seed input fills the first state, others 0."""
    seed = spec.init_input or next((i for i in spec.inputs if i in ("V0", "y0", "w0")), None)
    y0 = np.zeros(spec.n_states, dtype=float)
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
# Claret 2009 clinical TGI: dy/dt = kL*y - kD*E*exp(-lam*t)*y
# ----------------------------------------------------------------------------
def _claret_analytic(t, v):
    kL, kD, lam, E, y0 = v["kL"], v["kD"], v["lam"], v["E"], v["y0"]
    return y0 * np.exp(kL * t - (kD * E / lam) * (1.0 - np.exp(-lam * t)))


def _claret_rhs(t, y, v):
    return [(v["kL"] - v["kD"] * v["E"] * np.exp(-v["lam"] * t)) * y[0]]


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
}
