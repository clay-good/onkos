"""Round-trip validation — the discipline that keeps exports honest.

* analytic vs. SciPy ODE integration  -> ~1e-4 (ODE);
* SBML MathML re-parsed and evaluated -> ~1e-6 vs reference rhs (algebraic);
* NONMEM $THETA re-parsed             -> exact vs dataset values.
"""

import numpy as np
import onkos
import pytest
from onkos.export.nonmem import parse_nonmem_thetas, to_nonmem
from onkos.export.reference import eval_infix, integrate
from onkos.export.registry import get_kernel, kernel_values
from onkos.export.sbml import infix_to_mathml, mathml_to_infix, to_sbml

ODE_RECORDS = [
    "growth_laws.exponential",
    "growth_laws.logistic",
    "growth_laws.gompertz",
    "resistance.claret_2009.tgi",
    "tgi_metrics.wang_2009.biexponential",
]


def _vals(record):
    spec = get_kernel(record)
    vals = kernel_values(record)
    for inp in spec.inputs:
        if inp in ("V0", "y0"):
            vals[inp] = 80.0
        elif inp == "E":
            vals[inp] = 1.0
    return spec, vals


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_analytic_matches_ode_integration(rid):
    ds = onkos.load()
    spec, vals = _vals(ds[rid])
    t = np.linspace(0.0, 52.0, 261)
    analytic = np.asarray(spec.analytic(t, vals), dtype=float)
    ode = integrate(spec, t, vals, analytic[0])
    rel = np.max(np.abs(ode - analytic) / np.maximum(np.abs(analytic), 1e-9))
    assert rel < 1e-4, f"{rid}: analytic vs ODE rel err {rel}"


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_sbml_mathml_roundtrip(rid):
    ds = onkos.load()
    record = ds[rid]
    spec, vals = _vals(record)
    state = spec.states[0]
    infix = spec.rhs_infix[state]
    recovered = mathml_to_infix(infix_to_mathml(infix))

    rng = np.random.default_rng(0)
    for _ in range(8):
        y = float(rng.uniform(1.0, 150.0))
        tt = float(rng.uniform(0.0, 52.0))
        env = {**vals, state: y, "t": tt}
        lhs = float(eval_infix(recovered, env))
        rhs = float(spec.rhs(tt, [y], vals)[0])
        assert abs(lhs - rhs) < 1e-6, f"{rid}: MathML {lhs} != ref {rhs}"

    # the serialized SBML carries the parameter values
    xml = to_sbml(record, y0=80.0, drug_effect=1.0)
    for v in kernel_values(record).values():
        assert str(v) in xml


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_nonmem_theta_roundtrip(rid):
    ds = onkos.load()
    record = ds[rid]
    nm = to_nonmem(record, y0=80.0, drug_effect=1.0)
    thetas = parse_nonmem_thetas(nm)
    expected = list(kernel_values(record).values())
    # E and y0 inputs (when present) appended after kernel params
    assert thetas[: len(expected)] == expected


def test_clinical_use_in_every_export():
    ds = onkos.load()
    from onkos._const import CLINICAL_USE
    from onkos.export import to_pharmml, to_virtual_trial_json

    r = ds["resistance.claret_2009.tgi"]
    assert CLINICAL_USE in to_sbml(r)
    assert CLINICAL_USE in to_pharmml(r)
    assert CLINICAL_USE in to_virtual_trial_json(r)
    assert CLINICAL_USE in to_nonmem(r)
