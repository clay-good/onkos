"""Round-trip validation — the discipline that keeps exports honest.

* analytic vs. SciPy ODE integration  -> ~1e-4 (single-state closed forms);
* SBML MathML re-parsed and evaluated -> ~1e-6 vs reference rhs (per state, incl.
  the multi-state Simeoni system);
* NONMEM $THETA re-parsed             -> exact vs dataset values.
"""

import numpy as np
import onkos
import pytest
from onkos.export.nonmem import parse_nonmem_thetas, to_nonmem
from onkos.export.reference import eval_infix, integrate
from onkos.export.registry import get_kernel, kernel_values
from onkos.export.sbml import infix_to_mathml, mathml_to_infix, to_sbml

# Single-state records that carry a closed-form analytic solution.
ANALYTIC_RECORDS = [
    "growth_laws.exponential",
    "growth_laws.logistic",
    "growth_laws.gompertz",
    "resistance.claret_2009.tgi",
    "tgi_metrics.wang_2009.biexponential",
    "drug_effect.norton_simon.nsclc",
]

# Every ODE record (incl. multi-state Simeoni, the IO QSP, and the two-population
# resistance model) is round-tripped.
ODE_RECORDS = ANALYTIC_RECORDS + [
    "growth_laws.simeoni_exp_linear",
    "preclinical_translation.simeoni_2004.xenograft",
    "immuno_oncology.kuznetsov_1994.tumor_immune",
    "resistance.nsclc_first_line.two_population",
    "resistance.nsclc_first_line.acquired",
]


def _vals(record):
    spec = get_kernel(record)
    vals = kernel_values(record)
    for inp in spec.inputs:
        if inp in ("V0", "y0", "w0", "T0"):
            vals[inp] = 80.0
        elif inp == "E":
            vals[inp] = 1.0
    return spec, vals


@pytest.mark.parametrize("rid", ANALYTIC_RECORDS)
def test_analytic_matches_ode_integration(rid):
    ds = onkos.load()
    spec, vals = _vals(ds[rid])
    t = np.linspace(0.0, 52.0, 261)
    analytic = np.asarray(spec.analytic(t, vals), dtype=float)
    ode = integrate(spec, t, vals, analytic[0])
    # Scale-robust relative error: floor the denominator at 0.1% of the peak so
    # trajectories that decay toward zero (Norton-Simon Gompertz collapse) stay
    # well-conditioned instead of dividing by ~0.
    floor = 1e-3 * np.max(np.abs(analytic))
    rel = np.max(np.abs(ode - analytic) / np.maximum(np.abs(analytic), floor))
    assert rel < 1e-4, f"{rid}: analytic vs ODE rel err {rel}"


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_sbml_mathml_roundtrip(rid):
    """Each state's rate law survives infix -> MathML -> infix -> eval, matching
    the hand-written reference rhs at random points (per state)."""
    ds = onkos.load()
    record = ds[rid]
    spec, vals = _vals(record)
    n = spec.n_states
    recovered = {s: mathml_to_infix(infix_to_mathml(spec.rhs_infix[s])) for s in spec.states}

    rng = np.random.default_rng(0)
    for _ in range(8):
        yvec = rng.uniform(1.0, 50.0, size=n)
        tt = float(rng.uniform(0.0, 52.0))
        env = {**vals, **{s: float(yvec[i]) for i, s in enumerate(spec.states)}, "t": tt}
        ref = spec.rhs(tt, list(yvec), vals)
        for i, s in enumerate(spec.states):
            lhs = float(eval_infix(recovered[s], env))
            assert abs(lhs - float(ref[i])) < 1e-6, f"{rid}/{s}: MathML {lhs} != ref {ref[i]}"

    # the serialized SBML carries the parameter values and one species per state
    xml = to_sbml(record, y0=80.0, drug_effect=1.0)
    for v in kernel_values(record).values():
        assert str(v) in xml
    for s in spec.states:
        assert f'species id="{s}"' in xml


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_nonmem_theta_roundtrip(rid):
    ds = onkos.load()
    record = ds[rid]
    nm = to_nonmem(record, y0=80.0, drug_effect=1.0)
    thetas = parse_nonmem_thetas(nm)
    expected = list(kernel_values(record).values())
    # E and y0 inputs (when present) appended after kernel params
    assert thetas[: len(expected)] == expected
    # one compartment per state
    assert nm.count("COMP=(") == record_n_states(record)


def record_n_states(record):
    return get_kernel(record).n_states


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_rxode2_and_pumas_param_roundtrip(rid):
    """The R and Julia exports carry the dataset parameter values verbatim."""
    from onkos.export.pumas import parse_pumas_params, to_pumas
    from onkos.export.rxode2 import parse_rxode2_params, to_rxode2

    ds = onkos.load()
    record = ds[rid]
    expected = kernel_values(record)
    rx = parse_rxode2_params(to_rxode2(record, y0=80.0, drug_effect=1.0))
    pm = parse_pumas_params(to_pumas(record, y0=80.0, drug_effect=1.0))
    for k, v in expected.items():
        assert rx[k] == v, f"{rid}: rxode2 {k}"
        assert pm[k] == v, f"{rid}: pumas {k}"


@pytest.mark.parametrize("rid", ODE_RECORDS)
def test_cross_format_parameter_consistency(rid):
    """All model formats agree on the parameter values (one source of truth)."""
    from onkos.export.nonmem import parse_nonmem_thetas
    from onkos.export.pharmml_so import parse_so_estimates, to_pharmml_so
    from onkos.export.pumas import parse_pumas_params, to_pumas
    from onkos.export.rxode2 import parse_rxode2_params, to_rxode2

    ds = onkos.load()
    record = ds[rid]
    expected = kernel_values(record)
    keys = list(expected)

    so = {k: v for k, v in parse_so_estimates(to_pharmml_so(record)).items() if k in keys}
    rx = {k: v for k, v in parse_rxode2_params(to_rxode2(record)).items() if k in keys}
    pm = {k: v for k, v in parse_pumas_params(to_pumas(record)).items() if k in keys}
    nm = parse_nonmem_thetas(to_nonmem(record))[: len(keys)]
    assert so == expected
    assert rx == expected
    assert pm == expected
    assert nm == list(expected.values())


def test_so_carries_iiv_variance_and_validation():
    from onkos.export.pharmml_so import to_pharmml_so

    ds = onkos.load()
    so = to_pharmml_so(ds["resistance.claret_2009.tgi"])
    assert "randomEffectVariance" in so          # IIV as omega = ln(1+CV^2)
    assert "externalValidation" in so            # predictive_performance recorded
    assert "OS_C_index_external" in so


def test_clinical_use_in_every_export():
    ds = onkos.load()
    from onkos._const import CLINICAL_USE
    from onkos.export import to_pharmml, to_pharmml_so, to_virtual_trial_json

    for rid in ["resistance.claret_2009.tgi", "preclinical_translation.simeoni_2004.xenograft"]:
        r = ds[rid]
        assert CLINICAL_USE in to_sbml(r)
        assert CLINICAL_USE in to_pharmml(r)
        assert CLINICAL_USE in to_pharmml_so(r)
        assert CLINICAL_USE in to_virtual_trial_json(r)
        assert CLINICAL_USE in to_nonmem(r)


def test_omex_bundles_pharmml_so():
    import zipfile

    from onkos.export.combine import build_omex

    ds = onkos.load()
    out = build_omex(ds["resistance.claret_2009.tgi"], "/tmp/onkos_so_test.omex")
    names = zipfile.ZipFile(out).namelist()
    assert any(n.endswith(".so.xml") for n in names)
    assert any(n.endswith(".pharmml") for n in names)
