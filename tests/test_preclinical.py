"""Phase D — preclinical translation: Simeoni model + IVIVE."""

import numpy as np
import onkos
from onkos.export.reference import effect
from onkos.export.registry import get_kernel, kernel_values


def test_simeoni_unperturbed_is_exponential_then_linear():
    ds = onkos.load()
    t = np.linspace(0, 60, 241)  # days
    g = onkos.simulate(ds, "growth_laws.simeoni_exp_linear",
                       context={"tumor_type": "ovarian_xenograft"}, drug_effect=0.0, t=t)
    w = g.tumor_size
    early_logslope = (np.log(w[10]) - np.log(w[2])) / (t[10] - t[2])
    late_slope = (w[-1] - w[-10]) / (t[-1] - t[-10])
    assert abs(early_logslope - 0.20) < 0.03   # ~ lambda0 (exponential phase)
    assert abs(late_slope - 0.70) < 0.10       # ~ lambda1 (linear phase, g/day)


def test_simeoni_is_four_state_with_sum_observable():
    ds = onkos.load()
    spec = get_kernel(ds["preclinical_translation.simeoni_2004.xenograft"])
    assert spec.n_states == 4
    assert spec.observable == "x1 + x2 + x3 + x4"
    assert spec.analytic is None  # no closed form; integrated numerically


def test_drug_shrinks_xenograft_and_higher_dose_deeper():
    ds = onkos.load()
    t = np.linspace(0, 40, 161)
    ctx = {"tumor_type": "ovarian_xenograft"}
    rid = "preclinical_translation.simeoni_2004.xenograft"
    untreated = onkos.simulate(ds, rid, context=ctx, drug_effect=0.0, t=t)
    low = onkos.simulate(ds, rid, context=ctx, drug_effect=50.0, t=t)
    high = onkos.simulate(ds, rid, context=ctx, drug_effect=200.0, t=t)
    assert untreated.tumor_size[-1] > low.tumor_size[-1] > high.tumor_size[-1]
    assert high.metrics["depth_of_response"] > low.metrics["depth_of_response"]


def test_transit_chain_delays_cell_death():
    """A pulse of damage should keep total weight rising briefly (damaged cells
    are still counted in w as they traverse the transit chain) before falling."""
    ds = onkos.load()
    t = np.linspace(0, 30, 301)
    spec = get_kernel(ds["preclinical_translation.simeoni_2004.xenograft"])
    vals = kernel_values(ds["preclinical_translation.simeoni_2004.xenograft"])
    vals["w0"] = 0.2
    vals["E"] = 80.0
    from onkos.export.reference import integrate_observable
    w = integrate_observable(spec, t, vals)
    # weight keeps growing for a short period (delay), i.e. nadir is not at t=0
    assert int(np.argmin(w)) > 0


def test_preclinical_has_no_survival_curve_and_excluded_from_clinical_compare():
    ds = onkos.load()
    tr = onkos.simulate(ds, "preclinical_translation.simeoni_2004.xenograft",
                        context={"tumor_type": "ovarian_xenograft"}, drug_effect=100.0)
    assert tr.os_curve is None
    cmp = onkos.compare(ds, purpose="tgi", context={"tumor_type": "NSCLC", "line": "first"})
    assert not any("simeoni" in rid for rid, _ in cmp.excluded)
    assert not any("simeoni" in tr.record_id for tr in cmp.included)


def test_preclinical_to_human_floors_to_D():
    ds = onkos.load()
    tr = onkos.simulate(ds, "preclinical_translation.simeoni_2004.xenograft",
                        context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=100.0)
    assert tr.tier == "D"
    assert any("outside validated" in w for w in tr.warnings)


def test_ivive_power_translation():
    ds = onkos.load()
    r = ds["preclinical_translation.ivive_potency"]
    spec, v = get_kernel(r), kernel_values(r)
    # power=1 -> linear scaling: potency = scale * IC50
    assert float(effect(spec, 20.0, v)) == v["scale"] * 20.0
    # monotone increasing in the in-vitro metric
    assert float(effect(spec, 40.0, v)) > float(effect(spec, 10.0, v))


def test_simeoni_pk_driven_time_varying():
    ds = onkos.load()
    t = np.linspace(0, 40, 161)
    ctx = {"tumor_type": "ovarian_xenograft"}
    rid = "preclinical_translation.simeoni_2004.xenograft"
    # A raw declining concentration profile drives E(t) directly (no ER record),
    # integrating the multi-state system numerically.
    conc = 200.0 * np.exp(-0.1 * t)
    pk = onkos.simulate(ds, rid, context=ctx, exposure=conc, t=t)
    untreated = onkos.simulate(ds, rid, context=ctx, drug_effect=0.0, t=t)
    assert pk.tumor_size.shape == t.shape
    assert pk.tumor_size[-1] < untreated.tumor_size[-1]  # drug exposure shrinks the tumor
