"""Two-population (Goldie-Coldman) mechanistic resistance: a drug-sensitive clone
plus a pre-existing resistant clone, observed as one tumor.

Landmark suite (research spec mechanistic-resistance §5): the characteristic,
analytically-derivable properties of the published two-population model the kernel
implements, plus integration through the divergence view (the resistance mechanism as
a model-selection axis).
"""

import numpy as np
import onkos
from onkos.export.reference import KERNELS
from onkos.export.registry import get_kernel, kernel_values

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
RID = "resistance.nsclc_first_line.two_population"


def _kvals(ds, E=1.0):
    v = kernel_values(ds[RID])
    v["E"] = E
    return v


def test_kernel_registered_and_two_state():
    ds = onkos.load()
    assert "two_population_resistance" in KERNELS
    spec = get_kernel(ds[RID])
    assert spec.states == ["sensitive", "resistant"]
    assert spec.n_states == 2
    assert spec.observable == "sensitive + resistant"


def test_closed_form_matches_integration():
    """V(t) = V0·e^{(kg−kd·E)t} + R0·e^{kgr·t} (sensitive seeded by the baseline,
    resistant by R0)."""
    ds = onkos.load()
    t = np.linspace(0.0, 156.0, 313)
    tr = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.0, t=t)
    v = _kvals(ds)
    v0 = tr.tumor_size[0] - v["R0"]  # observed(0) = sensitive(0) + R0
    closed = v0 * np.exp((v["kg"] - v["kd"] * v["E"]) * t) + v["R0"] * np.exp(v["kgr"] * t)
    assert np.allclose(tr.tumor_size, closed, rtol=1e-4)


def test_no_resistance_reduces_to_sensitive_exponential_and_eradicates():
    """R0=0 with effective kill kd·E > kg ⇒ pure sensitive exponential, V→0 monotone."""
    ds = onkos.load()
    t = np.linspace(0.0, 104.0, 209)
    tr = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.0, t=t,
                        param_overrides={"R0": 0.0})
    assert np.all(np.diff(tr.tumor_size) < 1e-9)          # monotone decreasing
    assert tr.tumor_size[-1] < tr.tumor_size[0] * 1e-3    # eradicated


def test_resistant_clone_sets_the_late_time_growth_rate():
    """With R0>0 and kd·E>kg, the late-time log-slope of V tends to kgr."""
    ds = onkos.load()
    t = np.linspace(0.0, 312.0, 625)
    tr = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.0, t=t)
    slope = np.polyfit(t[-40:], np.log(tr.tumor_size[-40:]), 1)[0]
    assert np.isclose(slope, _kvals(ds)["kgr"], rtol=1e-3)


def test_nadir_then_regrowth():
    """A drug that kills the sensitive clone gives an interior nadir, then regrowth
    driven by the resistant clone."""
    ds = onkos.load()
    t = np.linspace(0.0, 208.0, 417)
    v = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.0, t=t).tumor_size
    i = int(np.argmin(v))
    assert 0 < i < len(v) - 1               # the minimum is interior (a real nadir)
    assert v[-1] > v[i] * 2                 # genuine regrowth after the nadir


def test_resistant_fraction_rises_monotonically_under_treatment():
    """R/V increases monotonically from R0/V0 under an effective kill."""
    ds = onkos.load()
    t = np.linspace(0.0, 156.0, 313)
    v = _kvals(ds)
    v0 = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.0, t=t).tumor_size[0] - v["R0"]
    s = v0 * np.exp((v["kg"] - v["kd"] * v["E"]) * t)
    r = v["R0"] * np.exp(v["kgr"] * t)
    frac = r / (s + r)
    assert np.all(np.diff(frac) >= -1e-12)   # monotone non-decreasing (saturates at 1.0)
    assert frac[-1] > frac[0]                 # and genuinely rises from R0/V0


def test_untreated_both_clones_grow():
    """E=0 ⇒ V(t) = V0·e^{kg t} + R0·e^{kgr t}; the tumor only grows."""
    ds = onkos.load()
    t = np.linspace(0.0, 104.0, 209)
    tr = onkos.simulate(ds, RID, context=NSCLC, drug_effect=0.0, t=t)
    assert np.all(np.diff(tr.tumor_size) > 0)


def test_appears_in_nsclc_divergence_as_a_resistance_mechanism_axis():
    """The mechanistic model joins the NSCLC divergence view alongside the
    phenomenological Claret model — two resistance mechanisms, one context."""
    ds = onkos.load()
    cmp = onkos.compare(ds, purpose="tgi", context=NSCLC, drug_effect=1.0)
    ids = {tr.record_id for tr in cmp.included}
    assert RID in ids
    assert "resistance.claret_2009.tgi" in ids        # the phenomenological contrast
    assert len(cmp.included) >= 4                      # claret + wang + norton-simon + two-pop


def test_out_of_context_floors_to_D():
    ds = onkos.load()
    tr = onkos.simulate(ds, RID, context={"tumor_type": "breast", "line": "first"},
                        drug_effect=1.0)
    assert tr.tier == "D"
    assert any("outside validated" in w for w in tr.warnings)


def test_resistant_burden_is_practically_unidentifiable():
    """R0 (the interpretable resistant fraction) is still poorly identified from a
    realistic trial — mechanistic does not mean measured (composes with v0.22)."""
    ds = onkos.load()
    res = onkos.identifiability(ds, RID, context=NSCLC)
    r0 = next(p for p in res.params if p.symbol == "R0")
    assert not r0.identifiable
