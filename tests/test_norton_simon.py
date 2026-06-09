"""Norton-Simon kill model (drug_effect subsystem; kill proportional to growth)."""

import numpy as np
import onkos
from onkos.export.reference import KERNELS
from onkos.export.registry import get_kernel, kernel_values

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
RID = "drug_effect.norton_simon.nsclc"


def test_kernel_registered_and_subsystem_filled():
    ds = onkos.load()
    assert "norton_simon" in KERNELS
    r = ds[RID]
    assert r.subsystem == "drug_effect"  # the previously-empty declared subsystem
    assert get_kernel(r).name == "norton_simon"


def test_kill_is_proportional_to_growth_rate():
    """Norton-Simon: a smaller (faster-growing Gompertz) tumor has a higher
    fractional kill rate than a large one near carrying capacity."""
    ds = onkos.load()
    v = kernel_values(ds[RID])
    v["E"] = 1.0
    spec = get_kernel(ds[RID])

    def fractional_kill(size):
        v["V0"] = size
        rate = spec.rhs(0.0, [size], v)[0]
        return -rate / size  # positive = shrinking

    assert fractional_kill(40.0) > fractional_kill(80.0) > fractional_kill(160.0)


def test_appears_in_nsclc_divergence_view():
    ds = onkos.load()
    cmp = onkos.compare(ds, purpose="tgi", context=NSCLC, drug_effect=1.0)
    assert RID in {tr.record_id for tr in cmp.included}
    assert len(cmp.included) >= 3  # claret + wang + norton-simon


def test_distinct_mechanism_from_log_kill():
    """Without resistance, Norton-Simon with kill > growth eradicates the tumor
    (deep, durable response), unlike the Claret resistance-driven regrowth."""
    ds = onkos.load()
    t = np.linspace(0, 156, 313)
    ns = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.0, t=t)
    claret = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0, t=t)
    # Norton-Simon tumor keeps shrinking (no regrowth); Claret regrows after nadir
    assert ns.tumor_size[-1] < ns.metrics["nadir_tumor_size"] + 1e-6
    assert claret.tumor_size[-1] > claret.metrics["nadir_tumor_size"] * 2


def test_out_of_context_floors_to_D():
    ds = onkos.load()
    tr = onkos.simulate(ds, RID, context={"tumor_type": "breast", "line": "first"}, drug_effect=1.0)
    assert tr.tier == "D"
    assert any("outside validated" in w for w in tr.warnings)


def test_higher_effect_deeper_response():
    ds = onkos.load()
    low = onkos.simulate(ds, RID, context=NSCLC, drug_effect=0.5)
    high = onkos.simulate(ds, RID, context=NSCLC, drug_effect=1.5)
    assert high.metrics["depth_of_response"] >= low.metrics["depth_of_response"] - 1e-9
    assert high.metrics["week8_relative_change"] < low.metrics["week8_relative_change"]
