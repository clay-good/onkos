"""Forward simulation and the virtual-trial divergence view."""

import numpy as np
import onkos
from onkos.tiers import propagate, worst_tier


def test_worst_tier_wins():
    assert worst_tier(["A", "B", "C"]) == "C"
    assert worst_tier(["A", "A"]) == "A"
    assert worst_tier([]) == "D"


def test_simulate_claret_in_context():
    ds = onkos.load()
    tr = onkos.simulate(
        ds, "resistance.claret_2009.tgi", context={"tumor_type": "NSCLC", "line": "first"}
    )
    assert tr.tumor_size.shape == tr.t.shape
    assert tr.os_curve is not None
    assert 0.0 < tr.os_curve[0] <= 1.0
    assert np.all(np.diff(tr.os_curve) <= 1e-9)  # survival is monotone non-increasing
    assert tr.median_os is not None
    assert tr.tier in ("B", "C")  # in-context: not floored to D


def test_out_of_context_floors_to_D():
    ds = onkos.load()
    # NSCLC-validated model applied to a colorectal context -> tier floor D + warning
    res = propagate(
        [ds["resistance.claret_2009.tgi"]], tumor_type="colorectal", line="first"
    )
    assert res.tier == "D"
    assert any("outside validated" in w for w in res.warnings)


def test_compare_excludes_out_of_context_models():
    ds = onkos.load()
    cmp = onkos.compare(
        ds, purpose="tgi", context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=1.0
    )
    included = {tr.record_id for tr in cmp.included}
    excluded = {rid for rid, _ in cmp.excluded}
    assert "resistance.claret_2009.tgi" in included
    assert "tgi_metrics.wang_2009.biexponential" in included
    assert "tgi_metrics.bruno_2020.breast_biexponential" in excluded


def test_divergence_is_measurable():
    ds = onkos.load()
    cmp = onkos.compare(
        ds, purpose="tgi", context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=1.0
    )
    assert cmp.os_divergence > 0.0
    lo, hi = cmp.median_os_range
    assert hi > lo  # model choice changes the survival prediction


def test_deeper_response_improves_survival():
    ds = onkos.load()
    weak = onkos.simulate(
        ds, "tgi_metrics.wang_2009.biexponential",
        context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=0.5,
    )
    strong = onkos.simulate(
        ds, "tgi_metrics.wang_2009.biexponential",
        context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=1.5,
    )
    assert strong.metrics["depth_of_response"] > weak.metrics["depth_of_response"]
    assert strong.median_os > weak.median_os
