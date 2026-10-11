"""Forward simulation and the virtual-trial divergence view."""

import numpy as np
import onkos
import pytest
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


def _ode_record_ids():
    from onkos.export.registry import get_kernel

    return [r.id for r in onkos.load() if r.kernel and get_kernel(r).kind == "ode"]


@pytest.mark.parametrize("rid", _ode_record_ids())
def test_time_grid_is_absolute_time_since_baseline(rid):
    """``t`` means weeks since baseline for every kernel — closed-form or numerically
    integrated. A grid that starts late, is unsorted, or repeats a time must return the
    same values as the reference grid (regression: integrated multi-state models used to
    treat the first grid point as baseline, silently shifting the trajectory)."""
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    ref = onkos.simulate(ds, rid, context=ctx, t=np.array([0.0, 4.0, 8.0, 24.0])).tumor_size
    late = onkos.simulate(ds, rid, context=ctx, t=np.array([4.0, 8.0, 24.0])).tumor_size
    shuffled = onkos.simulate(ds, rid, context=ctx, t=np.array([24.0, 0.0, 8.0, 8.0])).tumor_size
    baseline_only = onkos.simulate(ds, rid, context=ctx, t=np.array([0.0])).tumor_size
    assert np.allclose(late, ref[1:], rtol=1e-5)
    assert np.allclose(shuffled, ref[[3, 0, 2, 2]], rtol=1e-5)
    assert np.allclose(baseline_only, ref[:1])


def test_integrated_kernels_reject_negative_time():
    ds = onkos.load()
    with pytest.raises(ValueError, match="times must be >= 0"):
        onkos.simulate(ds, "growth_laws.simeoni_exp_linear", t=np.array([-4.0, 0.0, 8.0]))
