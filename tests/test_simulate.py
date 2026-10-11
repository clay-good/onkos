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
    integrated. A grid that starts late or repeats a time must return the
    same values as the reference grid (regression: integrated multi-state models used to
    treat the first grid point as baseline, silently shifting the trajectory)."""
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    ref = onkos.simulate(ds, rid, context=ctx, t=np.array([0.0, 4.0, 8.0, 24.0])).tumor_size
    late = onkos.simulate(ds, rid, context=ctx, t=np.array([4.0, 8.0, 24.0])).tumor_size
    repeated = onkos.simulate(ds, rid, context=ctx, t=np.array([0.0, 8.0, 8.0, 24.0])).tumor_size
    baseline_only = onkos.simulate(ds, rid, context=ctx, t=np.array([0.0])).tumor_size
    assert np.allclose(late, ref[1:], rtol=1e-5)
    assert np.allclose(repeated, ref[[0, 2, 2, 3]], rtol=1e-5)
    assert np.allclose(baseline_only, ref[:1])


def test_integrated_kernels_reject_negative_time():
    ds = onkos.load()
    with pytest.raises(ValueError, match="times must be >= 0"):
        onkos.simulate(ds, "growth_laws.simeoni_exp_linear", t=np.array([-4.0, 0.0, 8.0]))


CLARET = "resistance.claret_2009.tgi"
NSCLC = {"tumor_type": "NSCLC", "line": "first"}


@pytest.mark.parametrize(
    "t", [np.array([8.0, 0.0, 24.0]), np.array([]), np.array([0.0, np.nan]), np.zeros((2, 2))],
    ids=["unsorted", "empty", "nan", "2-D"],
)
def test_simulate_rejects_malformed_time_grids(t):
    """Metrics and survival assume an ordered grid; a shuffled one used to corrupt them
    silently (median OS off by weeks, burden AUC sign-flipped)."""
    with pytest.raises(ValueError, match="non-decreasing"):
        onkos.simulate(onkos.load(), CLARET, context=NSCLC, t=t)


@pytest.mark.parametrize(
    "t", [np.linspace(0.0, 4.0, 9), np.linspace(12.0, 104.0, 185)], ids=["ends-early", "starts-late"]
)
def test_week8_metrics_require_week8_on_the_grid(t):
    """A grid that does not reach week 8 must not pass off its nearest end as the week-8
    value, and a link driven by that covariate yields no curve (with a warning) instead of
    a curve drawn from the wrong number."""
    tr = onkos.simulate(onkos.load(), CLARET, context=NSCLC, t=t)
    assert np.isnan(tr.metrics["week8_relative_change"])
    assert np.isnan(tr.metrics["week8_tumor_size"])
    assert tr.os_curve is None and tr.pfs_curve is None
    assert sum("does not cover week 8" in w for w in tr.warnings) == 2  # the OS and PFS links


def test_tail_metric_link_survives_a_late_starting_grid():
    """Only week-8-driven links need week 8: the k_g link still yields its OS curve."""
    tr = onkos.simulate(
        onkos.load(), CLARET, context=NSCLC, t=np.linspace(12.0, 156.0, 289),
        survival_link="survival_link.nsclc_os_growth_rate",
    )
    assert tr.os_curve is not None and not any("week 8" in w for w in tr.warnings)
