"""Simulation-based logrank power — the check on the closed-form Schoenfeld core.

Follow-on to ``docs/specs/research/model-discriminability.md`` §7. Landmarks: hand-computed
logrank statistics, type-I error control, calibration against Schoenfeld under proportional
hazards, and the crossing-curve case where the closed form misjudges a pair.
"""

import numpy as np
import onkos
import pytest
from onkos.discriminability import (
    feasibility_class,
    horizon_hazard_ratio,
    logrank_statistic,
    required_events,
    sample_event_times,
    simulated_discriminability,
    simulated_power,
    simulated_required_events,
)

T = np.linspace(0.0, 312.0, 625)


# --- the logrank statistic, by hand ------------------------------------------


def test_logrank_matches_a_hand_computation():
    # group A dies at 1, 2; group B at 3, 4. At t=1: n=4, nA=2 -> E=1/2, V=1/4.
    # At t=2: n=3, nA=1 -> E=1/3, V=2/9. Later times carry no group-A risk.
    z = logrank_statistic([1, 2, 3, 4], [1, 1, 1, 1], [1, 1, 0, 0])
    expected = (2 - (1 / 2 + 1 / 3)) / np.sqrt(1 / 4 + 2 / 9)
    assert np.isclose(z, expected)


def test_logrank_handles_tied_event_times():
    # t=1: n=4, nA=2, d=3 (2 in A) -> E=3/2, V=3*(1/2)(1/2)(1)/3 = 1/4.
    z = logrank_statistic([1, 1, 1, 2], [1, 1, 1, 1], [1, 1, 0, 0])
    assert np.isclose(z, (2 - 1.5) / 0.5)


def test_logrank_is_antisymmetric_and_censoring_only_shapes_risk_sets():
    time = [1, 2, 3, 4, 5]
    event = [1, 0, 1, 1, 0]
    group = np.array([1, 1, 0, 0, 1], dtype=bool)
    z = logrank_statistic(time, event, group)
    assert np.isclose(logrank_statistic(time, event, ~group), -z)
    assert logrank_statistic([1, 2], [0, 0], [1, 0]) == 0.0  # no events, no information


# --- sampling and power -------------------------------------------------------


def test_event_times_follow_the_curve():
    rng = np.random.default_rng(1)
    s = np.exp(-T / 100.0)
    draws = sample_event_times(T, s, 20000, rng)
    assert np.isclose(np.mean(np.isinf(draws)), s[-1], atol=0.01)  # no event by horizon
    finite = draws[np.isfinite(draws)]
    assert np.isclose(np.median(finite), -100.0 * np.log((1 + s[-1]) / 2), rtol=0.03)


def test_identical_curves_hold_the_type_one_error():
    s = np.exp(-T / 80.0)
    r = simulated_power(s, s, T, n_per_arm=150, reps=600, seed=3)
    assert abs(r["power"] - 0.05) < 3 * np.sqrt(0.05 * 0.95 / 600)


def test_simulation_reproduces_schoenfeld_under_proportional_hazards():
    sa, sb = np.exp(-T / 60.0), np.exp(-T / 80.0)
    d = required_events(horizon_hazard_ratio(sa, sb))
    n = int(np.ceil(d / (2.0 - sa[-1] - sb[-1])))
    r = simulated_power(sa, sb, T, n_per_arm=n, reps=600, seed=0)
    assert abs(r["power"] - 0.8) < 0.06
    sim = simulated_required_events(sa, sb, T, reps=300, seed=0)
    assert sim["feasible"] and 0.8 <= d / sim["events"] <= 1.25


def test_too_few_events_cannot_reach_power():
    s = np.exp(-T / 80.0)
    r = simulated_required_events(s, s, T, reps=50, max_events=200)
    assert not r["feasible"] and r["events"] == float("inf")


def test_horizon_hr_ignores_curves_that_have_both_vanished():
    """Two PH Weibull curves that both fall below the floor before the horizon used to read
    HR = 1 (both clamped). The ratio is now read where both are still observable."""
    h0 = (T / 20.0) ** 1.3
    sa, sb = np.exp(-h0), np.exp(-2.0 * h0)
    assert sa[-1] < 1e-6 and sb[-1] < 1e-6
    assert np.isclose(horizon_hazard_ratio(sb, sa), 2.0, rtol=1e-6)


def test_feasibility_classes():
    assert feasibility_class(100) == "feasible"
    assert feasibility_class(1000) == "borderline"
    assert feasibility_class(float("inf")) == "infeasible"


# --- integration ---------------------------------------------------------------


def test_proportional_hazard_links_agree_with_schoenfeld():
    ds = onkos.load()
    sd = simulated_discriminability(ds, context={"tumor_type": "NSCLC", "line": "first"},
                                    reps=150)
    assert sd.pairs and not sd.misclassified_pairs
    ratios = [p["events_ratio"] for p in sd.pairs if p["events_ratio"] is not None]
    assert ratios and all(0.7 <= r <= 1.3 for r in ratios)


def test_crossing_joint_curves_break_the_closed_form():
    """Breast, joint link: Claret and two-population cross near week 266. The horizon HR
    (~0.63) promises a ~150-event trial; the logrank test has almost no power."""
    ds = onkos.load()
    sd = simulated_discriminability(ds, context={"tumor_type": "breast", "line": "first"},
                                    joint_alpha=1.0, reps=150)
    bad = {(p["record_a"], p["record_b"]) for p in sd.misclassified_pairs}
    assert ("resistance.breast_first_line.claret",
            "resistance.breast_first_line.two_population") in bad
    p = next(p for p in sd.pairs if p["record_b"] == "resistance.breast_first_line.two_population"
             and p["record_a"] == "resistance.breast_first_line.claret")
    assert p["schoenfeld_class"] == "feasible" and p["simulated_class"] == "infeasible"
    assert sd.to_dict()["NOT_FOR_CLINICAL_USE"] is True


def test_simulated_discriminability_rejects_bad_inputs():
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    with pytest.raises(ValueError):
        simulated_discriminability(ds, context=ctx, reps=5)
    with pytest.raises(ValueError):
        simulated_discriminability(ds, context=ctx, accrual=400.0)
    with pytest.raises(ValueError):
        simulated_discriminability(ds, context=ctx, joint_alpha=1.0,
                                   survival_link="survival_link.nsclc_os_cox")
    with pytest.raises(ValueError):
        logrank_statistic([1.0, np.nan], [1, 1], [1, 0])
