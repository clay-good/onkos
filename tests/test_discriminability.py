"""Model discriminability — can a trial even tell the competing models apart? (v0.38.0)

Given two models' population OS curves, the required number of events to distinguish them
by a logrank test is the Schoenfeld formula d = 4(z_{1-a/2}+z_{1-b})^2 / (ln HR)^2. A
divergence that needs tens of thousands of events is practically unidentifiable from the
trial — the model choice can only be assumed, not resolved. That reframes the whole
model-selection arc: the silent risks (resistance mechanism v0.24, origin v0.32) are silent
because the week-8 OS surrogate cannot distinguish them.

Closed-form landmarks pin the power core; the binding landmarks pin the (un)detectability
finding.
"""

import numpy as np
import onkos
from onkos.discriminability import (
    ModelDiscriminability,
    horizon_hazard_ratio,
    model_discriminability,
    required_events,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CLARET = "resistance.claret_2009.tgi"
TWO_POP = "resistance.nsclc_first_line.two_population"
ACQUIRED = "resistance.nsclc_first_line.acquired"
NORTON = "drug_effect.norton_simon.nsclc"


# ---- pure-core closed-form landmarks --------------------------------------

def test_schoenfeld_matches_the_textbook_value():
    """HR=0.5 at 80% power, two-sided alpha 0.05 needs ~65 events (the standard logrank
    sample-size benchmark)."""
    d = required_events(0.5, power=0.8, alpha=0.05)
    assert 63.0 < d < 68.0


def test_required_events_is_symmetric_in_hr_inverse():
    """Distinguishing HR and 1/HR needs the same events (depends only on |ln HR|)."""
    assert np.isclose(required_events(0.6), required_events(1.0 / 0.6))


def test_identical_curves_need_infinite_events():
    assert required_events(1.0) == float("inf")
    assert required_events(0.0) == float("inf")
    assert required_events(float("nan")) == float("inf")


def test_smaller_divergence_needs_more_events():
    """A hazard ratio closer to 1 (smaller divergence) needs more events — monotone."""
    assert required_events(0.5) < required_events(0.8) < required_events(0.95)


def test_more_power_or_smaller_alpha_needs_more_events():
    assert required_events(0.6, power=0.9) > required_events(0.6, power=0.8)
    assert required_events(0.6, alpha=0.01) > required_events(0.6, alpha=0.05)


def test_horizon_hazard_ratio_recovers_proportional_hazards():
    """For two genuinely proportional-hazards curves S_b = S_a**HR, the horizon ratio of
    cumulative hazards recovers HR exactly."""
    t = np.linspace(0.0, 200.0, 401)
    sa = np.exp(-((t / 80.0) ** 1.2))
    hr_true = 1.7
    sb = sa ** hr_true
    assert np.isclose(horizon_hazard_ratio(sb, sa), hr_true, rtol=1e-6)


def test_identical_curves_have_unit_hazard_ratio():
    t = np.linspace(0.0, 200.0, 401)
    sa = np.exp(-((t / 80.0) ** 1.2))
    assert np.isclose(horizon_hazard_ratio(sa, sa), 1.0)


# ---- the (un)detectability finding ----------------------------------------

def test_resistance_models_are_practically_indistinguishable_under_week8():
    """The v0.24/v0.32 finding, quantified: under the week-8 OS surrogate the resistance
    models (which diverge only in the regrowth tail) need an INFEASIBLE trial to tell apart,
    while a model that differs in early shrinkage is easily distinguished."""
    md = model_discriminability(onkos.load(), context=NSCLC)
    assert isinstance(md, ModelDiscriminability)

    def events(a, b):
        for p in md.pairs:
            if {p["record_a"], p["record_b"]} == {a, b}:
                return p["required_events"]
        raise KeyError

    # resistance mechanism / origin: undetectable under week-8
    assert events(CLARET, TWO_POP) > 3000
    assert events(CLARET, ACQUIRED) > 3000
    # but the complete responder (deep, eradicating) is easily distinguished from a resistance model
    assert events(CLARET, NORTON) < 500
    # and the within-resistance pairs need vastly more than the early-shrinkage-distinct pairs
    assert events(CLARET, TWO_POP) > 10 * events(CLARET, NORTON)


def test_several_pairs_are_flagged_indistinguishable():
    md = model_discriminability(onkos.load(), context=NSCLC)
    assert md.n_indistinguishable >= 3
    assert len(md.feasible_pairs) >= 1  # not everything is indistinguishable


def test_metric_choice_consequence_is_easily_detectable():
    """Contrast: the survival-METRIC choice (week-8 vs k_g for the SAME model) produces a
    large OS swing that a small trial detects — so the metric's consequence is identifiable
    even though the model choice it is blind to is not."""
    from onkos.discriminability import discriminating_events

    ds = onkos.load()
    t = np.linspace(0.0, 312.0, 625)
    w = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=t).os_curve
    g = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=t,
                       survival_link="survival_link.nsclc_os_growth_rate").os_curve
    assert required_events(horizon_hazard_ratio(w, g)) < 500
    # discriminating_events convenience wrapper agrees it is feasible between two models too
    de = discriminating_events(ds, CLARET, NORTON, context=NSCLC)
    assert de["required_events"] < 500


def test_inherits_tier_and_carries_clinical_use():
    md = model_discriminability(onkos.load(), context=NSCLC)
    assert md.tier in ("A", "B", "C", "D")
    import json

    d = json.loads(md.to_json())
    assert d["onkos:clinicalUse"].startswith("PROHIBITED")
    assert d["NOT_FOR_CLINICAL_USE"] is True
    # infinite required-events serialize as null, not a NaN/inf literal
    assert all(p["required_events"] is None or p["required_events"] > 0 for p in d["pairs"])
