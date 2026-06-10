"""Joint longitudinal–survival modeling — the current-value link (v0.34.0).

Every two-stage link (week-8, k_g, burden) collapses the trajectory to a static
covariate, so the hazard ratio is constant in time — a proportional hazard. The
joint current-value link makes the instantaneous hazard track the current tumor
size, producing a TIME-VARYING hazard ratio (a non-proportional hazard) the
two-stage links cannot represent.

Closed-form landmarks pin the pure core (the two limits where it must reduce to a
known model exactly); the binding landmarks pin the finding.
"""

import numpy as np
import onkos
from onkos.joint import (
    JointComparison,
    compare_joint_vs_two_stage,
    current_value_survival,
    joint_survival,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CLARET = "resistance.claret_2009.tgi"
TWO_POP = "resistance.nsclc_first_line.two_population"
NORTON = "drug_effect.norton_simon.nsclc"
WANG = "tgi_metrics.wang_2009.biexponential"
ACQUIRED = "resistance.nsclc_first_line.acquired"
T = np.linspace(0.0, 260.0, 521)
SHAPE, SCALE = 1.3, 60.0  # the NSCLC week-8 Weibull baseline


# ---- pure-core closed-form landmarks --------------------------------------

def test_constant_hazard_ratio_recovers_two_stage_exactly():
    """A constant HR telescopes to the two-stage Weibull-PH curve EXACTLY: the joint
    model is a strict generalization of proportional hazards, not an approximation."""
    for c in (0.3, 1.0, 2.5):
        hr = np.full_like(T, c)
        joint = current_value_survival(T, hr, shape=SHAPE, scale=SCALE)
        two_stage = np.exp(-((T / SCALE) ** SHAPE) * c)  # HR=c constant
        assert np.allclose(joint, two_stage, atol=1e-12)


def test_unit_hazard_ratio_recovers_weibull_baseline():
    """HR ≡ 1 (tumor held at baseline size) recovers the Weibull baseline survival."""
    joint = current_value_survival(T, np.ones_like(T), shape=SHAPE, scale=SCALE)
    baseline = np.exp(-((T / SCALE) ** SHAPE))
    assert np.allclose(joint, baseline, atol=1e-12)


def test_survival_is_monotone_decreasing_and_starts_at_one():
    hr = 1.0 + 0.5 * np.sin(T / 20.0) ** 2  # any positive HR
    s = current_value_survival(T, hr, shape=SHAPE, scale=SCALE)
    assert np.isclose(s[0], 1.0)
    assert np.all(np.diff(s) <= 1e-12)


def test_higher_hazard_ratio_lowers_survival_pointwise():
    s_lo = current_value_survival(T, np.full_like(T, 0.5), shape=SHAPE, scale=SCALE)
    s_hi = current_value_survival(T, np.full_like(T, 2.0), shape=SHAPE, scale=SCALE)
    assert np.all(s_hi <= s_lo + 1e-12)


# ---- binding & the finding -------------------------------------------------

def test_alpha_zero_removes_the_association():
    """alpha = 0 ⇒ HR ≡ 1 regardless of the tumor ⇒ the Weibull baseline survival."""
    j = joint_survival(onkos.load(), CLARET, context=NSCLC, alpha=0.0, t=T)
    assert np.allclose(j.hazard_ratio, 1.0)
    assert np.allclose(j.os_curve, np.exp(-((T / SCALE) ** SHAPE)), atol=1e-9)


def test_constant_size_trajectory_matches_a_two_stage_burden_link():
    """For a (hypothetical) constant-size trajectory the joint median equals the two-stage
    Weibull-PH median with x = log(size ratio), beta = alpha — the analytic bridge to v0.33.
    Verified directly on the pure core to avoid needing a flat-trajectory record."""
    c, alpha = 0.5, 1.0
    hr = np.full_like(T, c**alpha)
    joint = current_value_survival(T, hr, shape=SHAPE, scale=SCALE)
    two_stage = np.exp(-((T / SCALE) ** SHAPE) * np.exp(alpha * np.log(c)))
    assert np.allclose(joint, two_stage, atol=1e-12)


def test_regrowing_tumor_has_a_non_proportional_hazard():
    """The joint HR rises as a resistant clone regrows — ph_violation ≫ 1 — whereas a
    two-stage link has a constant HR (ph_violation == 1 by construction)."""
    j = joint_survival(onkos.load(), TWO_POP, context=NSCLC, alpha=1.0, t=T)
    assert j.hazard_ratio_at(8.0) < 1.0  # deep early response: hazard suppressed
    assert j.hazard_ratio[-1] > 5.0  # heavy regrowth tail: hazard elevated
    assert j.ph_violation > 5.0


def test_eradicating_responder_has_a_vanishing_tail_hazard():
    """A complete responder's current-value hazard collapses (HR → small) in the tail, so
    its joint OS exceeds every regrowing model's — the dynamic reward for durable response."""
    j = joint_survival(onkos.load(), NORTON, context=NSCLC, alpha=1.0, t=T)
    assert j.hazard_ratio[-1] < j.hazard_ratio_at(8.0)  # hazard keeps falling, not rising
    assert j.ph_violation < 1.0


def test_joint_link_re_ranks_relative_to_two_stage():
    """The current-value link re-orders the models versus the week-8 two-stage surrogate
    (it weights the regrowth tail), so two-stage-vs-joint is a real model-selection axis."""
    cmp = compare_joint_vs_two_stage(onkos.load(), context=NSCLC, alpha=1.0, t=T)
    assert isinstance(cmp, JointComparison)
    assert cmp.rank_discordant_pairs >= 1
    # the resistance models drive the largest PH violations:
    assert cmp.max_ph_violation > 10.0


def test_joint_demotes_the_heavier_regrowth_tail():
    """Matched at week 8 (both deep responders), the heavier-tail model is RANKED LOWER by
    the joint link than the phenomenological one — the opposite of the week-8 surrogate,
    which ranks the deep-early-shrinker on top."""
    ds = onkos.load()
    claret = joint_survival(ds, CLARET, context=NSCLC, alpha=1.0, t=T)
    two_pop = joint_survival(ds, TWO_POP, context=NSCLC, alpha=1.0, t=T)
    # week-8 two-stage ranks two-population above Claret ...
    assert two_pop.two_stage_median_os > claret.two_stage_median_os
    # ... the joint link inverts it (Claret's lighter tail wins):
    assert claret.median_os > two_pop.median_os


def test_joint_inherits_tier_and_transport_warnings():
    """A joint analysis never moves a tier; used out of context it floors to D + warns,
    exactly like the trajectory it summarizes."""
    ds = onkos.load()
    j = joint_survival(ds, CLARET, context={"tumor_type": "CRC", "line": "first"}, alpha=1.0, t=T)
    assert j.tier == "D"
    assert any("outside validated" in w for w in j.warnings)
    assert j.clinical_use.startswith("PROHIBITED")
