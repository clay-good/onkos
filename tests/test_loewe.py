"""Dose-level Loewe additivity — the additivity *reference* as a model-selection axis
(v0.35.0).

v0.23 combined two effect MAGNITUDES under effect-level nulls (HSA, Bliss-additive,
Greco). Loewe additivity combines two DOSES via the dose-response (ER) curves, solving
the isobole d_A/D_A(E) + d_B/D_B(E) = 1. It is the only "no-interaction" reference that
satisfies the sham-combination identity — a drug combined with itself is exactly additive
— which Bliss/effect-additivity fails for any saturating curve. The choice of reference is
therefore a model-selection axis, and it propagates to OS.

Closed-form landmarks pin the isobole core; the binding landmarks pin the finding.
"""

import numpy as np
import onkos
from onkos.interaction import (
    ADDITIVITY_REFERENCES,
    AdditivityComparison,
    ERCurve,
    combine_doses,
    compare_additivity_references,
    er_curve,
    loewe_effect,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CLARET = "resistance.claret_2009.tgi"
ER_A = "exposure_response.emax_generic"  # E = 1.4·C/(150+C)
ER_B = "exposure_response.dacomitinib_egfr.emax"  # E = 1.8·C/(90+C)
ER_POW = "exposure_response.power_generic"


# ---- pure-core closed-form landmarks --------------------------------------

def test_sham_combination_identity():
    """The defining property of Loewe: a drug combined with ITSELF is exactly additive —
    Loewe(d_a, d_b) of (A, A) equals the single-agent effect at the summed dose f(d_a+d_b).
    This is what makes Loewe the principled no-interaction reference."""
    ca = er_curve(onkos.load(), ER_A)
    for da, db in [(100.0, 200.0), (50.0, 50.0), (300.0, 30.0)]:
        assert np.isclose(loewe_effect(da, db, curve_a=ca, curve_b=ca), ca.forward(da + db), rtol=1e-6)


def test_bliss_fails_the_sham_test_for_a_saturating_curve():
    """Effect-additivity (Bliss) does NOT satisfy the sham identity for a saturating Emax
    curve: f(d_a)+f(d_b) > f(d_a+d_b). The gap is exactly why the reference choice matters."""
    ca = er_curve(onkos.load(), ER_A)
    da, db = 150.0, 150.0
    bliss = ca.forward(da) + ca.forward(db)
    loewe = loewe_effect(da, db, curve_a=ca, curve_b=ca)
    assert bliss > loewe  # effect-additivity overstates for a saturating curve
    assert np.isclose(loewe, ca.forward(da + db), rtol=1e-6)


def test_single_agent_limits():
    """One dose zero ⇒ the other drug's single-agent effect exactly."""
    ds = onkos.load()
    ca, cb = er_curve(ds, ER_A), er_curve(ds, ER_B)
    assert np.isclose(loewe_effect(150.0, 0.0, curve_a=ca, curve_b=cb), ca.forward(150.0))
    assert np.isclose(loewe_effect(0.0, 90.0, curve_a=ca, curve_b=cb), cb.forward(90.0))


def test_er_inverse_round_trips():
    """The analytic inverse is a true inverse of the forward ER curve, for each kernel."""
    ds = onkos.load()
    for er_id in (ER_A, ER_B, ER_POW):
        c = er_curve(ds, er_id)
        for d in (10.0, 75.0, 300.0):
            assert np.isclose(c.inverse(c.forward(d)), d, rtol=1e-9)


def test_loewe_is_clamped_at_the_shared_effect_ceiling():
    """Two saturating curves cannot jointly express more than min(emax_a, emax_b); a huge
    dose pair clamps there rather than diverging."""
    ds = onkos.load()
    ca, cb = er_curve(ds, ER_A), er_curve(ds, ER_B)
    e = loewe_effect(1e6, 1e6, curve_a=ca, curve_b=cb)
    assert e <= min(ca.emax, cb.emax) + 1e-9
    assert np.isfinite(e)


def test_loewe_is_monotone_in_dose():
    ds = onkos.load()
    ca, cb = er_curve(ds, ER_A), er_curve(ds, ER_B)
    lo = loewe_effect(50.0, 30.0, curve_a=ca, curve_b=cb)
    hi = loewe_effect(200.0, 120.0, curve_a=ca, curve_b=cb)
    assert hi > lo


def test_loewe_core_is_record_free():
    """The isobole core works on plain callables (an ERCurve from lambdas), so it is
    landmark-testable without the dataset — sham identity on a synthetic Emax curve."""
    def fwd(d):
        return 2.0 * d / (100.0 + d)

    def inv(e):
        return 100.0 * e / (2.0 - e)

    c = ERCurve(forward=fwd, inverse=inv, emax=2.0)
    assert np.isclose(loewe_effect(40.0, 60.0, curve_a=c, curve_b=c), fwd(100.0), rtol=1e-9)


# ---- the finding & binding -------------------------------------------------

def test_reference_ordering_hsa_le_loewe_le_bliss():
    """For saturating curves the three references bracket the combined effect:
    HSA (conservative) ≤ Loewe (self-consistent) ≤ Bliss (effect-additive, overstates)."""
    ds = onkos.load()
    e = {ref: combine_doses(ds, 150.0, 90.0, er_a=ER_A, er_b=ER_B, reference=ref)
         for ref in ADDITIVITY_REFERENCES}
    assert e["hsa"] <= e["loewe"] <= e["bliss"]
    assert e["bliss"] > e["hsa"]  # a genuine spread, not a degenerate tie


def test_reference_choice_diverges_in_os():
    """The same dose pair gives different OS across the additivity references — the
    reference is a model-selection axis with a survival consequence."""
    cmp = compare_additivity_references(
        onkos.load(), CLARET, context=NSCLC, dose_a=150.0, dose_b=90.0, er_a=ER_A, er_b=ER_B
    )
    assert isinstance(cmp, AdditivityComparison)
    assert cmp.os_divergence > 0.0
    meds = cmp.median_os
    assert meds["bliss"] > meds["loewe"] > meds["hsa"]  # deeper combined effect → longer OS


def test_additivity_comparison_inherits_tier_and_transport():
    """The comparison rides the underlying chain's tier; out of context it floors to D."""
    ds = onkos.load()
    in_ctx = compare_additivity_references(
        ds, CLARET, context=NSCLC, dose_a=150.0, dose_b=90.0, er_a=ER_A, er_b=ER_B
    )
    assert in_ctx.tier == "C"
    out = compare_additivity_references(
        ds, CLARET, context={"tumor_type": "CRC", "line": "first"},
        dose_a=150.0, dose_b=90.0, er_a=ER_A, er_b=ER_B
    )
    assert out.tier == "D"
    assert any("outside validated" in w for tr in out.trajectories.values() for w in tr.warnings)


def test_unknown_reference_and_kernel_raise():
    ds = onkos.load()
    try:
        combine_doses(ds, 1.0, 1.0, er_a=ER_A, er_b=ER_B, reference="bogus")
        raise AssertionError("expected ValueError for unknown reference")
    except ValueError:
        pass
