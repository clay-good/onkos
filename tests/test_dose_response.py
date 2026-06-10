"""Exposure-response model choice — the dose-extrapolation model-selection axis (v0.36.0).

Emax, sigmoid-Emax, and power dose-response shapes are re-anchored to agree at a
reference exposure (the studied dose), so they are indistinguishable there. The finding:
their predicted effect — and the resulting OS — diverges as you extrapolate to an
untested dose, sharpest on de-escalation. A dose-response model fit at one dose carries
an unquantified model-selection risk when used to pick another.

Closed-form landmarks pin the calibration (all shapes hit the anchor exactly); the
binding landmarks pin the extrapolation finding.
"""

import numpy as np
import onkos
from onkos.dose_response import (
    ER_SHAPE_RECORDS,
    ExtrapolationComparison,
    calibrated_er,
    compare_er_extrapolation,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CLARET = "resistance.claret_2009.tgi"
C_REF, E_REF = 150.0, 1.0


# ---- calibration landmarks -------------------------------------------------

def test_every_shape_hits_the_anchor_exactly():
    """Each re-anchored ER shape passes through (c_ref, e_ref) to machine precision — the
    construction that makes the curves indistinguishable at the studied dose."""
    ds = onkos.load()
    for er in ER_SHAPE_RECORDS:
        f = calibrated_er(ds, er, c_ref=C_REF, e_ref=E_REF)
        assert np.isclose(f(C_REF), E_REF, rtol=1e-12)


def test_calibration_preserves_shape_monotonicity():
    """A calibrated ER curve is increasing in dose (more drug -> more effect)."""
    ds = onkos.load()
    for er in ER_SHAPE_RECORDS:
        f = calibrated_er(ds, er, c_ref=C_REF, e_ref=E_REF)
        assert f(0.5 * C_REF) < f(C_REF) < f(2.0 * C_REF)


def test_shapes_diverge_off_the_anchor():
    """Anchored at one point, the saturating / unbounded / switch-like shapes give
    DIFFERENT effects at other doses — that is the whole risk."""
    ds = onkos.load()
    fs = [calibrated_er(ds, er, c_ref=C_REF, e_ref=E_REF) for er in ER_SHAPE_RECORDS]
    half = [f(0.5 * C_REF) for f in fs]
    assert max(half) - min(half) > 0.05  # a real spread below the anchor
    dbl = [f(2.0 * C_REF) for f in fs]
    assert max(dbl) - min(dbl) > 0.05  # and above it


def test_calibration_rejects_nonpositive_anchor():
    ds = onkos.load()
    for bad in [(0.0, 1.0), (150.0, 0.0), (-1.0, 1.0)]:
        try:
            calibrated_er(ds, ER_SHAPE_RECORDS[0], c_ref=bad[0], e_ref=bad[1])
            raise AssertionError("expected ValueError for non-positive anchor")
        except ValueError:
            pass


# ---- the extrapolation finding --------------------------------------------

def test_no_divergence_at_the_studied_dose():
    """The control: at the anchor dose the ER-model choice has ZERO OS consequence —
    so any divergence elsewhere is attributable to extrapolation, not to the models
    being different overall."""
    cmp = compare_er_extrapolation(onkos.load(), CLARET, context=NSCLC, c_ref=C_REF, e_ref=E_REF)
    assert isinstance(cmp, ExtrapolationComparison)
    assert np.isclose(cmp.reference_os_divergence, 0.0, atol=1e-9)
    # effects all equal at the anchor:
    anchor = min(cmp.rows, key=lambda r: abs(r["dose"] - C_REF))
    assert np.isclose(anchor["effect_divergence"], 0.0, atol=1e-9)


def test_divergence_appears_on_extrapolation():
    """Off the anchor the OS prediction depends on the ER model — a real, nonzero spread."""
    cmp = compare_er_extrapolation(onkos.load(), CLARET, context=NSCLC, c_ref=C_REF, e_ref=E_REF)
    assert cmp.max_os_divergence > 5.0  # weeks of OS riding on the ER-shape choice alone
    assert cmp.os_divergence_at(0.25 * C_REF) > cmp.reference_os_divergence


def test_downward_extrapolation_diverges_most_in_os():
    """De-escalation (lower dose) sits on the steep part of the effect->OS relationship,
    so the ER-model choice matters most there — the clinically pointed case."""
    cmp = compare_er_extrapolation(onkos.load(), CLARET, context=NSCLC, c_ref=C_REF, e_ref=E_REF)
    down = cmp.os_divergence_at(0.25 * C_REF)
    up = cmp.os_divergence_at(4.0 * C_REF)
    assert down > up


def test_effect_divergence_grows_below_the_anchor():
    """Effect spread is ~0 at the anchor and larger as you move away from it (down)."""
    cmp = compare_er_extrapolation(onkos.load(), CLARET, context=NSCLC, c_ref=C_REF, e_ref=E_REF)
    assert cmp.effect_divergence_at(0.25 * C_REF) > cmp.effect_divergence_at(C_REF)


def test_single_model_has_no_divergence():
    """With one ER shape there is nothing to disagree about — divergence is identically 0."""
    cmp = compare_er_extrapolation(
        onkos.load(), CLARET, context=NSCLC, er_ids=["exposure_response.emax_generic"],
        c_ref=C_REF, e_ref=E_REF,
    )
    assert cmp.max_os_divergence == 0.0


def test_inherits_tier_and_transport():
    """A dose-response analysis rides the chain's tier; out of context it floors to D."""
    ds = onkos.load()
    assert compare_er_extrapolation(ds, CLARET, context=NSCLC, c_ref=C_REF, e_ref=E_REF).tier == "C"
    out = compare_er_extrapolation(
        ds, CLARET, context={"tumor_type": "CRC", "line": "first"}, c_ref=C_REF, e_ref=E_REF
    )
    assert out.tier == "D"


def test_to_json_carries_clinical_use():
    cmp = compare_er_extrapolation(onkos.load(), CLARET, context=NSCLC, c_ref=C_REF, e_ref=E_REF)
    import json

    d = json.loads(cmp.to_json())
    assert d["onkos:clinicalUse"].startswith("PROHIBITED")
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert d["reference_os_divergence"] == 0.0
