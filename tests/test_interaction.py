"""Drug-combination interaction landmark suite — the combination layer *is* the
standard interaction nulls and a monotone interaction index, not an unconstrained
synergy knob.

Mirrors ``test_combine.py`` / ``test_identifiability.py``: closed-form properties of
the combination rules themselves (research spec §5), plus integration checks through
``compare_interactions`` / ``simulate_combination``.
"""

import numpy as np
import onkos
from onkos.interaction import (
    INTERACTION_MODELS,
    bliss_fraction,
    combine_effects,
    compare_interactions,
    simulate_combination,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


# --- pure interaction math -------------------------------------------------


def test_greco_zero_psi_is_the_additive_null():
    for a, b in [(0.6, 0.4), (0.2, 1.1), (0.0, 0.5)]:
        assert np.isclose(combine_effects(a, b, model="greco", psi=0.0),
                          combine_effects(a, b, model="additive"))
        assert np.isclose(combine_effects(a, b, model="additive"), a + b)


def test_bliss_independence_equals_additive_for_loglinear_kill():
    """Combining the two log-linear survival fractions (Bliss) gives effect a+b."""
    for a, b in [(0.6, 0.4), (1.2, 0.3), (0.05, 2.0)]:
        assert np.isclose(bliss_fraction(a, b), 1.0 - np.exp(-(a + b)))


def test_interaction_model_ordering():
    """hsa <= additive, and greco(-) <= additive <= greco(+) for non-negative effects."""
    for a, b in [(0.6, 0.4), (0.9, 0.9), (0.1, 2.0)]:
        hsa = combine_effects(a, b, model="hsa")
        add = combine_effects(a, b, model="additive")
        syn = combine_effects(a, b, model="greco", psi=0.5)
        ant = combine_effects(a, b, model="greco", psi=-0.5)
        assert hsa <= add + 1e-12
        assert ant <= add + 1e-12 <= syn + 1e-12


def test_combined_effect_is_monotone_in_psi():
    prev = None
    for psi in np.linspace(-1.0, 2.0, 13):
        e = combine_effects(0.6, 0.6, model="greco", psi=psi)
        if prev is not None:
            assert e >= prev - 1e-12
        prev = e


def test_single_agent_degeneracy_no_manufactured_interaction():
    """If one agent is inactive, every model reduces to monotherapy."""
    for model in INTERACTION_MODELS:
        assert np.isclose(combine_effects(0.7, 0.0, model=model, psi=0.9), 0.7)
        assert np.isclose(combine_effects(0.0, 0.7, model=model, psi=0.9), 0.7)


def test_symmetry():
    for model in INTERACTION_MODELS:
        assert np.isclose(combine_effects(0.6, 0.3, model=model, psi=0.4),
                          combine_effects(0.3, 0.6, model=model, psi=0.4))


def test_monotone_in_each_effect():
    base = combine_effects(0.5, 0.5, model="additive")
    assert combine_effects(0.8, 0.5, model="additive") >= base
    assert combine_effects(0.5, 0.8, model="greco", psi=0.5) >= \
        combine_effects(0.5, 0.5, model="greco", psi=0.5)


def test_antagonism_is_floored_at_zero():
    """A strongly antagonistic interaction never yields a negative combined effect."""
    assert combine_effects(0.6, 0.6, model="greco", psi=-100.0) == 0.0


def test_unknown_model_rejected():
    try:
        combine_effects(0.5, 0.5, model="loewe_dose")
        raise AssertionError("expected ValueError for an unknown interaction model")
    except ValueError as e:
        assert "unknown interaction model" in str(e)


# --- integration through the simulation bridge -----------------------------


def test_interaction_divergence_positive_with_synergy_zero_without_partner():
    """ψ≠0 with both agents active ⇒ the interaction assumption moves the OS curve;
    an inactive partner ⇒ no divergence (no manufactured interaction)."""
    ds = onkos.load()
    cmp = compare_interactions(ds, "resistance.claret_2009.tgi", context=NSCLC,
                               effect_a=0.6, effect_b=0.6, psi=0.5)
    assert cmp.os_divergence > 0.0
    lo, hi = cmp.median_os_range
    assert hi > lo  # the interaction assumption alone spreads median OS

    mono = compare_interactions(ds, "resistance.claret_2009.tgi", context=NSCLC,
                                effect_a=0.6, effect_b=0.0, psi=0.5)
    assert np.isclose(mono.os_divergence, 0.0)
    assert len(set(round(v, 9) for v in mono.combined_effects.values())) == 1


def test_combination_tier_passes_through_underlying_chain():
    """The interaction assumption cannot raise the tier; it equals the tier of a plain
    simulation at the combined effect."""
    ds = onkos.load()
    tr = simulate_combination(ds, "resistance.claret_2009.tgi", context=NSCLC,
                              effect_a=0.6, effect_b=0.6, interaction="additive")
    plain = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.2)
    assert tr.tier == plain.tier
    # An out-of-context transport still floors to D through the combination bridge.
    out = simulate_combination(ds, "resistance.claret_2009.tgi",
                               context={"tumor_type": "melanoma", "line": "first"},
                               effect_a=0.6, effect_b=0.6, interaction="additive")
    assert out.tier == "D"


def test_synergy_warning_is_attached_only_when_assumed():
    ds = onkos.load()
    syn = simulate_combination(ds, "resistance.claret_2009.tgi", context=NSCLC,
                               effect_a=0.6, effect_b=0.6, interaction="greco", psi=0.5)
    assert any("synergy_is_an_assumption" in w for w in syn.warnings)
    add = simulate_combination(ds, "resistance.claret_2009.tgi", context=NSCLC,
                               effect_a=0.6, effect_b=0.6, interaction="additive")
    assert not any("synergy_is_an_assumption" in w for w in add.warnings)


def test_result_carries_clinical_use_and_synergy_note():
    ds = onkos.load()
    d = compare_interactions(ds, "resistance.claret_2009.tgi", context=NSCLC,
                             effect_a=0.6, effect_b=0.6, psi=0.5).to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "PROHIBITED" in d["onkos:clinicalUse"]
    assert "combination trial" in d["synergy_note"]
    assert set(d["combined_effects"]) == set(d["median_os_weeks"])
