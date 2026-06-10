"""Survival-metric choice — which on-treatment TGI metric drives the OS link is
itself a model-selection axis.

The default link reads the early week-8 change (a shrinkage surrogate); the
non-default k_g link reads the post-nadir growth-rate constant (the tail-sensitive,
more-prognostic Stein/Bruno quantity). Switching the metric can re-rank — even
invert — which model looks better, completing the v0.24 finding that a week-8
surrogate is nearly blind to the resistance-model choice.
"""

import numpy as np
import onkos

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
KG_LINK = "survival_link.nsclc_os_growth_rate"
CLARET = "resistance.claret_2009.tgi"
TWO_POP = "resistance.nsclc_first_line.two_population"
NORTON = "drug_effect.norton_simon.nsclc"
WANG = "tgi_metrics.wang_2009.biexponential"
T = np.linspace(0.0, 312.0, 625)


def _os(rid, link=None):
    return onkos.simulate(onkos.load(), rid, context=NSCLC, drug_effect=1.0, t=T,
                          survival_link=link).median_os


def test_kg_link_is_non_default_and_metric_configured():
    ds = onkos.load()
    link = ds[KG_LINK]
    assert link.structure.get("default") is False              # opt-in, like the Cox link
    assert link.structure.get("link_metric") == "tumor_growth_rate_kg"
    # It must NOT join the default OS divergence (default links only).
    tr = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=T)
    assert "OS" in tr.survival  # the default week-8 link still drives OS


def test_metric_default_is_backward_compatible():
    """A link with no declared link_metric still reads the week-8 change, so the
    default OS curves are unchanged by the metric-configurability."""
    ds = onkos.load()
    default_os = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=T).median_os
    explicit_week8 = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=T,
                                    survival_link="survival_link.nsclc_os_week8").median_os
    assert np.isclose(default_os, explicit_week8)


def test_growth_rate_link_gives_a_different_os_than_week8():
    """The two metrics are not the same survival model — k_g moves the prediction."""
    assert not np.isclose(_os(CLARET), _os(CLARET, KG_LINK), rtol=0.05)


def test_metric_choice_inverts_the_resistance_model_ranking():
    """Under the week-8 surrogate the mechanistic (deeper early shrinkage) model looks
    BETTER than the phenomenological one; under the k_g metric (faster regrowth) it
    looks WORSE. Which survival metric you pick flips the resistance-model ranking —
    the v0.24 finding made consequential."""
    assert _os(TWO_POP) > _os(CLARET)                  # week-8: two-population wins
    assert _os(TWO_POP, KG_LINK) < _os(CLARET, KG_LINK)  # k_g: phenomenological wins


def test_complete_responder_is_penalized_by_week8_but_rewarded_by_kg():
    """Norton-Simon eradicates (slow early shrinkage, no regrowth). The week-8
    surrogate undervalues it; the k_g metric (no regrowth -> baseline hazard)
    correctly makes it the longest survivor."""
    assert _os(NORTON) < _os(CLARET)                   # week-8: undervalued
    assert _os(NORTON, KG_LINK) > _os(CLARET, KG_LINK)  # k_g: rewarded


def test_undefined_growth_rate_maps_to_baseline_hazard():
    """A model with no regrowth has k_g = nan, which maps to the no-effect covariate
    (x=0) — a finite, best-case (baseline) survival, never a nan curve."""
    ds = onkos.load()
    m = onkos.simulate(ds, NORTON, context=NSCLC, drug_effect=1.0, t=T).metrics
    assert not np.isfinite(m["tumor_growth_rate_kg"])  # eradication: no regrowth
    g = onkos.simulate(ds, NORTON, context=NSCLC, drug_effect=1.0, t=T, survival_link=KG_LINK)
    assert np.all(np.isfinite(g.os_curve))
    assert g.median_os is not None


def test_growth_rate_link_transports_and_tiers_like_any_link():
    """Used outside its validated context the k_g link floors the composed tier to D."""
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.crc_first_line.claret",
                        context={"tumor_type": "CRC", "line": "first"}, drug_effect=1.0, t=T,
                        survival_link=KG_LINK)
    assert tr.tier == "D"
    assert any("outside validated" in w for w in tr.warnings)
