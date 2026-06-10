"""RECIST best-response & ORR landmark suite, plus the ORR -> OS conditional-surrogacy
result (research spec recist-orr-surrogate §5).

The classifier is RECIST 1.1 arithmetic over a trajectory; the rates are population
fractions over the IIV ensemble; the surrogate discordance is a statement about models.
"""

import numpy as np
import onkos
from onkos.response import (
    best_response,
    objective_response_rate,
    response_vs_survival,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
KG_LINK = "survival_link.nsclc_os_growth_rate"
T = np.linspace(0.0, 156.0, 313)


# --- RECIST best-response classification -----------------------------------


def test_pr_boundary_is_thirty_percent():
    t = np.linspace(0, 52, 105)
    assert best_response(t, np.linspace(100, 70, 105)) == "PR"      # exactly -30%
    assert best_response(t, np.linspace(100, 71, 105)) == "SD"      # -29%, no progression


def test_cr_threshold_is_near_complete():
    t = np.linspace(0, 52, 105)
    assert best_response(t, np.linspace(100, 4, 105)) == "CR"       # 96% shrinkage
    assert best_response(t, np.linspace(100, 6, 105)) == "PR"       # 94% -> PR, not CR


def test_pd_requires_no_pr_and_regrowth_from_nadir():
    t = np.linspace(0, 52, 105)
    # shrink 10% to nadir then regrow well past +20% of nadir, never reaching PR.
    v = np.concatenate([np.linspace(100, 90, 50), np.linspace(90, 140, 55)])
    assert best_response(t, v) == "PD"


def test_pr_that_later_regrows_is_still_pr():
    """Best overall response is the best timepoint: a PR that progresses stays PR."""
    t = np.linspace(0, 104, 209)
    v = np.concatenate([np.linspace(100, 60, 100), np.linspace(60, 130, 109)])  # PR at nadir
    assert best_response(t, v) == "PR"


def test_flat_trajectory_is_stable_disease():
    t = np.linspace(0, 52, 105)
    assert best_response(t, np.full(105, 100.0)) == "SD"


def test_monotone_growth_is_progressive_disease():
    t = np.linspace(0, 52, 105)
    assert best_response(t, np.linspace(100, 200, 105)) == "PD"


# --- population rates -------------------------------------------------------


def test_rates_form_a_simplex_and_order():
    rr = objective_response_rate(onkos.load(), "resistance.claret_2009.tgi",
                                 context=NSCLC, n=200)
    assert np.isclose(sum(rr.distribution.values()), 1.0)
    assert 0.0 <= rr.orr <= rr.dcr <= 1.0


def test_orr_is_monotone_in_drug_effect():
    ds = onkos.load()
    low = objective_response_rate(ds, "resistance.claret_2009.tgi", context=NSCLC,
                                  drug_effect=0.5, n=200)
    high = objective_response_rate(ds, "resistance.claret_2009.tgi", context=NSCLC,
                                   drug_effect=1.5, n=200)
    assert high.orr >= low.orr - 1e-9


def test_zero_iiv_gives_degenerate_orr():
    """A model with no reported IIV yields identical samples ⇒ ORR is 0 or 1."""
    ds = onkos.load()
    rr = objective_response_rate(ds, "tgi_metrics.wang_2009.biexponential", context=NSCLC,
                                 n=50)
    # wang biexp carries IIV; instead verify the degeneracy property on a no-IIV growth law.
    rr0 = objective_response_rate(ds, "growth_laws.exponential", context=NSCLC, n=50)
    assert rr0.orr in (0.0, 1.0)
    assert rr.distribution["CR"] + rr.distribution["PR"] == rr.orr


def test_result_carries_clinical_use():
    rr = objective_response_rate(onkos.load(), "resistance.claret_2009.tgi",
                                 context=NSCLC, n=80)
    d = rr.to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True and "PROHIBITED" in d["onkos:clinicalUse"]
    assert "distribution" in d and "median_os_weeks" in d


def test_tier_passes_through_and_floors_on_transport():
    ds = onkos.load()
    rr = objective_response_rate(ds, "resistance.claret_2009.tgi", context=NSCLC, n=60)
    assert rr.tier == "C"
    out = objective_response_rate(ds, "resistance.claret_2009.tgi",
                                  context={"tumor_type": "melanoma", "line": "first"}, n=60)
    assert out.tier == "D"


# --- the ORR -> OS conditional surrogate -----------------------------------


def test_orr_predicts_os_under_week8_but_not_under_kg():
    """The headline: ORR ranks OS perfectly under the shrinkage-based week-8 link, and
    badly under the tail-sensitive k_g link — surrogacy is conditional on the survival
    mechanism."""
    ds = onkos.load()
    week8 = response_vs_survival(ds, context=NSCLC, t=T, n=250)
    assert week8.orr_predicts_os                     # 0 discordant pairs
    assert week8.discordant_fraction == 0.0

    kg = response_vs_survival(ds, context=NSCLC, survival_link=KG_LINK, t=T, n=250)
    assert not kg.orr_predicts_os                    # ORR mis-ranks tail-driven survival
    assert kg.discordant_fraction > 0.0


def test_surrogate_result_carries_clinical_use_and_pairs():
    ds = onkos.load()
    d = response_vs_survival(ds, context=NSCLC, survival_link=KG_LINK, t=T, n=200).to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert d["total_pairs"] >= 1 and 0.0 <= d["discordant_fraction"] <= 1.0
    assert len(d["rows"]) >= 3


# --- cross-context generalization (v0.29) ----------------------------------


def test_orr_surrogate_generalizes_across_solid_tumor_contexts():
    """v0.29 breadth: the conditional ORR -> OS surrogacy is not an NSCLC artifact. In
    every first-line solid-tumor context, ORR ranks OS faithfully under the shrinkage-based
    week-8 link but mis-ranks it under the tail-sensitive k_g link — because each context
    now carries a mechanistic two-population model (high ORR, fast resistant regrowth)."""
    ds = onkos.load()
    t = np.linspace(0.0, 312.0, 625)
    for tt in ("NSCLC", "breast", "CRC", "HCC", "melanoma"):
        ctx = {"tumor_type": tt, "line": "first"}
        kg_link = f"survival_link.{tt.lower()}_os_growth_rate"
        week8 = response_vs_survival(ds, context=ctx, t=t, n=200)
        kg = response_vs_survival(ds, context=ctx, survival_link=kg_link, t=t, n=200)
        assert week8.orr_predicts_os, f"{tt}: expected ORR concordant under week-8"
        assert kg.discordant_fraction > 0.0, f"{tt}: expected ORR discordant under k_g"
