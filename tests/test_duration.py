"""Duration-of-response landmark suite (research spec duration-of-response §5):
DoR is the durability dimension ORR omits, and the ORR-surrogate failure is a
durability failure (the highest-ORR model is the shortest-DoR one).
"""

import numpy as np
import onkos
from onkos.response import best_response, objective_response_rate, response_episode

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
KG_LINK = "survival_link.nsclc_os_growth_rate"
T = np.linspace(0.0, 208.0, 417)


# --- response_episode: category + DoR from one trajectory ------------------


def test_episode_category_matches_best_response():
    t = np.linspace(0, 104, 209)
    for v in (
        np.linspace(100, 60, 209),                                   # PR
        np.full(209, 100.0),                                         # SD
        np.linspace(100, 200, 209),                                  # PD
        np.concatenate([np.linspace(100, 60, 100), np.linspace(60, 140, 109)]),  # PR then PD
    ):
        assert response_episode(t, v)[0] == best_response(t, v)


def test_non_responder_has_nan_dor():
    t = np.linspace(0, 104, 209)
    assert np.isnan(response_episode(t, np.full(209, 100.0))[1])        # SD
    assert np.isnan(response_episode(t, np.linspace(100, 200, 209))[1])  # PD


def test_closed_form_dor():
    """A PR that onsets early and progresses later: DoR = t_PD − t_PR."""
    t = np.linspace(0, 104, 209)
    # shrink to 50% by week ~26, hold, then regrow past +20% of nadir by the end.
    v = np.concatenate([np.linspace(100, 50, 53), np.full(80, 50.0), np.linspace(50, 90, 76)])
    cat, dor = response_episode(t, v)
    assert cat == "PR"
    onset = next(t[k] for k in range(len(t)) if v[k] <= 70.0)
    nadir = v.min()
    prog = next(t[k] for k in range(int(np.argmin(v)), len(t)) if v[k] >= 1.2 * nadir)
    assert np.isclose(dor, prog - onset)


def test_durable_response_is_censored():
    """A response that never regrows ⇒ DoR is nan (right-censored), not zero."""
    t = np.linspace(0, 104, 209)
    v = np.concatenate([np.linspace(100, 40, 100), np.full(109, 40.0)])  # PR, no progression
    cat, dor = response_episode(t, v)
    assert cat == "PR" and np.isnan(dor)


def test_slower_regrowth_gives_longer_dor():
    ds = onkos.load()
    # two-population kgr (resistant growth) drives regrowth; a smaller kgr -> longer DoR.
    fast = objective_response_rate(ds, "resistance.nsclc_first_line.two_population",
                                   context=NSCLC, t=T, n=300)
    slow = objective_response_rate(ds, "resistance.claret_2009.tgi", context=NSCLC, t=T, n=300)
    assert slow.median_dor_weeks > fast.median_dor_weeks


# --- population DoR honesty -------------------------------------------------


def test_dor_fields_are_consistent():
    rr = objective_response_rate(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC,
                                 t=T, n=300)
    assert 0.0 <= rr.dor_censored_fraction <= 1.0
    n_resp = round((rr.distribution["CR"] + rr.distribution["PR"]) * rr.n)
    assert abs(rr.n_responders - n_resp) <= 1               # responders == CR+PR count
    assert rr.median_dor_weeks is None or rr.median_dor_weeks > 0


def test_heavy_censoring_is_flagged():
    """The eradicating model produces many durable (censored) responders."""
    rr = objective_response_rate(onkos.load(), "drug_effect.norton_simon.nsclc",
                                 context=NSCLC, t=T, n=300)
    if rr.dor_censored_fraction >= 0.5:
        assert any("dor_heavily_censored" in w for w in rr.warnings)


def test_dor_in_dict_and_tier_passthrough():
    ds = onkos.load()
    d = objective_response_rate(ds, "resistance.claret_2009.tgi", context=NSCLC, t=T, n=120).to_dict()
    assert "median_dor_weeks" in d and "dor_censored_fraction" in d and "n_responders" in d
    out = objective_response_rate(ds, "resistance.claret_2009.tgi",
                                  context={"tumor_type": "melanoma", "line": "first"}, t=T, n=80)
    assert out.tier == "D"


# --- depth != durability, and the surrogate mechanism ----------------------


def test_highest_orr_model_has_shorter_dor_than_a_lower_orr_model():
    """Depth is not durability: the highest-ORR NSCLC model is *not* the most durable."""
    ds = onkos.load()
    two_pop = objective_response_rate(ds, "resistance.nsclc_first_line.two_population",
                                      context=NSCLC, t=T, n=400)
    claret = objective_response_rate(ds, "resistance.claret_2009.tgi", context=NSCLC, t=T, n=400)
    assert two_pop.orr >= claret.orr                       # higher (or equal) breadth
    assert two_pop.median_dor_weeks < claret.median_dor_weeks  # but shorter durability


def test_durability_tracks_survival_where_breadth_does_not():
    """Under the tail-sensitive k_g link, the highest-ORR model has the SHORTEST OS (the
    v0.27 discordance) — yet the longest-OS model is MORE durable (longer DoR) than that
    highest-ORR model. Durability tracks survival where breadth (ORR) inverts it: DoR is
    the mechanism of the surrogate failure."""
    ds = onkos.load()
    rs = onkos.response_vs_survival(ds, context=NSCLC, survival_link=KG_LINK, t=T, n=300)
    rows = [r for r in rs.rows if r["median_dor_weeks"] is not None]
    top_orr = max(rows, key=lambda r: r["orr"])
    longest_os = max(rows, key=lambda r: r["median_os_weeks"])
    shortest_os = min(rows, key=lambda r: r["median_os_weeks"])
    # The most-responsive model is the worst survivor (ORR inverts OS) ...
    assert top_orr["record_id"] == shortest_os["record_id"]
    # ... but the best survivor is more durable than it (DoR aligns with OS).
    assert longest_os["median_dor_weeks"] > top_orr["median_dor_weeks"]
