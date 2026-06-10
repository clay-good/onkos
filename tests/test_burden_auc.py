"""Integrated tumor burden — a third TGI→OS bridge metric (v0.33.0).

The default OS link reads the early week-8 change (depth-only, blind to the tail);
the non-default k_g link reads the post-nadir regrowth slope (tail-only, blind to
depth). The burden-AUC metric `log_burden_auc` is the time-averaged log relative
tumor size over the horizon (the AUC of the log-size curve) — the one summary that
integrates BOTH depth and tail. It re-ranks the model set a third way, and it
repairs a depth-blind pathology of the pure-tail k_g metric.

Closed-form landmarks pin the metric; the OS-ranking landmarks pin the finding.
"""

import numpy as np
import onkos
from onkos.metrics import _BURDEN_FLOOR, extract_tgi_metrics

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
BURDEN_LINK = "survival_link.nsclc_os_burden_auc"
KG_LINK = "survival_link.nsclc_os_growth_rate"
WEEK8_LINK = "survival_link.nsclc_os_week8"
CLARET = "resistance.claret_2009.tgi"
TWO_POP = "resistance.nsclc_first_line.two_population"
NORTON = "drug_effect.norton_simon.nsclc"
WANG = "tgi_metrics.wang_2009.biexponential"
ACQUIRED = "resistance.nsclc_first_line.acquired"
T = np.linspace(0.0, 260.0, 521)


def _os(rid, link=None):
    return onkos.simulate(onkos.load(), rid, context=NSCLC, drug_effect=1.0, t=T,
                          survival_link=link).median_os


# ---- closed-form metric landmarks -----------------------------------------

def test_baseline_trajectory_has_zero_burden():
    """A tumor held at baseline (v ≡ y0) has log_burden_auc = 0 exactly — the same
    zero-point convention as week8_relative_change, so a no-effect tumor maps to the
    baseline hazard."""
    t = np.linspace(0.0, 100.0, 201)
    v = np.full_like(t, 5.0)
    m = extract_tgi_metrics(t, v, y0=5.0)
    assert np.isclose(m["log_burden_auc"], 0.0, atol=1e-12)


def test_constant_size_burden_is_log_of_ratio():
    """A trajectory held at c·y0 (constant) has log_burden_auc = log(c) exactly."""
    t = np.linspace(0.0, 100.0, 201)
    for c in (0.5, 0.1, 2.0):
        m = extract_tgi_metrics(t, np.full_like(t, c * 5.0), y0=5.0)
        assert np.isclose(m["log_burden_auc"], np.log(c), atol=1e-12)


def test_burden_is_monotone_in_size():
    """Larger constant size ⇒ larger (less negative) integrated burden."""
    t = np.linspace(0.0, 100.0, 201)
    small = extract_tgi_metrics(t, np.full_like(t, 1.0), y0=5.0)["log_burden_auc"]
    large = extract_tgi_metrics(t, np.full_like(t, 4.0), y0=5.0)["log_burden_auc"]
    assert large > small


def test_eradication_is_floored_and_finite():
    """A trajectory reaching v = 0 yields a finite burden bounded below by log(floor),
    not −∞ — the integral is stable under complete response."""
    t = np.linspace(0.0, 100.0, 201)
    v = np.zeros_like(t)            # immediate, total eradication
    m = extract_tgi_metrics(t, v, y0=5.0)
    assert np.isfinite(m["log_burden_auc"])
    assert np.isclose(m["log_burden_auc"], np.log(_BURDEN_FLOOR), atol=1e-12)
    # A real (gradual) complete responder also stays finite and ≥ the floor.
    mn = onkos.simulate(onkos.load(), NORTON, context=NSCLC, drug_effect=1.0, t=T).metrics
    assert np.isfinite(mn["log_burden_auc"])
    assert mn["log_burden_auc"] >= np.log(_BURDEN_FLOOR) - 1e-9


def test_burden_grows_with_horizon_over_a_regrowing_tail():
    """For a monotonically growing tumor the integrated burden increases with the
    horizon — the documented cumulative-summary property (more time large ⇒ more
    burden)."""
    g = 0.02
    short_t = np.linspace(0.0, 100.0, 201)
    long_t = np.linspace(0.0, 300.0, 601)
    short = extract_tgi_metrics(short_t, 5.0 * np.exp(g * short_t), y0=5.0)["log_burden_auc"]
    long = extract_tgi_metrics(long_t, 5.0 * np.exp(g * long_t), y0=5.0)["log_burden_auc"]
    assert long > short
    assert np.isclose(short, g * 100.0 / 2.0, rtol=1e-3)   # ∫g·t/T = g·T/2


# ---- the bridge-metric finding --------------------------------------------

def test_burden_link_is_non_default_and_metric_configured():
    ds = onkos.load()
    link = ds[BURDEN_LINK]
    assert link.structure.get("default") is False                 # opt-in, like k_g/Cox
    assert link.structure.get("link_metric") == "log_burden_auc"


def test_default_view_is_unchanged_by_the_new_link():
    """The non-default burden link does not join the default OS divergence; the default
    week-8 OS prediction is identical to selecting the week-8 link explicitly."""
    ds = onkos.load()
    default_os = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=T).median_os
    week8_os = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=T,
                              survival_link=WEEK8_LINK).median_os
    assert np.isclose(default_os, week8_os)


def test_burden_metric_is_tail_sensitive_where_week8_is_blind():
    """Matched as near-equal week-8 responders, the deep-then-regrow model (two-pop)
    carries MORE integrated burden than the eradicating complete responder — the metric
    sees the tail the week-8 change misses."""
    ds = onkos.load()
    tp = onkos.simulate(ds, TWO_POP, context=NSCLC, drug_effect=1.0, t=T).metrics
    ns = onkos.simulate(ds, NORTON, context=NSCLC, drug_effect=1.0, t=T).metrics
    assert tp["log_burden_auc"] > ns["log_burden_auc"]


def test_burden_metric_is_depth_sensitive_where_kg_is_blind():
    """The minimal responder (Wang: nadir ~75% of baseline) carries MORE burden than a
    deep responder, so it is ranked LAST by burden — even though k_g ranks it high on a
    slow regrowth slope. Burden repairs the pure-tail metric's depth-blindness."""
    ds = onkos.load()
    wang = onkos.simulate(ds, WANG, context=NSCLC, drug_effect=1.0, t=T).metrics
    claret = onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=1.0, t=T).metrics
    assert wang["log_burden_auc"] > claret["log_burden_auc"]      # Wang never got small
    # k_g, by contrast, ranks Wang ABOVE the deep-but-doomed resistance models:
    assert _os(WANG, KG_LINK) > _os(TWO_POP, KG_LINK)
    # burden puts Wang last among them:
    assert _os(WANG, BURDEN_LINK) < _os(TWO_POP, BURDEN_LINK)


def test_burden_link_gives_a_third_distinct_ranking():
    """The burden-AUC OS ranking of the NSCLC model set differs from BOTH the week-8 and
    the k_g ranking — which bridge metric you pick remains a live model-selection axis."""
    models = [CLARET, TWO_POP, NORTON, WANG, ACQUIRED]

    def order(link):
        scored = sorted(models, key=lambda r: _os(r, link) or -1.0, reverse=True)
        return tuple(scored)

    week8, kg, burden = order(WEEK8_LINK), order(KG_LINK), order(BURDEN_LINK)
    assert burden != week8
    assert burden != kg
    # The complete responder inverts from mid-pack (week-8) to first (burden):
    assert week8.index(NORTON) > 0
    assert burden.index(NORTON) == 0


def test_burden_link_is_eligible_and_transports_like_any_link():
    """It joins the model-selection-budget V_link factor (4 OS links for NSCLC/first),
    and used outside its validated context it floors the composed tier to D."""
    from onkos.budget import eligible_survival_links
    ds = onkos.load()
    assert BURDEN_LINK in eligible_survival_links(ds, NSCLC, "OS")
    assert len(eligible_survival_links(ds, NSCLC, "OS")) == 4     # week8, cox, k_g, burden
    tr = onkos.simulate(ds, "resistance.crc_first_line.claret",
                        context={"tumor_type": "CRC", "line": "first"}, drug_effect=1.0, t=T,
                        survival_link=BURDEN_LINK)
    assert tr.tier == "D"
    assert any("outside validated" in w for w in tr.warnings)
