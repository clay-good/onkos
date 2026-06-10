"""Early-surrogate readout timing — *when* you read the surrogate is a model-selection
axis (v0.37.0).

ctDNA molecular response is modeled as proportional to tumor burden, so the modeled
distinction from a RECIST-size readout is purely the landmark TIME. An earlier landmark
sits closer to the nadir, before resistant regrowth, so it over-rewards deep-but-doomed
early responders; the discordance with a tail-aware durable-benefit ranking falls as the
landmark moves later.

Closed-form landmarks pin the pure core; the binding landmarks pin the timing finding.
"""

import numpy as np
import onkos
from onkos.early_surrogate import (
    SurrogateTiming,
    discordant_pairs,
    landmark_response,
    surrogate_timing_fidelity,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


# ---- pure-core landmarks ---------------------------------------------------

def test_landmark_response_recovers_week8_metric():
    """landmark_response at week 8 equals the curated week8_relative_change metric — the
    generalization recovers the special case the rest of Onkos uses."""
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0)
    assert np.isclose(landmark_response(tr.t, tr.tumor_size, 8.0),
                      tr.metrics["week8_relative_change"], rtol=1e-9)


def test_landmark_response_is_record_free_and_signed():
    """Works on a synthetic trajectory; 0 at baseline, negative under shrinkage."""
    t = np.linspace(0.0, 100.0, 201)
    flat = np.full_like(t, 5.0)
    assert landmark_response(t, flat, 50.0) == 0.0
    shrink = 5.0 * np.exp(-0.05 * t)
    assert landmark_response(t, shrink, 20.0) < 0.0


def test_landmark_response_tracks_deepening_then_regrowth():
    """For a shrink-then-regrow trajectory the response is deepest near the nadir and less
    deep at a late landmark once regrowth sets in — the reason an early readout flatters
    a doomed responder."""
    t = np.linspace(0.0, 200.0, 401)
    v = 5.0 * (0.05 + np.exp(-0.3 * t) + 0.02 * t)  # deep dip then linear regrowth
    nadir_week = float(t[np.argmin(v)])
    early = landmark_response(t, v, nadir_week)
    late = landmark_response(t, v, 180.0)
    assert early < late  # earlier landmark shows the deeper (more favorable) response


def test_discordant_pairs_counts_inversions():
    assert discordant_pairs(["a", "b", "c"], ["a", "b", "c"]) == 0
    assert discordant_pairs(["a", "b", "c"], ["c", "b", "a"]) == 3
    assert discordant_pairs(["a", "b", "c"], ["b", "a", "c"]) == 1


# ---- the timing finding ----------------------------------------------------

def test_fidelity_improves_with_a_later_landmark():
    """The core finding: reading the surrogate later agrees BETTER with the durable-benefit
    ranking than reading it early — earliness trades against fidelity."""
    st = surrogate_timing_fidelity(onkos.load(), context=NSCLC)
    assert isinstance(st, SurrogateTiming)
    assert st.earliest_discordance > st.latest_discordance


def test_early_readout_is_substantially_discordant():
    """At the ctDNA-era early readout the surrogate ranking is mostly inverted relative to
    durable benefit (it over-rewards deep-but-doomed early responders)."""
    st = surrogate_timing_fidelity(onkos.load(), context=NSCLC)
    early = st.discordance_at(2.0)
    assert early >= 0.7 * st.total_pairs  # the earliest landmark is badly discordant


def test_discordance_is_weakly_monotone_in_landmark():
    """Discordance never INCREASES as the landmark moves later (it falls or holds) — the
    monotone earliness-fidelity trade-off."""
    st = surrogate_timing_fidelity(onkos.load(), context=NSCLC)
    d = [r["discordant_pairs"] for r in st.rows]
    assert all(d[i + 1] <= d[i] for i in range(len(d) - 1))


def test_durable_benefit_reference_demotes_the_fast_doomed_models():
    """The tail-aware reference ranks the mechanistic-resistance models (fast deep response,
    fast regrowth) near the BOTTOM — exactly the models an early readout puts on top."""
    st = surrogate_timing_fidelity(onkos.load(), context=NSCLC)
    ref = st.reference_ranking
    assert ref.index("resistance.nsclc_first_line.acquired") >= len(ref) - 2
    assert ref.index("resistance.nsclc_first_line.two_population") >= len(ref) - 3
    # and the early-landmark ranking puts a fast-doomed model on top:
    early_top = st.rows[0]["ranking"][0]
    assert early_top in ("resistance.nsclc_first_line.acquired",
                         "resistance.nsclc_first_line.two_population")


def test_reproduces_across_contexts():
    """The earliness-fidelity trade-off is not an NSCLC artifact: it holds in the other
    solid-tumor contexts too (each has a two-population model + a k_g link since v0.29)."""
    ds = onkos.load()
    for tt in ("breast", "CRC", "HCC", "melanoma"):
        st = surrogate_timing_fidelity(ds, context={"tumor_type": tt, "line": "first"})
        assert st.earliest_discordance >= st.latest_discordance
        assert st.total_pairs >= 1


def test_inherits_tier_and_carries_clinical_use():
    st = surrogate_timing_fidelity(onkos.load(), context=NSCLC)
    assert st.tier in ("A", "B", "C", "D")
    import json

    d = json.loads(st.to_json())
    assert d["onkos:clinicalUse"].startswith("PROHIBITED")
    assert d["NOT_FOR_CLINICAL_USE"] is True
