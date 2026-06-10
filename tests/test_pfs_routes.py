"""PFS two-route landmark suite (research spec pfs-endpoint §5): a PFS number has two
legitimate routes — the statistical week-8-keyed survival link and the mechanistic RECIST
time-to-progression read off the tumor trajectory — and for shrink-then-regrow resistance
dynamics they invert the model ranking, in every solid-tumor context.
"""

import numpy as np
import onkos
from onkos.response import (
    pfs_route_divergence,
    progression_free_survival,
    time_to_progression,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CONTEXTS = [
    {"tumor_type": "NSCLC", "line": "first"},
    {"tumor_type": "breast", "line": "first"},
    {"tumor_type": "CRC", "line": "first"},
    {"tumor_type": "HCC", "line": "first"},
    {"tumor_type": "melanoma", "line": "first"},
]
CLARET = "resistance.claret_2009.tgi"
TWO_POP = "resistance.nsclc_first_line.two_population"


# --- time_to_progression: the mechanistic route on one trajectory -----------


def test_closed_form_ttp():
    """A shrink-then-regrow trajectory crosses +20% of its nadir at a known time."""
    t = np.linspace(0, 104, 209)
    # shrink to 50 by ~week 26, hold, then regrow; nadir = 50, PD level = 60.
    v = np.concatenate([np.linspace(100, 50, 53), np.full(80, 50.0), np.linspace(50, 90, 76)])
    ttp = time_to_progression(t, v)
    # first POST-NADIR time v >= 1.2 * nadir(=50) = 60, on the regrowth limb.
    nadir_i = int(np.argmin(v))
    k = next(i for i in range(nadir_i, len(t)) if v[i] >= 60.0)
    assert np.isclose(ttp, t[k])


def test_running_nadir_rule():
    """Progression is measured against the running nadir, not baseline: a deep-shrink-then-
    regrow trajectory progresses when it rises 20% above its *nadir*, far below baseline."""
    t = np.linspace(0, 104, 209)
    v = np.concatenate([np.linspace(100, 40, 105), np.linspace(40, 60, 104)])  # nadir 40 -> 60 = +50%
    ttp = time_to_progression(t, v)
    assert np.isfinite(ttp)
    # it progresses well before the SLD ever returns to baseline (100).
    assert v[next(i for i in range(len(t)) if t[i] >= ttp)] < 100.0


def test_monotone_growth_progresses_against_baseline():
    t = np.linspace(0, 104, 209)
    v = np.linspace(100, 200, 209)  # nadir = baseline at t=0; +20% at v=120
    ttp = time_to_progression(t, v)
    k = next(i for i in range(len(t)) if v[i] >= 120.0)
    assert np.isclose(ttp, t[k])


def test_durable_non_progressor_is_censored():
    """A trajectory that shrinks and never regrows past +20% of its nadir is censored (nan)."""
    t = np.linspace(0, 104, 209)
    v = np.linspace(100, 40, 209)  # monotone shrink, never regrows
    assert np.isnan(time_to_progression(t, v))


# --- progression_free_survival: both routes over the ensemble ---------------


def test_both_routes_present_for_a_pfs_context():
    ds = onkos.load()
    pf = progression_free_survival(ds, CLARET, context=NSCLC, n=120)
    assert pf.has_pfs_link
    assert pf.median_ttp_weeks is not None and np.isfinite(pf.median_ttp_weeks)
    assert pf.median_pfs_link_weeks is not None and np.isfinite(pf.median_pfs_link_weeks)
    assert pf.route_ratio is not None


def test_landmark_rate_is_a_probability():
    ds = onkos.load()
    pf = progression_free_survival(ds, CLARET, context=NSCLC, n=120)
    assert 0.0 <= pf.mechanistic_pfs_rate <= 1.0
    assert 0.0 <= pf.ttp_censored_fraction <= 1.0


def test_clinical_use_flag_on_payload():
    ds = onkos.load()
    pf = progression_free_survival(ds, CLARET, context=NSCLC, n=60)
    d = pf.to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "PROHIBITED" in d["onkos:clinicalUse"]


def test_out_of_context_floors_tier():
    """Transporting a model outside its validated context floors the tier to D, exactly as
    for OS/ORR (worst-input-wins rides through the PFS endpoint)."""
    ds = onkos.load()
    pf = progression_free_survival(
        ds, CLARET, context={"tumor_type": "melanoma", "line": "first"}, n=60
    )
    assert pf.tier == "D"


# --- pfs_route_divergence: the route is a model-selection axis ---------------


def test_nsclc_route_inversion():
    """The two-population model has a SHORTER mechanistic PFS than Claret but a LONGER (or
    equal) statistical PFS — the routes disagree on which model looks better."""
    ds = onkos.load()
    div = pfs_route_divergence(ds, context=NSCLC, n=200)
    by_id = {r["record_id"]: r for r in div.rows}
    claret, twopop = by_id[CLARET], by_id[TWO_POP]
    assert twopop["median_ttp_weeks"] < claret["median_ttp_weeks"]          # mech: two-pop shorter
    assert twopop["median_pfs_link_weeks"] >= claret["median_pfs_link_weeks"]  # stat: two-pop longer
    assert div.discordant_pairs >= 1
    assert not div.routes_agree


def test_dataset_wide_route_discordance():
    """The route disagreement reproduces in every solid-tumor context — not an NSCLC
    artifact (every context has both a PFS link and a two-population model)."""
    ds = onkos.load()
    for ctx in CONTEXTS:
        div = pfs_route_divergence(ds, context=ctx, n=150)
        assert div.discordant_pairs >= 1, ctx
        assert 0.0 <= div.discordant_fraction <= 1.0


def test_divergence_rows_carry_tier_and_both_routes():
    ds = onkos.load()
    div = pfs_route_divergence(ds, context=NSCLC, n=80)
    assert div.total_pairs >= 1
    for r in div.rows:
        assert r["tier"] in {"A", "B", "C", "D"}
        assert "median_ttp_weeks" in r and "median_pfs_link_weeks" in r
