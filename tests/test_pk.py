"""PK bridge: exposure metrics + the full PK -> exposure -> TGI -> survival chain."""

import numpy as np
import onkos
from onkos import pk

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
PKARGS = dict(ka=0.5, ke=0.05, v=10.0, f=0.8)


def test_c_avg_matches_closed_form():
    m = pk.steady_state_metrics(dose=300, tau=24, **PKARGS)
    cl = PKARGS["ke"] * PKARGS["v"]
    assert m["c_avg"] == PKARGS["f"] * 300 / (cl * 24)
    assert m["auc_tau"] == PKARGS["f"] * 300 / cl


def test_c_avg_is_dose_linear():
    a = pk.steady_state_metrics(dose=100, tau=24, **PKARGS)["c_avg"]
    b = pk.steady_state_metrics(dose=300, tau=24, **PKARGS)["c_avg"]
    assert abs(b / a - 3.0) < 1e-9


def test_metrics_ordering_and_accumulation():
    m = pk.steady_state_metrics(dose=300, tau=12, **PKARGS)
    assert m["c_min"] < m["c_avg"] < m["c_max"]
    assert m["accumulation"] > 1.0  # multiple dosing accumulates with this half-life


def test_single_dose_profile_rises_then_falls():
    t = np.linspace(0, 200, 400)
    c = pk.concentration_profile(300, 1e9, 1, **PKARGS, t=t)  # one dose (tau huge)
    peak = int(np.argmax(c))
    assert 0 < peak < len(t) - 1
    # Bateman peak time t_max = ln(ka/ke)/(ka-ke)
    tmax = np.log(PKARGS["ka"] / PKARGS["ke"]) / (PKARGS["ka"] - PKARGS["ke"])
    assert abs(t[peak] - tmax) < 1.0


def test_multiple_dose_accumulates():
    t = np.linspace(0, 240, 600)
    one = pk.concentration_profile(300, 24, 1, **PKARGS, t=t)
    many = pk.concentration_profile(300, 24, 10, **PKARGS, t=t)
    assert np.max(many) > np.max(one)


def test_from_profile_ingests_and_interpolates():
    t = np.linspace(0, 100, 201)
    c = pk.from_profile([0, 10, 100], [0.0, 50.0, 10.0], t)
    assert c[0] == 0.0
    assert abs(c[np.argmin(np.abs(t - 10))] - 50.0) < 1e-6
    # held at endpoints outside the supplied range
    assert pk.from_profile([10, 20], [5.0, 7.0], np.array([0.0]))[0] == 5.0


def test_full_chain_dose_response():
    """Higher dose -> higher C_avg -> deeper response -> longer OS (the go/no-go chain)."""
    ds = onkos.load()
    er = "exposure_response.emax_generic"
    out = {}
    for dose in (1000, 2000, 4000):
        cavg = pk.steady_state_metrics(dose=dose, tau=24, **PKARGS)["c_avg"]
        tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC,
                            exposure=cavg, exposure_response=er)
        out[dose] = (tr.metrics["depth_of_response"], tr.median_os)
    assert out[2000][0] > out[1000][0]      # deeper response
    assert out[4000][0] > out[2000][0]
    assert out[4000][1] > out[1000][1]      # longer median OS


def test_time_varying_ingested_profile_drives_tgi():
    ds = onkos.load()
    t = np.linspace(0, 104, 209)
    profile = pk.from_profile([0, 8, 52, 104], [0, 300, 200, 120], t)
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC,
                        exposure=profile, exposure_response="exposure_response.emax_generic", t=t)
    assert tr.tumor_size.shape == t.shape
    assert tr.metrics["depth_of_response"] > 0
