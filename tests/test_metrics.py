"""TGI-metric extraction (Stein/Bruno panel) — model-agnostic, parameter-recovering."""

import numpy as np
import onkos
from onkos.metrics import extract_tgi_metrics

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
T = np.linspace(0, 156, 313)


def _metrics(rid, drug_effect=1.0):
    return onkos.simulate(onkos.load(), rid, context=NSCLC, drug_effect=drug_effect, t=T).metrics


def test_panel_keys_present():
    m = _metrics("resistance.claret_2009.tgi")
    for k in (
        "depth_of_response",
        "nadir_tumor_size",
        "time_to_nadir_weeks",
        "tumor_growth_rate_kg",
        "tumor_shrinkage_rate_ks",
        "time_to_growth_weeks",
        "duration_of_response_weeks",
        "week8_relative_change",
    ):
        assert k in m


def test_biexp_recovers_generating_rates():
    """The extractor recovers the biexponential kernel's generating kg and ks."""
    ds = onkos.load()
    rec = ds["tgi_metrics.wang_2009.biexponential"]
    kg_gen, ks_gen = rec["kg"].central, rec["ks"].central  # ks scaled by E=1
    m = _metrics("tgi_metrics.wang_2009.biexponential")
    assert abs(m["tumor_growth_rate_kg"] - kg_gen) / kg_gen < 0.20
    assert abs(m["tumor_shrinkage_rate_ks"] - ks_gen) / ks_gen < 0.20


def test_claret_growth_rate_recovers_kL():
    ds = onkos.load()
    kL = ds["resistance.claret_2009.tgi"]["kL"].central
    m = _metrics("resistance.claret_2009.tgi")
    assert abs(m["tumor_growth_rate_kg"] - kL) / kL < 0.15


def test_pure_growth_has_no_shrinkage_rate():
    m = _metrics("growth_laws.exponential", drug_effect=0.0)
    assert np.isnan(m["tumor_shrinkage_rate_ks"])
    assert np.isnan(m["time_to_growth_weeks"]) or m["time_to_nadir_weeks"] == 0.0
    assert abs(m["tumor_growth_rate_kg"] - 0.02) < 0.005  # the exponential kg


def test_unit_trajectory_definitions():
    # synthetic shrink-then-grow: clear nadir, response, and progression
    t = np.linspace(0, 100, 401)
    y0 = 100.0
    v = y0 * (np.exp(-0.15 * t) + np.exp(0.02 * t) - 1.0)
    m = extract_tgi_metrics(t, v, y0)
    assert 0 < m["time_to_nadir_weeks"] < 100
    assert m["depth_of_response"] > 0.30  # deep enough to be a PR
    assert np.isfinite(m["duration_of_response_weeks"])  # PR reached then progressed
    assert m["tumor_growth_rate_kg"] > 0 and m["tumor_shrinkage_rate_ks"] > 0


def test_no_response_gives_nan_duration():
    # shallow shrinkage (<30%) never reaches RECIST PR -> no duration of response
    m = _metrics("tgi_metrics.wang_2009.biexponential", drug_effect=0.5)
    assert m["depth_of_response"] < 0.30
    assert np.isnan(m["duration_of_response_weeks"])
