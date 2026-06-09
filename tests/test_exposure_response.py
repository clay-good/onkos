"""Exposure-response: transform math, ER-driven simulation, PK-driven dynamics."""

import numpy as np
import onkos
import pytest
from onkos.export.reference import KERNELS, effect
from onkos.export.registry import get_kernel, kernel_values

ER_RECORDS = [
    "exposure_response.emax_generic",
    "exposure_response.sigmoid_emax_generic",
    "exposure_response.power_generic",
    "exposure_response.dacomitinib_egfr.emax",
]


def test_er_records_present_and_bound():
    ds = onkos.load()
    for rid in ER_RECORDS:
        r = ds[rid]
        assert r.purpose == "exposure_response"
        assert get_kernel(r).kind == "exposure_response"


def test_emax_half_maximal_at_ec50():
    ds = onkos.load()
    r = ds["exposure_response.emax_generic"]
    spec, v = get_kernel(r), kernel_values(ds["exposure_response.emax_generic"])
    e = float(effect(spec, v["EC50"], v))
    assert e == pytest.approx(v["Emax"] / 2.0, rel=1e-9)


def test_emax_monotone_and_bounded():
    ds = onkos.load()
    r = ds["exposure_response.emax_generic"]
    spec, v = get_kernel(r), kernel_values(r)
    c = np.linspace(0, 5000, 50)
    e = effect(spec, c, v)
    assert np.all(np.diff(e) >= 0)          # monotone increasing
    assert np.all(e < v["Emax"])            # bounded by Emax
    assert e[0] == pytest.approx(0.0)       # zero exposure -> zero effect


def test_sigmoid_steeper_than_emax_below_ec50():
    ds = onkos.load()
    plain = ds["exposure_response.emax_generic"]
    sig = ds["exposure_response.sigmoid_emax_generic"]
    vp, vs = kernel_values(plain), kernel_values(sig)
    c = vp["EC50"] / 3.0  # well below EC50: Hill>1 suppresses effect more
    ep = float(effect(get_kernel(plain), c, vp))
    es = float(effect(get_kernel(sig), c, vs))
    assert es < ep


def test_power_no_ceiling():
    ds = onkos.load()
    r = ds["exposure_response.power_generic"]
    spec, v = get_kernel(r), kernel_values(r)
    assert float(effect(spec, 2000.0, v)) > float(effect(spec, 200.0, v))


def test_constant_exposure_drives_tgi():
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    low = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                         exposure=50.0, exposure_response="exposure_response.dacomitinib_egfr.emax")
    high = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                          exposure=400.0, exposure_response="exposure_response.dacomitinib_egfr.emax")
    # more exposure -> stronger kill -> deeper response and longer median OS
    assert high.metrics["depth_of_response"] > low.metrics["depth_of_response"]
    assert high.median_os > low.median_os


def test_time_varying_pk_uses_ode_integration():
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0, 104, 209)
    # A declining PK profile should land between a constant-high and constant-low effect.
    decaying = 300.0 * np.exp(-0.02 * t)
    er = "exposure_response.emax_generic"
    tv = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                        exposure=decaying, exposure_response=er, t=t)
    assert tv.tumor_size.shape == t.shape
    assert np.all(tv.tumor_size > 0)
    # constant at the initial concentration shrinks at least as much as the decaying profile
    const_hi = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                              exposure=300.0, exposure_response=er, t=t)
    assert const_hi.metrics["depth_of_response"] >= tv.metrics["depth_of_response"] - 1e-9


def test_out_of_context_er_floors_to_D():
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi",
                        context={"tumor_type": "NSCLC", "line": "first"},
                        exposure=200.0,
                        exposure_response="exposure_response.dacomitinib_egfr.emax")
    assert tr.tier in ("B", "C")  # all in context
    # breast context: both the TGI model and the NSCLC-only ER leave their envelope
    trd = onkos.simulate(ds, "resistance.claret_2009.tgi",
                         context={"tumor_type": "breast", "line": "first"},
                         exposure=200.0,
                         exposure_response="exposure_response.dacomitinib_egfr.emax")
    assert trd.tier == "D"
    assert any("outside validated" in w for w in trd.warnings)


def test_er_kernels_registered():
    for name in ("er_emax", "er_sigmoid_emax", "er_power"):
        assert name in KERNELS
