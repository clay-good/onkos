"""Acquired (drug-induced) resistance landmark suite (research spec acquired-resistance §4):
the resistance ORIGIN — acquired switching vs a pre-existing subclone — is a model-selection
axis. Matched on kg/kd/kgr, the two origins agree at week-8 and on the week-8 OS surrogate but
diverge in the regrowth tail (shallower nadir, earlier progression for the acquired model).
"""

import numpy as np
import onkos
from onkos.export.reference import KERNELS
from onkos.response import time_to_progression

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
ACQUIRED = "resistance.nsclc_first_line.acquired"
PREEXISTING = "resistance.nsclc_first_line.two_population"
T = np.linspace(0.0, 156.0, 313)


# --- kernel-level (pure rhs) landmarks --------------------------------------


def test_recovers_two_population_at_alpha_zero():
    """With α=0 the acquired rhs equals the pre-existing two-population rhs (same kinetics)."""
    acq = KERNELS["acquired_resistance"].rhs
    two = KERNELS["two_population_resistance"].rhs
    v = {"kg": 0.021, "kd": 0.3, "kgr": 0.025, "E": 1.0}
    y = [40.0, 5.0]
    a = acq(0.0, y, {**v, "alpha": 0.0})
    b = two(0.0, y, v)
    assert np.allclose(a, b)


def test_switching_flux_conserves_total_rate():
    """The acquired influx leaves S and enters R one-for-one, so the total dV/dt is independent
    of α at a fixed state (the α terms cancel in dS+dR)."""
    acq = KERNELS["acquired_resistance"].rhs
    v = {"kg": 0.021, "kd": 0.3, "kgr": 0.025, "E": 1.0}
    y = [40.0, 5.0]
    no_switch = sum(acq(0.0, y, {**v, "alpha": 0.0}))
    switch = sum(acq(0.0, y, {**v, "alpha": 0.05}))
    assert np.isclose(no_switch, switch)


def test_no_drug_no_acquired_resistance():
    """With E=0 the switching flux α·E·S vanishes — no resistance is generated."""
    acq = KERNELS["acquired_resistance"].rhs
    v = {"kg": 0.021, "kd": 0.3, "kgr": 0.025, "alpha": 0.05, "E": 0.0}
    ds_, dr = acq(0.0, [40.0, 0.0], v)
    assert dr == 0.0  # R stays at 0 when there is no drug and none pre-exists
    assert np.isclose(ds_, 0.021 * 40.0)  # pure sensitive-clone growth


# --- simulate-level: the origin divergence ----------------------------------


def test_acquired_nadir_is_shallower_than_pre_existing():
    """Matched on kg/kd/kgr, drug-driven switching limits the depth of response — the acquired
    model's nadir is less deep than the pre-existing model's."""
    ds = onkos.load()
    acq = onkos.simulate(ds, ACQUIRED, context=NSCLC, drug_effect=1.0, t=T).tumor_size
    pre = onkos.simulate(ds, PREEXISTING, context=NSCLC, drug_effect=1.0, t=T).tumor_size
    assert acq.min() > pre.min()


def test_acquired_progresses_earlier():
    """The acquired model reaches RECIST progression (v0.30 mechanistic TTP) earlier than the
    pre-existing model — the origin is visible in the tail."""
    ds = onkos.load()
    acq = onkos.simulate(ds, ACQUIRED, context=NSCLC, drug_effect=1.0, t=T).tumor_size
    pre = onkos.simulate(ds, PREEXISTING, context=NSCLC, drug_effect=1.0, t=T).tumor_size
    assert time_to_progression(T, acq) < time_to_progression(T, pre)


def test_week8_os_surrogate_is_blind_to_origin():
    """At week 8 and on the week-8-driven OS link the two origins agree — the surrogate cannot
    see the origin (the silent-transport risk this model surfaces)."""
    ds = onkos.load()
    acq = onkos.simulate(ds, ACQUIRED, context=NSCLC, drug_effect=1.0, t=T)
    pre = onkos.simulate(ds, PREEXISTING, context=NSCLC, drug_effect=1.0, t=T)
    assert abs(acq.median_os - pre.median_os) < 8.0  # OS surrogate barely distinguishes them
    assert abs(acq.tumor_size[16] - pre.tumor_size[16]) < 10.0  # ~week 8, both deep responders


# --- guardrails -------------------------------------------------------------


def test_tier_and_out_of_context_floor():
    ds = onkos.load()
    assert onkos.simulate(ds, ACQUIRED, context=NSCLC, drug_effect=1.0).tier == "C"
    floored = onkos.simulate(
        ds, ACQUIRED, context={"tumor_type": "melanoma", "line": "first"}, drug_effect=1.0
    )
    assert floored.tier == "D"  # transported outside its validated context


def test_record_in_nsclc_compare_set():
    ds = onkos.load()
    ids = [r.record_id for r in onkos.compare(ds, purpose="tgi", context=NSCLC).included]
    assert ACQUIRED in ids
