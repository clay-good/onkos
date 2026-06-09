"""Variance-based parameter sensitivity (which IIV drives the prediction)."""

import numpy as np
import onkos
import pytest
from onkos.sensitivity import sensitivity

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


def test_contributions_sum_to_one_and_are_ranked():
    res = sensitivity(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC, n=300, seed=0)
    contribs = [p.contribution for p in res.indices]
    assert abs(sum(contribs) - 1.0) < 1e-6
    assert contribs == sorted(contribs, reverse=True)
    assert 0.0 <= res.r_squared <= 1.0


def test_kill_rate_dominates_median_os_for_claret():
    res = sensitivity(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC,
                      target="median_os_weeks", n=400, seed=0)
    assert res.dominant.symbol == "kD"
    assert res.dominant.contribution > 0.5
    # higher kill -> longer OS: positive SRC on median OS
    assert res.dominant.src > 0


def test_direction_flips_for_week8_change():
    # higher kill -> deeper (more negative) week-8 change
    res = sensitivity(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC,
                      target="week8_relative_change", n=400, seed=0)
    kd = next(p for p in res.indices if p.symbol == "kD")
    assert kd.src < 0


def test_reproducible_with_seed():
    a = sensitivity(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC, n=200, seed=5)
    b = sensitivity(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC, n=200, seed=5)
    assert [p.symbol for p in a.indices] == [p.symbol for p in b.indices]
    assert np.allclose([p.src for p in a.indices], [p.src for p in b.indices])


def test_only_iiv_parameters_are_analyzed():
    # the Claret model has 3 IIV params (kL, kD, lambda); growth-law Vmax has none
    res = sensitivity(onkos.load(), "resistance.claret_2009.tgi", context=NSCLC, n=120, seed=0)
    assert {p.symbol for p in res.indices} == {"kL", "kD", "lambda"}


def test_record_without_iiv_raises():
    # growth laws carry no IIV CV -> nothing to analyze
    with pytest.raises(ValueError, match="no IIV parameters"):
        sensitivity(onkos.load(), "growth_laws.gompertz", context=NSCLC, n=10)
