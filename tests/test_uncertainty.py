"""Monte-Carlo parameter-uncertainty propagation (IIV CV -> prediction bands)."""

import numpy as np
import onkos
from onkos.uncertainty import simulate_ensemble

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


def test_bands_are_ordered():
    ds = onkos.load()
    ens = simulate_ensemble(ds, "resistance.claret_2009.tgi", context=NSCLC, n=120, seed=0)
    b = ens.tumor_size
    assert np.all(b.lo <= b.median + 1e-9)
    assert np.all(b.median <= b.hi + 1e-9)
    assert ens.os_curve is not None
    o = ens.os_curve
    assert np.all(o.lo <= o.hi + 1e-9)


def test_reproducible_with_seed():
    ds = onkos.load()
    a = simulate_ensemble(ds, "resistance.claret_2009.tgi", context=NSCLC, n=80, seed=7)
    b = simulate_ensemble(ds, "resistance.claret_2009.tgi", context=NSCLC, n=80, seed=7)
    assert np.allclose(a.tumor_size.median, b.tumor_size.median)
    assert np.allclose(a.tumor_size.hi, b.tumor_size.hi)


def test_no_iiv_gives_degenerate_band():
    ds = onkos.load()
    # growth laws carry no IIV CV -> every sample identical -> zero-width band
    ens = simulate_ensemble(ds, "growth_laws.gompertz", context=NSCLC, drug_effect=0.0, n=40)
    assert np.allclose(ens.tumor_size.lo, ens.tumor_size.hi)


def test_higher_cv_widens_the_band():
    """A record whose kill/resistance CVs are larger yields a wider week-8 band."""
    ds = onkos.load()
    import copy

    rid = "resistance.claret_2009.tgi"
    wide = simulate_ensemble(ds, rid, context=NSCLC, n=200, seed=1)

    # Build a low-IIV clone of the same record (Parameter is a frozen dataclass).
    low = copy.deepcopy(ds[rid])
    for p in low.parameters:
        if p.iiv_cv_percent:
            object.__setattr__(p, "iiv_cv_percent", 5.0)
    ds._records[rid] = low  # type: ignore[attr-defined]
    narrow = simulate_ensemble(ds, rid, context=NSCLC, n=200, seed=1)

    def width(ens):
        m = ens.metrics["week8_relative_change"]
        return m["hi"] - m["lo"]

    assert width(wide) > width(narrow)


def test_tier_and_warnings_preserved():
    ds = onkos.load()
    # out-of-context still floors to D, independent of sampling
    ens = simulate_ensemble(
        ds, "resistance.claret_2009.tgi", context={"tumor_type": "colorectal", "line": "first"}, n=20
    )
    assert ens.tier == "D"
    assert any("outside validated" in w for w in ens.warnings)


def test_median_os_interval_brackets_point_estimate():
    ds = onkos.load()
    point = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0)
    ens = simulate_ensemble(ds, "resistance.claret_2009.tgi", context=NSCLC, n=300, seed=0)
    mos = ens.metrics["median_os_weeks"]
    assert mos["lo"] < mos["hi"]
    # the deterministic point estimate lies within a wide CI
    assert mos["lo"] - 5 <= point.median_os <= mos["hi"] + 5
