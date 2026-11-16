"""Readout noise for the early-surrogate timing analysis (earliness vs precision).

Follow-on to ``docs/specs/research/early-surrogate-timing.md`` §7. Closed-form landmarks for
the noise model, then the binding's limits: vanishing noise reproduces the noise-free
counts, overwhelming noise turns every pair into a coin flip.
"""

import numpy as np
import onkos
import pytest
from onkos.early_surrogate import (
    ranking_flip_probability,
    readout_se,
    surrogate_timing_fidelity,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


def test_readout_se_closed_form():
    s = np.sqrt(np.log1p(0.3**2))
    assert np.isclose(readout_se(0.3, 50), 2 * s / np.sqrt(50))
    assert readout_se(0.0, 10) == 0.0
    assert np.isclose(readout_se(0.3, 200), readout_se(0.3, 50) / 2)  # 4x patients, half SE


def test_flip_probability_limits_and_symmetry():
    assert ranking_flip_probability(0.0, 0.1) == 0.5
    assert ranking_flip_probability(0.2, 0.0) == 0.0
    assert ranking_flip_probability(0.0, 0.0) == 0.5
    assert np.isclose(ranking_flip_probability(0.1, 0.1), 0.15865525, atol=1e-7)  # Phi(-1)
    assert ranking_flip_probability(-0.3, 0.2) == ranking_flip_probability(0.3, 0.2)
    assert ranking_flip_probability(0.1, 0.2) > ranking_flip_probability(0.3, 0.2)


def test_noise_inputs_are_validated():
    for cv in (-0.1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            readout_se(cv, 10)
    for n in (0, 2.5):
        with pytest.raises(ValueError):
            readout_se(0.3, n)


def test_vanishing_noise_reproduces_the_noise_free_counts():
    ds = onkos.load()
    st = surrogate_timing_fidelity(ds, context=NSCLC, noise_cv=1e-9, n_patients=50)
    for r in st.rows:
        assert np.isclose(r["expected_discordant_pairs"], r["discordant_pairs"])
        assert r["expected_flips"] < 1e-6


def test_flips_grow_with_noise_and_never_pass_a_coin_flip():
    """Each pair's flip probability rises toward 1/2 as the noise grows, so the expected
    flips per landmark increase monotonically and stay at or below half the pairs."""
    ds = onkos.load()
    prev = None
    for cv in (0.1, 0.5, 2.0, 50.0):
        st = surrogate_timing_fidelity(ds, context=NSCLC, noise_cv=cv, n_patients=10)
        flips = [r["expected_flips"] for r in st.rows]
        assert all(f <= st.total_pairs / 2 + 1e-9 for f in flips)
        if prev is not None:
            assert all(f >= p - 1e-9 for f, p in zip(flips, prev))
        prev = flips


def test_early_readouts_are_the_least_reproducible():
    """ctDNA-like noise (CV 0.3, 50 patients per arm): the week-2 readout flips far more
    model pairs than the week-8 one, and the week-52 readout is essentially exact."""
    ds = onkos.load()
    st = surrogate_timing_fidelity(ds, context=NSCLC, noise_cv=0.3, n_patients=50)
    flips = {r["week"]: r["expected_flips"] for r in st.rows}
    assert flips[2.0] > 3 * flips[8.0]
    assert flips[52.0] < 0.05
    d = st.to_dict()
    assert d["noise_cv"] == 0.3 and "expected_flips" in d["rows"][0]


def test_default_output_is_unchanged():
    ds = onkos.load()
    d = surrogate_timing_fidelity(ds, context=NSCLC).to_dict()
    assert "noise_cv" not in d and "expected_flips" not in d["rows"][0]
