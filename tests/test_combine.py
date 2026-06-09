"""Model-averaging landmark suite — the combination estimator *is* the law of
total variance and a convex forecast combination, not an unconstrained fit.

Mirrors ``test_landmarks.py``: closed-form properties of the combination math
itself (research spec §6), proven on constructed inputs with known moments, plus
a few integration checks through ``Comparison.model_average``.
"""

import numpy as np
import onkos
from onkos.combine import (
    SCHEMES,
    average_curve,
    compute_weights,
    decompose,
    tier_scores,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


def _cmp(ctx=NSCLC, t=None):
    if t is None:
        t = np.linspace(0.0, 156.0, 313)
    return onkos.compare(onkos.load(), purpose="tgi", context=ctx, drug_effect=1.0, t=t)


# --- pure combination math -------------------------------------------------


def test_equal_weight_identity():
    """Q̄ under equal weights is the arithmetic mean of the per-model means."""
    means = [40.0, 70.0, 95.0]
    w, _ = compute_weights("equal", ["A", "B", "C"], [None, None, None])
    point, *_ = decompose(means, [1.0, 2.0, 3.0], w)
    assert np.isclose(point, np.mean(means))


def test_weights_sum_to_one_every_scheme():
    tiers = ["A", "B", "C"]
    evidence = [0.66, 0.60, None]
    for scheme in SCHEMES:
        w, _ = compute_weights(scheme, tiers, evidence)
        assert np.isclose(w.sum(), 1.0)
        assert np.all(w >= 0)


def test_zero_weight_model_is_inert():
    """A zero-weight model leaves Q̄, S̄, and both variances unchanged."""
    means, varis = [40.0, 70.0, 1e6], [1.0, 2.0, 9.0]
    w2 = np.array([0.5, 0.5])
    w3 = np.array([0.5, 0.5, 0.0])
    a = decompose(means[:2], varis[:2], w2)
    b = decompose(means, varis, w3)
    assert np.allclose(a, b)

    curves = np.array([np.linspace(1, 0.2, 10), np.linspace(1, 0.5, 10), np.zeros(10)])
    s2, btw2, _ = average_curve(curves[:2], w2)
    s3, btw3, _ = average_curve(curves, w3)
    assert np.allclose(s2, s3) and np.allclose(btw2, btw3)


def test_identical_components_have_zero_between_any_weights():
    means = [55.0, 55.0, 55.0]
    for w in (np.array([1.0, 1, 1]), np.array([0.7, 0.2, 0.1]), np.array([0.0, 0, 1])):
        _, _, between, frac = decompose(means, [3.0, 7.0, 2.0], w)
        assert np.isclose(between, 0.0)
        assert np.isclose(frac, 0.0)


def test_law_of_total_variance_on_lognormal_mixture():
    """Total = WITHIN + BETWEEN to ≤1e-9 on a mixture of lognormals with known
    analytic moments (the decomposition is the law, not a coincidence)."""
    rng = np.random.default_rng(0)
    mus = np.array([0.3, 0.9, -0.2])
    sigmas = np.array([0.4, 0.7, 0.5])
    w = _normalize(rng.random(3))
    means = np.exp(mus + sigmas**2 / 2)
    varis = (np.exp(sigmas**2) - 1) * np.exp(2 * mus + sigmas**2)
    point, within, between, _ = decompose(means, varis, w)
    # True variance of the mixture distribution: E[Q²] − E[Q]².
    true_total = float(np.sum(w * (varis + means**2)) - point**2)
    assert abs(true_total - (within + between)) < 1e-9


def test_convex_hull_bound():
    """Q̄ ∈ [min, max] of the per-model means; never extrapolates."""
    means = [40.0, 70.0, 95.0]
    for w in (np.array([1.0, 1, 1]), np.array([0.8, 0.1, 0.1]), np.array([0.1, 0.1, 0.8])):
        point, *_ = decompose(means, [1.0, 1, 1], w)
        assert min(means) - 1e-9 <= point <= max(means) + 1e-9


def test_curve_convex_hull_bound():
    curves = np.array([np.linspace(1, 0.2, 20), np.linspace(1, 0.6, 20), np.linspace(1, 0.4, 20)])
    s_bar, _, _ = average_curve(curves, np.array([0.2, 0.5, 0.3]))
    assert np.all(s_bar <= curves.max(axis=0) + 1e-12)
    assert np.all(s_bar >= curves.min(axis=0) - 1e-12)


def test_monotone_reweighting_moves_point_toward_that_model():
    """Raising one model's weight moves Q̄ monotonically toward its mean."""
    means = np.array([40.0, 95.0])
    target = means[1]
    prev = None
    for a in np.linspace(0.0, 1.0, 11):
        point, *_ = decompose(means, [1.0, 1.0], np.array([1 - a, a]))
        if prev is not None:
            assert point >= prev - 1e-12  # moving toward the larger mean
        prev = point
    assert np.isclose(prev, target)


def test_tier_scores_ratio_is_four_two_one():
    s = tier_scores()
    assert s["A"] / s["C"] == 4.0 and s["B"] / s["C"] == 2.0
    assert s["A"] > s["B"] > s["C"] > s["D"]


def test_evidence_scheme_falls_back_when_no_validation():
    w, notes = compute_weights("evidence", ["C", "C"], [None, None])
    assert np.allclose(w, [0.5, 0.5])
    assert any("unavailable" in n for n in notes)


def _normalize(w):
    w = np.asarray(w, float)
    return w / w.sum()


# --- integration through Comparison.model_average --------------------------


def test_survival_function_validity():
    """S̄(t) is a valid survival function: S̄(0)=1, monotone non-increasing, in [0,1]."""
    ma = _cmp().model_average(target="median_os_weeks", endpoint="OS", weights="equal", n=60)
    s = ma.curve
    assert np.isclose(s[0], 1.0, atol=1e-6)
    assert np.all(np.diff(s) <= 1e-9)
    assert np.all((s >= -1e-12) & (s <= 1.0 + 1e-12))


def test_averaged_tier_is_worst_included():
    cmp = _cmp()
    ma = cmp.model_average(n=40)
    worst = max(tr.tier for tr in cmp.included)  # 'A'<'B'<'C'<'D' lexicographically
    assert ma.tier == worst


def test_single_eligible_model_is_flagged_and_fraction_zero():
    cmp = _cmp()
    cmp.included = cmp.included[:1]  # degenerate set, M=1
    ma = cmp.model_average(n=40)
    assert ma.model_selection_fraction == 0.0
    assert np.isclose(ma.between_var, 0.0)
    assert any("single_eligible_model" in w for w in ma.warnings)
    # S̄ ≡ S_1 (the lone model's own within-mean curve)
    assert np.all(np.isfinite(ma.curve))


def test_identical_components_zero_between_integration():
    cmp = _cmp()
    one = cmp.included[0]
    cmp.included = [one, one, one]  # three copies of one model
    ma = cmp.model_average(n=40)
    assert np.isclose(ma.between_var, 0.0, atol=1e-6)
    assert np.isclose(ma.model_selection_fraction, 0.0, atol=1e-6)


def test_result_carries_clinical_use_and_fraction_inseparably():
    ma = _cmp().model_average(n=40)
    d = ma.to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "PROHIBITED" in d["onkos:clinicalUse"]
    assert "model_selection_fraction" in d  # never reported without its disagreement
    assert d["onkos:modelSelectionUncertainty"]["fraction"] == ma.model_selection_fraction
    assert "NOT posterior" in d["weights_meaning"]


def test_decomposition_table_spans_schemes():
    dec = _cmp().uncertainty_decomposition(target="median_os_weeks", n=40)
    assert set(dec) == set(SCHEMES)
    for row in dec.values():
        assert {"point", "within_var", "between_var", "model_selection_fraction"} <= set(row)
        assert 0.0 <= row["model_selection_fraction"] <= 1.0
