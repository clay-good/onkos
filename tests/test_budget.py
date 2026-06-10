"""Model-selection budget landmark suite — the decomposition *is* a balanced two-way
variance-component (ANOVA / first-order Sobol over the structural factors), not an
unconstrained attribution.

Mirrors ``test_combine.py``: closed-form properties of the variance-component algebra on
constructed grids, plus integration checks through ``model_selection_budget``.
"""

import numpy as np
import onkos
from onkos.budget import model_selection_budget, variance_components

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


# --- pure variance-component algebra ---------------------------------------


def test_components_sum_to_total_and_fractions_to_one():
    rng = np.random.default_rng(0)
    means = rng.normal(50, 10, size=(4, 3))
    within = rng.uniform(1, 20, size=(4, 3))
    c = variance_components(means, within)
    assert np.isclose(c["within"] + c["v_model"] + c["v_link"] + c["v_inter"], c["total"])
    assert np.isclose(sum(c["fractions"].values()), 1.0)


def test_all_components_non_negative():
    rng = np.random.default_rng(1)
    for _ in range(50):
        means = rng.normal(0, 5, size=(rng.integers(1, 5), rng.integers(1, 5)))
        within = rng.uniform(0, 5, size=means.shape)
        c = variance_components(means, within)
        for k in ("within", "v_model", "v_link", "v_inter", "between", "total"):
            assert c[k] >= -1e-12, k


def test_single_link_collapses_to_v021_within_between():
    """L=1 ⇒ V_link = V_inter = 0 and V_model = BETWEEN — the v0.21 within/between split."""
    means = np.array([[40.0], [70.0], [95.0]])      # 3 models, 1 link
    within = np.array([[2.0], [3.0], [4.0]])
    c = variance_components(means, within)
    assert np.isclose(c["v_link"], 0.0)
    assert np.isclose(c["v_inter"], 0.0)
    assert np.isclose(c["v_model"], c["between"])
    # BETWEEN equals the population variance of the per-model means (the v0.21 quantity).
    assert np.isclose(c["between"], np.var([40.0, 70.0, 95.0]))


def test_single_model_collapses_symmetrically():
    means = np.array([[40.0, 70.0, 95.0]])          # 1 model, 3 links
    c = variance_components(means, np.ones((1, 3)))
    assert np.isclose(c["v_model"], 0.0)
    assert np.isclose(c["v_inter"], 0.0)
    assert np.isclose(c["v_link"], c["between"])


def test_identical_cells_zero_between():
    means = np.full((3, 4), 55.0)
    c = variance_components(means, np.full((3, 4), 7.0))
    assert np.isclose(c["between"], 0.0)
    assert np.isclose(c["fractions"]["parameter"], 1.0)


def test_pure_link_main_effect():
    """Cell means depend only on the link ⇒ all between variance is V_link."""
    col = np.array([30.0, 60.0, 90.0])
    means = np.vstack([col, col, col])              # identical rows
    c = variance_components(means, np.ones((3, 3)))
    assert np.isclose(c["v_model"], 0.0)
    assert np.isclose(c["v_inter"], 0.0)
    assert np.isclose(c["v_link"], c["between"])


def test_additive_layout_has_zero_interaction():
    """μ_ab = r_a + c_b ⇒ no interaction term."""
    r = np.array([10.0, 20.0, 35.0])[:, None]
    c_ = np.array([1.0, 4.0, 9.0, 16.0])[None, :]
    means = r + c_
    comp = variance_components(means, np.zeros_like(means))
    assert np.isclose(comp["v_inter"], 0.0, atol=1e-9)


def test_matches_direct_two_way_anova_sum_of_squares():
    """The component variances equal a balanced two-way SS decomposition ÷ N."""
    rng = np.random.default_rng(2)
    means = rng.normal(0, 3, size=(5, 4))
    M, L = means.shape
    grand = means.mean()
    ss_a = L * np.sum((means.mean(axis=1) - grand) ** 2)
    ss_b = M * np.sum((means.mean(axis=0) - grand) ** 2)
    ss_total = np.sum((means - grand) ** 2)
    ss_ab = ss_total - ss_a - ss_b
    c = variance_components(means, np.zeros_like(means))
    N = M * L
    assert np.isclose(c["v_model"], ss_a / N)
    assert np.isclose(c["v_link"], ss_b / N)
    assert np.isclose(c["v_inter"], ss_ab / N)


def test_zero_within_is_purely_structural():
    means = np.array([[40.0, 60.0], [50.0, 90.0]])
    c = variance_components(means, np.zeros((2, 2)))
    assert np.isclose(c["within"], 0.0)
    assert np.isclose(c["fractions"]["parameter"], 0.0)


def test_grand_mean_in_convex_hull():
    rng = np.random.default_rng(3)
    means = rng.normal(50, 12, size=(4, 3))
    c = variance_components(means, np.ones((4, 3)))
    assert means.min() - 1e-9 <= c["grand_mean"] <= means.max() + 1e-9


# --- integration through model_selection_budget ----------------------------


def test_nsclc_budget_has_all_four_components_and_a_dominant_axis():
    ds = onkos.load()
    b = model_selection_budget(ds, context=NSCLC, endpoint="OS", n=80)
    assert len(b.models) >= 3 and len(b.links) >= 2     # the rich 4x3-ish grid
    assert np.isclose(sum(b.fractions.values()), 1.0)
    assert all(v >= -1e-12 for v in b.fractions.values())
    assert b.dominant in b.fractions
    assert 0.0 <= b.structural_fraction <= 1.0
    assert b.within > 0                                  # parameter IIV exists


def test_budget_tier_is_worst_included_and_carries_clinical_use():
    ds = onkos.load()
    b = model_selection_budget(ds, context=NSCLC, endpoint="OS", n=60)
    worst = max(tr.tier for tr in onkos.compare(ds, purpose="tgi", context=NSCLC).included)
    assert b.tier == worst
    d = b.to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True and "PROHIBITED" in d["onkos:clinicalUse"]
    assert "structural_fraction" in d and "dominant" in d


def test_single_link_context_flags_and_zeroes_v_link():
    """A context with one OS survival link reduces to the v0.21 split and is flagged."""
    ds = onkos.load()
    b = model_selection_budget(ds, context={"tumor_type": "breast", "line": "first"},
                               endpoint="OS", n=60)
    assert len(b.links) == 1
    assert np.isclose(b.v_link, 0.0) and np.isclose(b.v_inter, 0.0)
    assert any("single_survival_link" in w for w in b.warnings)
