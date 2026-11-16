"""The exposure-response shape as a third factor of the model-selection budget.

Follow-on to ``docs/specs/research/exposure-response-extrapolation.md`` §7. The three-factor
split is the functional ANOVA of :func:`anova_components`; landmarks first, then the
binding, whose key property is that the ER share vanishes exactly at the studied dose.
"""

import numpy as np
import onkos
import pytest
from onkos.budget import anova_components, model_selection_budget, variance_components

NSCLC_2L = {"tumor_type": "NSCLC", "line": "second"}  # small grid: 2 models x 1 link
N = 20


# --- functional ANOVA ---------------------------------------------------------


def test_two_factor_anova_is_the_existing_decomposition():
    means = np.random.default_rng(0).normal(60, 15, size=(4, 3))
    c = variance_components(means, np.zeros_like(means))
    a = anova_components(means)
    assert np.isclose(a[(0,)], c["v_model"])
    assert np.isclose(a[(1,)], c["v_link"])
    assert np.isclose(a[(0, 1)], c["v_inter"])


def test_three_factor_effects_are_orthogonal_and_non_negative():
    means = np.random.default_rng(1).normal(60, 15, size=(4, 3, 2))
    a = anova_components(means)
    assert len(a) == 7
    assert all(v >= 0 for v in a.values())
    assert np.isclose(sum(a.values()), means.var())


def test_an_additive_grid_has_no_interactions():
    x, y, z = np.array([1.0, 4.0, 9.0]), np.array([0.0, 2.0]), np.array([5.0, -1.0, 3.0, 7.0])
    means = x[:, None, None] + y[None, :, None] + z[None, None, :]
    a = anova_components(means)
    assert np.isclose(a[(0,)], x.var()) and np.isclose(a[(1,)], y.var())
    assert np.isclose(a[(2,)], z.var())
    assert all(np.isclose(v, 0.0) for k, v in a.items() if len(k) > 1)


def test_a_constant_third_factor_contributes_nothing():
    base = np.random.default_rng(2).normal(size=(3, 2))
    a = anova_components(np.repeat(base[:, :, None], 4, axis=2))
    assert all(np.isclose(v, 0.0) for k, v in a.items() if 2 in k)
    assert np.isclose(a[(0, 1)], variance_components(base, np.zeros_like(base))["v_inter"])


# --- binding --------------------------------------------------------------------


def test_er_share_is_exactly_zero_at_the_studied_dose():
    """Every ER shape is anchored to the same effect at c_ref, so at er_dose == c_ref the
    third factor is constant and the budget is the two-factor one."""
    ds = onkos.load()
    two = model_selection_budget(ds, context=NSCLC_2L, n=N)
    three = model_selection_budget(ds, context=NSCLC_2L, n=N, er_dose=150.0, c_ref=150.0)
    assert len(three.er_models) == 3
    assert np.isclose(three.v_er, 0.0) and np.isclose(three.v_er_total, 0.0)
    assert np.allclose(three.cell_means[:, :, 0], two.cell_means)
    for k in ("parameter", "tgi_model", "survival_link", "interaction"):
        assert np.isclose(three.fractions[k], two.fractions[k])
    assert "er_models" not in two.to_dict()


def test_er_share_appears_on_extrapolation():
    ds = onkos.load()
    b = model_selection_budget(ds, context=NSCLC_2L, n=N, er_dose=37.5)
    assert b.v_er > 0 and b.er_total_fraction >= b.fractions["er_model"] > 0
    assert np.isclose(sum(b.fractions.values()), 1.0)
    assert np.isclose(b.structural_fraction + b.fractions["parameter"], 1.0)
    d = b.to_dict()
    assert d["er_dose"] == 37.5 and "v_er" in d["components"] and d["er_total_fraction"] > 0


def test_er_factor_rejects_bad_inputs():
    ds = onkos.load()
    for bad in (0.0, -5.0, float("nan")):
        with pytest.raises(ValueError, match="positive"):
            model_selection_budget(ds, context=NSCLC_2L, n=2, er_dose=bad)
    with pytest.raises(ValueError, match="drug_effect"):
        model_selection_budget(ds, context=NSCLC_2L, n=2, er_dose=50.0, drug_effect=0.0)
    with pytest.raises(ValueError, match="at least one"):
        model_selection_budget(ds, context=NSCLC_2L, n=2, er_dose=50.0, er_ids=[])
