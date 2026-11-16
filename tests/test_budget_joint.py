"""The joint (current-value) link as a level of the budget's survival-link factor.

Follow-on to ``docs/specs/research/joint-survival.md`` §7: the budget decomposed only the
proportional-hazards links; ``joint_alpha`` adds the non-proportional current-value link
as one more level. Closed-form checks first, then integration through
``model_selection_budget``.
"""

import warnings

import numpy as np
import onkos
import pytest
from onkos.budget import _joint_medians, joint_link_label, model_selection_budget
from onkos.joint import current_value_hazard_ratio, current_value_survival

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
N = 30


# --- closed form -------------------------------------------------------------


def test_constant_trajectory_joint_median_is_the_weibull_ph_median():
    """A trajectory held at c x baseline has HR = c^alpha forever, so the joint median is
    the proportional-hazards Weibull median scale * (ln 2 / HR)^(1/shape)."""
    t = np.linspace(0.0, 400.0, 4001)
    shape, scale, alpha, c = 1.3, 80.0, 1.0, 0.5
    tumor = np.full((3, t.size), 40.0 * c)
    tumor[:, 0] = 40.0  # baseline; the size steps to c x baseline after the first step
    med = _joint_medians(tumor, t, alpha, shape, scale)
    # HR jumps from 1 to c over the first step (0.1 wk); away from it the curve is PH.
    expected = scale * (np.log(2.0) / c ** alpha) ** (1.0 / shape)
    assert np.allclose(med, expected, rtol=2e-3)


def test_hazard_ratio_cap_keeps_extreme_alpha_finite():
    v = np.array([10.0, 5.0, 1e-6, 40.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no overflow RuntimeWarning
        hr = current_value_hazard_ratio(v, 1e9)
        s = current_value_survival(np.array([0.0, 1.0, 2.0, 3.0]), hr, shape=1.3, scale=80.0)
    assert np.all(np.isfinite(hr)) and hr[0] == 1.0
    assert np.all(np.isfinite(s)) and np.all((s >= 0) & (s <= 1))


def test_label_is_stable():
    assert joint_link_label(1.0) == "joint:current_value(alpha=1)"
    assert joint_link_label(0.5) == "joint:current_value(alpha=0.5)"


# --- integration -------------------------------------------------------------


def test_default_budget_is_unchanged_and_alpha_zero_adds_a_flat_baseline_column():
    """alpha=0 makes HR = 1: every model's joint cell is the baseline Weibull median with
    no parameter variance, and the proportional-hazards cells are untouched."""
    ds = onkos.load()
    base = model_selection_budget(ds, context=NSCLC, n=N)
    b0 = model_selection_budget(ds, context=NSCLC, n=N, joint_alpha=0.0)
    assert b0.models == base.models
    assert b0.links == base.links + [joint_link_label(0.0)]
    assert np.allclose(b0.cell_means[:, :-1], base.cell_means)
    col = b0.cell_means[:, -1]
    weibull = ds["survival_link.nsclc_os_week8"]
    shape = weibull["weibull_shape"].central
    scale = weibull["weibull_scale"].central
    assert np.allclose(col, scale * np.log(2.0) ** (1.0 / shape), rtol=5e-3)
    assert np.allclose(col, col[0])
    assert not any(joint_link_label(0.0) in w for w in base.warnings)


def test_joint_level_drops_the_curative_model_with_a_censoring_warning():
    """Under a current-value hazard the eradicating Norton-Simon model's median lies
    beyond the horizon for most samples. The cell is treated as missing (not averaged
    over the non-responders), the model is dropped, and the joint level is kept."""
    ds = onkos.load()
    b = model_selection_budget(ds, context=NSCLC, n=N, joint_alpha=1.0)
    assert b.links[-1] == joint_link_label(1.0)
    assert "drug_effect.norton_simon.nsclc" not in b.models
    assert len(b.models) >= 3
    assert any(w.startswith("censored_median") for w in b.warnings)
    assert np.all(np.isfinite(b.cell_means))
    assert np.isclose(sum(b.fractions.values()), 1.0)
    assert b.v_link > 0


def test_joint_alpha_rejects_pfs_and_non_finite():
    ds = onkos.load()
    with pytest.raises(ValueError, match="OS link"):
        model_selection_budget(ds, context=NSCLC, endpoint="PFS", n=2, joint_alpha=1.0)
    for bad in (float("nan"), float("inf")):
        with pytest.raises(ValueError, match="finite"):
            model_selection_budget(ds, context=NSCLC, n=2, joint_alpha=bad)
