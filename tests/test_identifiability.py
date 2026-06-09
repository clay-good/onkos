"""Practical-identifiability landmark suite — the analyzer *is* the Fisher-
information Cramér-Rao bound and the Brun collinearity index, not an unconstrained
precision guess.

Mirrors ``test_combine.py``: closed-form properties of the information algebra
itself (research spec §5), proven on constructed inputs, plus integration checks
through ``onkos.identifiability`` on real records.
"""

import numpy as np
import onkos
from onkos.identify import (
    collinearity_index,
    crlb_rse,
    fisher_information,
    identifiability,
)

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


# --- pure information algebra ----------------------------------------------


def test_information_is_additive_over_observations():
    """M(A ∪ B) = M(A) + M(B) — Fisher information adds over independent rows."""
    rng = np.random.default_rng(0)
    a = rng.normal(size=(4, 3))
    b = rng.normal(size=(5, 3))
    both = np.vstack([a, b])
    assert np.allclose(fisher_information(both), fisher_information(a) + fisher_information(b))


def test_more_observations_never_raise_rse():
    """Adding an observation can only add (PSD) information ⇒ no RSE rises."""
    rng = np.random.default_rng(1)
    s = rng.normal(size=(6, 3))
    theta = np.array([1.0, 2.0, 0.5])
    _, _, rse_sub = crlb_rse(s[:4], theta)
    _, _, rse_full = crlb_rse(s, theta)
    assert np.all(rse_full <= rse_sub + 1e-9)


def test_residual_error_scaling_scales_rse_linearly():
    """Scaling every σ by c (i.e. S̃ by 1/c) scales every RSE by c."""
    rng = np.random.default_rng(2)
    s = rng.normal(size=(8, 3))
    theta = np.array([1.0, 2.0, 3.0])
    _, _, rse1 = crlb_rse(s, theta)
    _, _, rse2 = crlb_rse(s / 2.0, theta)  # σ doubled
    assert np.allclose(rse2, 2.0 * rse1)


def test_structural_nonidentifiability_is_infinite_not_regularized():
    """Two identical (parallel) columns ⇒ singular FIM ⇒ inf RSE and inf γ_K."""
    col = np.linspace(1, 2, 5)
    s = np.column_stack([col, col, np.ones(5)])  # cols 0,1 identical
    _, _, rse = crlb_rse(s, np.array([1.0, 1.0, 1.0]))
    assert np.all(np.isinf(rse))
    assert np.isinf(collinearity_index(s))


def test_orthogonal_design_has_collinearity_one():
    """Orthogonal weighted sensitivity columns ⇒ γ_K = 1 (its lower bound)."""
    s = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
    assert np.isclose(collinearity_index(s), 1.0)


def test_collinearity_index_is_scale_invariant():
    """Rescaling one parameter's column leaves γ_K unchanged (it is normalized)."""
    rng = np.random.default_rng(3)
    s = rng.normal(size=(7, 3))
    s2 = s.copy()
    s2[:, 1] *= 1000.0
    assert np.isclose(collinearity_index(s), collinearity_index(s2))


def test_collinearity_index_floor_is_one():
    rng = np.random.default_rng(4)
    for _ in range(20):
        s = rng.normal(size=(6, 3))
        assert collinearity_index(s) >= 1.0 - 1e-9


def test_crlb_consistency_with_direct_inverse():
    """RSE = sqrt(diag((S̃ᵀS̃)⁻¹)) / |θ|, computed independently."""
    rng = np.random.default_rng(5)
    s = rng.normal(size=(10, 3))
    theta = np.array([2.0, 0.5, 4.0])
    _, _, rse = crlb_rse(s, theta)
    direct = np.sqrt(np.diag(np.linalg.inv(s.T @ s))) / np.abs(theta)
    assert np.allclose(rse, direct)


# --- integration through onkos.identifiability -----------------------------


def test_exponential_crlb_matches_closed_form():
    """For y = V0·e^{kt} with proportional error, RSE_k = σ_prop / (|k|·sqrt(Σ t²))."""
    ds = onkos.load()
    sched = [2.0, 4.0, 8.0, 16.0, 24.0]
    sp = 0.2
    res = identifiability(ds, "growth_laws.exponential", context=NSCLC, schedule=sched,
                          sigma_prop=sp)
    p = res.params[0]
    closed = sp / (abs(p.central) * np.sqrt(np.sum(np.asarray(sched) ** 2))) * 100.0
    assert np.isclose(p.rse_percent, closed, rtol=1e-3)
    assert np.isclose(res.collinearity_index, 1.0)  # a single parameter is orthogonal to itself


def test_resistance_term_is_a_flat_likelihood_artifact():
    """The Claret λ term: high stored CV AND practically unidentifiable on a short
    schedule ⇒ the cv_is_identifiability_artifact flag fires; the model as a whole is
    not practically identifiable under a realistic RECIST cadence."""
    ds = onkos.load()
    res = identifiability(ds, "resistance.claret_2009.tgi", context=NSCLC)
    assert not res.practically_identifiable
    lam = next(p for p in res.params if p.symbol == "lambda")
    assert lam.iiv_cv_percent and lam.iiv_cv_percent >= 50
    assert not lam.identifiable
    assert any("cv_is_identifiability_artifact" in w and "lambda" in w for w in res.warnings)


def test_lengthening_the_schedule_improves_precision():
    """A superset schedule cannot worsen any parameter's predicted RSE (end-to-end
    monotonic-information property)."""
    ds = onkos.load()
    short = identifiability(ds, "resistance.claret_2009.tgi", context=NSCLC,
                            schedule=[0, 6, 12, 18, 24])
    long = identifiability(ds, "resistance.claret_2009.tgi", context=NSCLC,
                           schedule=[0, 6, 12, 18, 24, 48, 72, 104])
    by_symbol = {p.symbol: p.rse_percent for p in short.params}
    for p in long.params:
        assert p.rse_percent <= by_symbol[p.symbol] + 1e-6


def test_tier_passes_through_unchanged():
    """Identifiability cannot move a tier — the result carries the record's own
    propagated tier (and an out-of-context transport still tiers to D)."""
    ds = onkos.load()
    in_ctx = identifiability(ds, "resistance.claret_2009.tgi", context=NSCLC)
    assert in_ctx.tier == ds["resistance.claret_2009.tgi"].tier
    # Transported outside its validated context: propagation floors the tier to D,
    # and the analyzer passes that through rather than masking it.
    out = identifiability(ds, "resistance.claret_2009.tgi",
                          context={"tumor_type": "melanoma", "line": "first"})
    assert out.tier == "D"


def test_non_ode_kernel_is_rejected():
    """Identifiability of a *trajectory* is undefined for survival / transform kernels."""
    ds = onkos.load()
    survival_id = next(r.id for r in ds if r.purpose == "survival_link")
    try:
        identifiability(ds, survival_id)
        raise AssertionError("expected ValueError for a non-ODE kernel")
    except ValueError as e:
        assert "ODE" in str(e) or "dynamic" in str(e)


def test_result_carries_clinical_use_and_named_design():
    ds = onkos.load()
    d = identifiability(ds, "resistance.claret_2009.tgi", context=NSCLC).to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "PROHIBITED" in d["onkos:clinicalUse"]
    assert "schedule_weeks" in d["design"] and "sigma_prop" in d["design"]
    assert "practically_identifiable" in d  # the verdict travels with the design
