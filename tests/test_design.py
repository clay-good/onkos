"""D-optimal trial-design landmark suite (research spec optimal-design §5): the best
sampling schedule a fixed budget allows, from the v0.22 design Fisher information. The
payload: the optimal design rescues the *circumstantially* flat parameter but the deeply
*structurally* flat growth rate stays unidentifiable under the best possible schedule.
"""

import numpy as np
import onkos
import pytest
from onkos.design import d_optimal_rows, optimal_schedule
from onkos.identify import fisher_information

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CLARET = "resistance.claret_2009.tgi"
BIEXP = "tgi_metrics.wang_2009.biexponential"


# --- d_optimal_rows: the pure selection core --------------------------------


def test_closed_form_selection():
    """On a matrix with two large orthogonal rows and two tiny ones, the D-optimal 2-subset
    is the two large orthogonal rows (they maximize det of the summed outer products)."""
    ss = np.array([[10.0, 0.0], [0.0, 10.0], [1.0, 0.0], [0.0, 1.0]])
    assert sorted(d_optimal_rows(ss, 2)) == [0, 1]


def test_additivity_of_information():
    """The Fisher information of a schedule is the sum of its rows' outer products — the
    property the optimizer relies on."""
    ss = np.array([[2.0, 1.0], [1.0, 3.0], [0.0, 2.0]])
    rows = [0, 2]
    direct = fisher_information(ss[rows])
    summed = sum(np.outer(ss[i], ss[i]) for i in rows)
    assert np.allclose(direct, summed)


def test_seed_rows_are_kept():
    ss = np.array([[10.0, 0.0], [0.0, 10.0], [1.0, 1.0], [5.0, 5.0]])
    chosen = d_optimal_rows(ss, 3, seed_rows=(2,))
    assert 2 in chosen and len(chosen) == 3


def test_budget_bounds():
    ss = np.eye(3)
    with pytest.raises(ValueError):
        d_optimal_rows(ss, 4)  # more than the candidate rows
    with pytest.raises(ValueError):
        d_optimal_rows(ss, 0)


# --- optimal_schedule: bound to a record ------------------------------------


def test_d_efficiency_at_least_one():
    """The reported optimal is the better of greedy/uniform, so it is never less informative
    than uniform — D-efficiency >= 1 by construction."""
    od = optimal_schedule(onkos.load(), CLARET, context=NSCLC, n_samples=7, horizon=48.0)
    assert od.d_efficiency >= 1.0 - 1e-9
    assert od.optimal.log_det_fim >= od.uniform.log_det_fim - 1e-9


def test_budget_and_baseline():
    od = optimal_schedule(onkos.load(), CLARET, context=NSCLC, n_samples=7, horizon=48.0)
    assert len(od.optimal.schedule) == 7
    assert all(0.0 <= t <= 48.0 for t in od.optimal.schedule)
    assert od.optimal.schedule[0] == 0.0  # baseline anchored (include_baseline default)


def test_structural_flat_survives_the_best_design():
    """The deeply flat growth rate kL stays unidentifiable even under the D-optimal design —
    a structural flat, not a design failure (the robust verdict; the borderline lambda rescue
    is budget-sensitive and not asserted)."""
    od = optimal_schedule(onkos.load(), CLARET, context=NSCLC, n_samples=7, horizon=48.0)
    assert "kL" in od.structurally_flat
    assert od.optimal.rse_percent["kL"] >= 50.0
    # the optimal design still tightens it (more informative than uniform).
    assert od.optimal.rse_percent["kL"] <= od.uniform.rse_percent["kL"]
    # kD is comfortably identifiable under both.
    assert od.optimal.rse_percent["kD"] < 50.0


def test_method_works_when_the_model_allows():
    """The 2-parameter biexponential is fully identifiable under the optimal design — the
    control proving the kL failure is the model's structure, not the optimizer."""
    od = optimal_schedule(onkos.load(), BIEXP, context=NSCLC, n_samples=7, horizon=48.0)
    assert od.structurally_flat == []
    assert od.d_efficiency >= 1.0 - 1e-9


def test_tier_passthrough_and_clinical_flag():
    od = optimal_schedule(onkos.load(), CLARET, context=NSCLC, n_samples=7, horizon=48.0)
    assert od.tier == "C"  # the record's propagated tier — design cannot move it
    d = od.to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "PROHIBITED" in d["onkos:clinicalUse"]


def test_non_ode_record_raises():
    """A survival link is not a trajectory to design for."""
    with pytest.raises(ValueError):
        optimal_schedule(onkos.load(), "survival_link.nsclc_os_week8", context=NSCLC)


def test_underdetermined_budget_raises():
    """A budget below the parameter count is rank-deficient and rejected."""
    with pytest.raises(ValueError):
        optimal_schedule(onkos.load(), CLARET, context=NSCLC, n_samples=2)  # 3 params
