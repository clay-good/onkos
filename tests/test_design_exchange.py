"""Fedorov exchange refinement of the D-optimal schedule.

Follow-on to ``docs/specs/research/optimal-design.md`` §8: greedy forward selection is a
heuristic; exchange is guaranteed no worse, ends at a local optimum under single swaps,
and recovers the brute-force optimum in most small cases where greedy misses it.
"""

import itertools

import numpy as np
import onkos
import pytest
from onkos.design import _logdet_fim, d_optimal_rows, exchange_refine, optimal_schedule


def _brute_force(ss, n, fixed=(0,)):
    free = [r for r in range(ss.shape[0]) if r not in fixed]
    return max((sorted(list(fixed) + list(c)) for c in itertools.combinations(free, n - len(fixed))),
               key=lambda rows: _logdet_fim(ss[rows], ss.shape[1]))


def test_exchange_never_loses_information_and_keeps_fixed_rows():
    for seed in range(40):
        ss = np.random.default_rng(seed).normal(size=(12, 3))
        greedy = sorted(d_optimal_rows(ss, 5, seed_rows=(0,)))
        refined = exchange_refine(ss, greedy, fixed_rows=(0,))
        assert 0 in refined and len(refined) == 5 == len(set(refined))
        assert _logdet_fim(ss[refined], 3) >= _logdet_fim(ss[greedy], 3) - 1e-12


def test_exchange_ends_at_a_local_optimum():
    ss = np.random.default_rng(7).normal(size=(10, 3))
    rows = exchange_refine(ss, sorted(d_optimal_rows(ss, 4, seed_rows=(0,))), fixed_rows=(0,))
    base = _logdet_fim(ss[rows], 3)
    for i, out in enumerate(rows):
        if out == 0:
            continue
        for cand in set(range(10)) - set(rows):
            swapped = rows[:i] + [cand] + rows[i + 1:]
            assert _logdet_fim(ss[swapped], 3) <= base + 1e-9


def test_exchange_recovers_the_optimum_greedy_misses():
    ss = np.random.default_rng(2).normal(size=(9, 3))  # a case where greedy is suboptimal
    greedy = sorted(d_optimal_rows(ss, 4, seed_rows=(0,)))
    best = _brute_force(ss, 4)
    assert _logdet_fim(ss[greedy], 3) < _logdet_fim(ss[best], 3) - 1e-6
    assert exchange_refine(ss, greedy, fixed_rows=(0,)) == best


def test_exchange_closes_most_greedy_gaps():
    gaps = recovered = 0
    for seed in range(150):
        ss = np.random.default_rng(seed).normal(size=(9, 3))
        greedy = sorted(d_optimal_rows(ss, 4, seed_rows=(0,)))
        best = _brute_force(ss, 4)
        if _logdet_fim(ss[greedy], 3) < _logdet_fim(ss[best], 3) - 1e-6:
            gaps += 1
            refined = exchange_refine(ss, greedy, fixed_rows=(0,))
            recovered += np.isclose(_logdet_fim(ss[refined], 3), _logdet_fim(ss[best], 3))
    assert gaps >= 5 and recovered / gaps >= 0.8


def test_refined_claret_design_is_at_least_as_informative_and_kl_stays_flat():
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    greedy = optimal_schedule(ds, "resistance.claret_2009.tgi", context=ctx, refine=False)
    refined = optimal_schedule(ds, "resistance.claret_2009.tgi", context=ctx)
    assert refined.optimal.log_det_fim >= greedy.optimal.log_det_fim
    assert refined.d_efficiency >= greedy.d_efficiency >= 1.0
    assert "kL" not in refined.optimal.identifiable  # structural flat, not a design failure
    assert any("Fedorov" in w for w in refined.warnings)


def test_exchange_rejects_malformed_input():
    with pytest.raises(ValueError):
        exchange_refine(np.ones(4), [0])
    with pytest.raises(ValueError):
        exchange_refine(np.ones((4, 2)), [])
    with pytest.raises(ValueError):
        exchange_refine(np.ones((4, 2)), [0, 9])
