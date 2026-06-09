"""Progression-free survival (PFS) — the second survival endpoint."""

import numpy as np
import onkos
from onkos.filter import filter_records

NSCLC = {"tumor_type": "NSCLC", "line": "first"}
CONTEXTS = ["NSCLC", "breast", "CRC", "HCC", "melanoma"]


def _endpoints(ds, tumor):
    out = set()
    for r in filter_records(ds, purpose="survival_link"):
        if r.derivation_context and r.derivation_context.tumor_type == tumor:
            out.add(r.structure.get("endpoint", "OS"))
    return out


def test_every_context_has_os_and_pfs_links():
    ds = onkos.load()
    for tt in CONTEXTS:
        assert _endpoints(ds, tt) == {"OS", "PFS"}, tt


def test_simulation_produces_both_endpoints():
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0)
    assert set(tr.survival) == {"OS", "PFS"}
    assert tr.os_curve is not None and tr.pfs_curve is not None


def test_pfs_is_shorter_than_os():
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0)
    assert tr.median_pfs < tr.median_os


def test_pfs_curve_is_monotone_non_increasing():
    ds = onkos.load()
    tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=NSCLC, drug_effect=1.0)
    assert np.all(np.diff(tr.pfs_curve) <= 1e-9)


def test_compare_reports_pfs_divergence():
    ds = onkos.load()
    cmp = onkos.compare(ds, purpose="tgi", context=NSCLC, drug_effect=1.0)
    assert cmp.pfs_divergence > 0.0
    lo, hi = cmp.median_pfs_range
    assert hi > lo
    # every included model carries both endpoints
    assert all(set(tr.survival) == {"OS", "PFS"} for tr in cmp.included)


def test_deeper_response_lengthens_pfs():
    ds = onkos.load()
    weak = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential", context=NSCLC, drug_effect=0.5)
    strong = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential", context=NSCLC, drug_effect=1.5)
    assert strong.median_pfs > weak.median_pfs


def test_sensitivity_can_target_pfs():
    ds = onkos.load()
    res = onkos.sensitivity(ds, "resistance.claret_2009.tgi", context=NSCLC,
                            target="median_pfs_weeks", n=150, seed=0)
    assert res.dominant is not None and res.r_squared > 0
