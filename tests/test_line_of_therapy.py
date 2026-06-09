"""Line-of-therapy dimension: second-line context + line-aware survival matching."""

import numpy as np
import onkos
from onkos.filter import filter_records
from onkos.simulate import _find_survival_links

T = np.linspace(0.0, 208.0, 417)


def _ctx_endpoints(ds, line):
    out = set()
    for r in filter_records(ds, purpose="survival_link"):
        dc = r.derivation_context
        if dc and dc.tumor_type == "NSCLC" and dc.line_of_therapy == line:
            out.add(r.structure.get("endpoint", "OS"))
    return out


def test_nsclc_has_first_and_second_line_links():
    ds = onkos.load()
    assert _ctx_endpoints(ds, "first") == {"OS", "PFS"}
    assert _ctx_endpoints(ds, "second") == {"OS", "PFS"}


def test_survival_matching_is_line_aware():
    ds = onkos.load()
    first = {r.id for r in _find_survival_links(ds, "NSCLC", "first")}
    second = {r.id for r in _find_survival_links(ds, "NSCLC", "second")}
    assert first and second
    assert first.isdisjoint(second)  # no shared link across lines (no silent borrow)
    assert all("2l" in rid or "_pfs_week8" in rid or "_os_week8" in rid for rid in second)


def test_second_line_compare_has_two_models():
    ds = onkos.load()
    cmp = onkos.compare(ds, purpose="tgi", context={"tumor_type": "NSCLC", "line": "second"},
                        drug_effect=1.0, t=T)
    ids = {tr.record_id for tr in cmp.included}
    assert len(cmp.included) >= 2
    # the first-line-only Claret model is excluded for second line
    assert "resistance.claret_2009.tgi" not in ids
    assert "resistance.nsclc_second_line.claret" in ids


def test_second_line_survival_is_shorter_than_first():
    ds = onkos.load()
    # same model (validated for both lines), different line -> shorter 2L survival
    first = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential",
                           context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=1.0, t=T)
    second = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential",
                            context={"tumor_type": "NSCLC", "line": "second"}, drug_effect=1.0, t=T)
    assert second.median_os < first.median_os
    assert second.median_pfs < first.median_pfs


def test_unsupported_line_gets_no_survival_curve():
    ds = onkos.load()
    # no third-line links exist -> honest empty survival, not a borrowed one
    tr = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential",
                        context={"tumor_type": "NSCLC", "line": "third"}, drug_effect=1.0, t=T)
    assert tr.survival == {}


def test_second_line_baseline_is_more_advanced():
    ds = onkos.load()
    b1 = ds["tumor_type_baselines.nsclc_first_line"]["baseline_tumor_size"].central
    b2 = ds["tumor_type_baselines.nsclc_second_line"]["baseline_tumor_size"].central
    assert b2 > b1
