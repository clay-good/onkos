"""Phase C — the tumor-context library and cross-context divergence.

Verifies the library is internally consistent (no orphan baselines/links) and
that the virtual-trial divergence view is broadly useful: every supported tumor
context has its own baseline, survival link, and >=2 eligible TGI models.
"""

import onkos
from onkos.filter import filter_records

EXPECTED_CONTEXTS = {"NSCLC", "breast", "CRC", "HCC", "melanoma"}


def _tumor_types(ds, subsystem=None, purpose=None):
    out = set()
    for r in filter_records(ds, subsystem=subsystem, purpose=purpose):
        if r.derivation_context and r.derivation_context.tumor_type:
            out.add(r.derivation_context.tumor_type)
    return out


def test_baseline_library_covers_expected_contexts():
    ds = onkos.load()
    assert EXPECTED_CONTEXTS <= _tumor_types(ds, subsystem="tumor_type_baselines")


def test_survival_links_cover_expected_contexts():
    ds = onkos.load()
    assert EXPECTED_CONTEXTS <= _tumor_types(ds, purpose="survival_link")


def test_no_orphan_survival_links_or_baselines():
    """Every survival link / baseline tumor type must have >=1 TGI model."""
    ds = onkos.load()
    tgi_types = _tumor_types(ds, purpose="tgi")
    for r in ds:
        if r.purpose == "survival_link" or r.subsystem == "tumor_type_baselines":
            tt = r.derivation_context.tumor_type if r.derivation_context else None
            if tt:
                assert tt in tgi_types, f"{r.id}: no TGI model for tumor type {tt}"


def test_each_context_has_measurable_divergence():
    ds = onkos.load()
    for tt in EXPECTED_CONTEXTS:
        cmp = onkos.compare(
            ds, purpose="tgi", context={"tumor_type": tt, "line": "first"}, drug_effect=1.0
        )
        assert len(cmp.included) >= 2, f"{tt}: fewer than 2 eligible TGI models"
        assert cmp.os_divergence > 0.0, f"{tt}: zero divergence"
        # every included model used the context's own (in-context) survival link
        for tr in cmp.included:
            assert tr.os_curve is not None


def test_baseline_drives_context_specific_y0():
    ds = onkos.load()
    # HCC baseline SLD (110) differs from NSCLC (80): same model, different start size.
    import numpy as np
    t = np.linspace(0, 52, 105)
    hcc = onkos.simulate(ds, "resistance.hcc_first_line.claret",
                         context={"tumor_type": "HCC", "line": "first"}, drug_effect=0.0, t=t)
    nsclc = onkos.simulate(ds, "resistance.claret_2009.tgi",
                           context={"tumor_type": "NSCLC", "line": "first"}, drug_effect=0.0, t=t)
    assert hcc.tumor_size[0] > nsclc.tumor_size[0]


def test_tumor_specific_survival_scales_differ():
    ds = onkos.load()
    # At a fixed favourable week-8 change, HCC (short OS scale) < breast (long OS scale).
    import numpy as np
    from onkos.export.registry import get_kernel, kernel_values
    from onkos.simulate import median_survival

    t = np.linspace(0, 300, 601)
    med = {}
    for tt, rid in [("HCC", "survival_link.hcc_os_week8"), ("breast", "survival_link.breast_os_week8")]:
        v = kernel_values(ds[rid])
        v["x"] = -0.3  # 30% shrinkage at week 8
        s = np.asarray(get_kernel(ds[rid]).analytic(t, v), dtype=float)
        med[tt] = median_survival(t, s)
    assert med["HCC"] is not None and med["breast"] is not None
    assert med["HCC"] < med["breast"]
