"""Serializable virtual-trial result (Comparison.to_dict/to_json) + dashboard."""

import json
import py_compile
from pathlib import Path

import numpy as np
import onkos
from onkos._data import dataset_dir

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


def _cmp():
    return onkos.compare(onkos.load(), purpose="tgi", context=NSCLC, drug_effect=1.0)


def test_to_dict_structure_and_endpoints():
    d = _cmp().to_dict()
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "PROHIBITED" in d["onkos:clinicalUse"]
    assert d["context"] == NSCLC
    assert d["n_included"] == len(d["included"]) >= 2
    for m in d["included"]:
        assert {"id", "tier", "median_os_weeks", "median_pfs_weeks", "warnings"} <= set(m)


def test_to_dict_divergences_match_properties():
    cmp = _cmp()
    d = cmp.to_dict()
    assert d["os_divergence"] == cmp.os_divergence
    assert d["pfs_divergence"] == cmp.pfs_divergence
    assert d["median_os_range"] == list(cmp.median_os_range)


def test_to_json_is_serializable_and_round_trips():
    s = _cmp().to_json()
    d = json.loads(s)
    assert isinstance(d["included"], list) and d["excluded"]  # NSCLC has excluded models


def test_include_curves_adds_arrays():
    d = _cmp().to_dict(include_curves=True)
    m = d["included"][0]
    assert len(m["tumor_size"]) == len(m["t"])
    assert "OS" in m["survival"] and "PFS" in m["survival"]
    # plain JSON types, not numpy
    assert isinstance(m["tumor_size"][0], float)
    json.dumps(d)  # must not raise


def test_no_curves_by_default_keeps_payload_small():
    m = _cmp().to_dict()["included"][0]
    assert "tumor_size" not in m and "survival" not in m


def test_curves_are_finite():
    d = _cmp().to_dict(include_curves=True)
    for m in d["included"]:
        assert np.all(np.isfinite(m["tumor_size"]))


def test_dashboard_compiles():
    """The Streamlit dashboard stays syntactically valid against the API."""
    app = dataset_dir().parent / "dashboard" / "app.py"
    if app.exists():  # source checkout only
        py_compile.compile(str(app), doraise=True)


def test_dashboard_uses_only_public_api():
    app = Path(dataset_dir().parent / "dashboard" / "app.py")
    if app.exists():
        src = app.read_text()
        # the dashboard must not reach into private simulate internals
        assert "_survival_vals" not in src and "_find_survival" not in src
