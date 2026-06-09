"""Phase E — immuno-oncology QSP (hypothesis-tier, NOT FOR PREDICTION)."""

import json
from pathlib import Path

import numpy as np
import onkos
from onkos._data import dataset_dir
from onkos.export import to_nonmem, to_sbml, to_virtual_trial_json
from onkos.export.annotate import PREDICTION_PROHIBITED

CANONICAL = "immuno_oncology.kuznetsov_1994.tumor_immune"
POORLY = "immuno_oncology.poorly_immunogenic.hypothesis"
T = np.linspace(0, 200, 801)


def _final(rid, E):
    tr = onkos.simulate(onkos.load(), rid, context={"tumor_type": "x", "y0": 10.0}, drug_effect=E, t=T)
    return float(tr.tumor_size[-1])


def test_io_records_are_hypothesis_tier_D():
    ds = onkos.load()
    io = [r for r in ds if r.subsystem == "immuno_oncology"]
    assert len(io) >= 2
    for r in io:
        assert r.tier == "D", f"{r.id}: must be tier D"
        assert all(p.tier == "D" for p in r.parameters), f"{r.id}: all params must be tier D"


def test_immune_control_vs_escape():
    # Canonical (immunogenic) controls the tumor; the poorly-immunogenic variant escapes.
    assert _final(CANONICAL, 0.0) < 100.0          # controlled, far below capacity (1/beta=500)
    assert _final(POORLY, 0.0) > 400.0             # escape, near carrying capacity


def test_checkpoint_blockade_rescues_escape():
    # An immunotherapy effect above the bistable threshold flips escape -> control.
    assert _final(POORLY, 12.0) < 0.25 * _final(POORLY, 0.0)


def test_io_excluded_from_clinical_compare():
    ds = onkos.load()
    cmp = onkos.compare(ds, purpose="tgi", context={"tumor_type": "NSCLC", "line": "first"})
    assert not any(t.record_id.startswith("immuno_oncology") for t in cmp.included)
    assert not any(rid.startswith("immuno_oncology") for rid, _ in cmp.excluded)


def test_io_has_no_survival_curve():
    ds = onkos.load()
    tr = onkos.simulate(ds, CANONICAL, context={"tumor_type": "immunogenic_tumor", "y0": 10.0})
    assert tr.os_curve is None


def test_do_not_use_for_prediction_in_exports():
    ds = onkos.load()
    r = ds[CANONICAL]
    assert PREDICTION_PROHIBITED in to_sbml(r)
    assert PREDICTION_PROHIBITED in to_nonmem(r)
    vt = json.loads(to_virtual_trial_json(r))
    assert vt["DO_NOT_USE_FOR_PREDICTION"] is True
    assert vt["onkos:predictionStatus"] == PREDICTION_PROHIBITED
    # a clinical (non-hypothesis) record must NOT carry the prediction prohibition
    clin = json.loads(to_virtual_trial_json(ds["resistance.claret_2009.tgi"]))
    assert clin["DO_NOT_USE_FOR_PREDICTION"] is False


def test_validate_enforces_tier_D_for_io():
    """The validator rejects an immuno_oncology record that is not tier D."""
    from onkos.validate import validate_dataset

    assert validate_dataset() == []  # the shipped dataset is clean

    # construct a temp dataset dir with a mis-tiered IO record and confirm it fails
    base = dataset_dir()
    bad = json.loads((base / "records" / f"{CANONICAL}.json").read_text())
    bad["tier"] = "B"
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "records").mkdir()
        (d / "schema").mkdir()
        (d / "citations").mkdir()
        (d / "schema" / "record.schema.json").write_text(
            (base / "schema" / "record.schema.json").read_text()
        )
        (d / "records" / f"{CANONICAL}.json").write_text(json.dumps(bad))
        (d / "citations" / "kuznetsov-1994-tumor-immune.json").write_text(
            (base / "citations" / "kuznetsov-1994-tumor-immune.json").read_text()
        )
        errs = validate_dataset(str(d))
        assert any("must be tier D" in e for e in errs)
