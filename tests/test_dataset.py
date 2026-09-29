"""Dataset integrity: schema validity, referential integrity, invariants."""

import json
import shutil

import onkos
from onkos._data import dataset_dir
from onkos.validate import validate_dataset


def test_dataset_validates():
    assert validate_dataset() == []


def test_load_nonempty():
    ds = onkos.load()
    assert len(ds) >= 7
    assert "resistance.claret_2009.tgi" in ds


def test_record_tier_is_worst_parameter_tier():
    ds = onkos.load()
    for r in ds:
        if not r.parameters:
            continue
        worst = max(p.tier for p in r.parameters)  # 'A' < 'B' < 'C' < 'D' lexicographically
        assert r.tier >= worst, f"{r.id}: record tier {r.tier} better than worst param {worst}"


def test_every_record_has_resolvable_primary_citation():
    ds = onkos.load()
    for r in ds:
        assert r.primary_citation is not None, r.id
        assert r.primary_citation.doi, r.id


def test_high_uncertainty_terms_carry_iiv():
    ds = onkos.load()
    claret = ds["resistance.claret_2009.tgi"]
    assert claret["lambda"].iiv_cv_percent is not None
    assert claret["lambda"].iiv_cv_percent > 50  # resistance term is poorly identified


def test_parameter_access_by_symbol():
    ds = onkos.load()
    claret = ds["resistance.claret_2009.tgi"]
    assert claret["kL"].central > 0
    assert "lambda" in claret


def _minimal_dataset(tmp_path):
    base = dataset_dir()
    for name in ("records", "schema", "citations"):
        (tmp_path / name).mkdir()
    shutil.copy2(base / "schema" / "record.schema.json", tmp_path / "schema")
    shutil.copy2(
        base / "records" / "growth_laws.exponential.json", tmp_path / "records"
    )
    citation = base / "citations" / "stein-2008-grc.json"
    shutil.copy2(citation, tmp_path / "citations")
    return citation


def test_validate_rejects_malformed_citation_json(tmp_path):
    citation = _minimal_dataset(tmp_path)
    (tmp_path / "citations" / citation.name).write_text("{")

    errors = validate_dataset(str(tmp_path))

    assert any("invalid citation JSON" in error for error in errors)


def test_validate_rejects_non_object_citation(tmp_path):
    citation = _minimal_dataset(tmp_path)
    (tmp_path / "citations" / citation.name).write_text("[]")

    errors = validate_dataset(str(tmp_path))

    assert any("citation must be a JSON object" in error for error in errors)


def test_validate_rejects_citation_filename_key_mismatch(tmp_path):
    citation = _minimal_dataset(tmp_path)
    data = json.loads(citation.read_text())
    data["key"] = "different-key"
    (tmp_path / "citations" / citation.name).write_text(json.dumps(data))

    errors = validate_dataset(str(tmp_path))

    assert any("filename does not match citation key" in error for error in errors)
    assert any("unknown citation key 'stein-2008-grc'" in error for error in errors)


def test_validate_rejects_duplicate_citation_keys(tmp_path):
    citation = _minimal_dataset(tmp_path)
    shutil.copy2(citation, tmp_path / "citations" / "duplicate.json")

    errors = validate_dataset(str(tmp_path))

    assert any("duplicate citation key 'stein-2008-grc'" in error for error in errors)
