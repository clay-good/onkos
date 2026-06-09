"""Evidence-based tier audit (spec §5: "tier assignment is partly numeric")."""

import copy

import onkos
from onkos.audit import audit_tiers, evidence_ceiling, inflated_records, is_evidence_tiered
from onkos.models import TIER_ORDER


def test_shipped_dataset_has_no_tier_inflation():
    assert inflated_records(onkos.load()) == []


def test_only_clinical_tgi_and_survival_are_audited():
    ds = onkos.load()
    audited = {f.record_id for f in audit_tiers(ds)}
    for r in ds:
        if r.subsystem in ("preclinical_translation", "immuno_oncology") or r.purpose in (
            "growth", "exposure_response", "translation",
        ):
            assert r.id not in audited


def test_high_cv_resistance_term_caps_ceiling_at_C():
    # Claret's resistance term has ~96% CV -> ceiling C even with an external check.
    ceiling = evidence_ceiling(onkos.load()["resistance.claret_2009.tgi"])
    assert ceiling == "C"


def test_well_identified_external_model_ceiling_is_B():
    # Survival links carry an external C-index and no IIV -> evidence supports B.
    ceiling = evidence_ceiling(onkos.load()["survival_link.nsclc_os_week8"])
    assert ceiling == "B"


def test_no_external_validation_ceiling_is_C():
    ds = onkos.load()
    r = copy.deepcopy(ds["survival_link.nsclc_os_week8"])
    r.predictive_performance = []  # strip the external metric (Record is mutable)
    assert evidence_ceiling(r) == "C"


def test_validate_rejects_tier_inflation(tmp_path):
    """A clinical record claiming a tier above its evidence ceiling fails validation."""
    import json

    from onkos._data import dataset_dir
    from onkos.validate import validate_dataset

    base = dataset_dir()
    rid = "survival_link.nsclc_os_week8"  # ceiling B
    rec = json.loads((base / "records" / f"{rid}.json").read_text())
    rec["tier"] = "A"  # inflate above the ceiling

    (tmp_path / "records").mkdir()
    (tmp_path / "schema").mkdir()
    (tmp_path / "citations").mkdir()
    (tmp_path / "schema" / "record.schema.json").write_text(
        (base / "schema" / "record.schema.json").read_text()
    )
    (tmp_path / "records" / f"{rid}.json").write_text(json.dumps(rec))
    (tmp_path / "citations" / "wang-2009-tgi.json").write_text(
        (base / "citations" / "wang-2009-tgi.json").read_text()
    )
    errs = validate_dataset(str(tmp_path))
    assert any("exceeds the evidence ceiling" in e for e in errs)


def test_ceiling_never_better_than_assigned_in_dataset():
    for f in audit_tiers(onkos.load()):
        assert TIER_ORDER[f.assigned] >= TIER_ORDER[f.ceiling]  # assigned no better than ceiling


def test_is_evidence_tiered_predicate():
    ds = onkos.load()
    assert is_evidence_tiered(ds["resistance.claret_2009.tgi"])
    assert not is_evidence_tiered(ds["growth_laws.gompertz"])
    assert not is_evidence_tiered(ds["immuno_oncology.kuznetsov_1994.tumor_immune"])


def test_cli_audit(capsys):
    from onkos.cli import main

    assert main(["audit"]) == 0
    out = capsys.readouterr().out
    assert "tier audit" in out and "inflated" in out
