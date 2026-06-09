"""JSON-Schema validation of the dataset (structural integrity gate)."""

from __future__ import annotations

import json
from pathlib import Path

from jsonschema import Draft7Validator

from ._data import dataset_dir


def _schema(base: Path) -> dict:
    return json.loads((base / "schema" / "record.schema.json").read_text())


def validate_dataset(path: str | None = None) -> list[str]:
    """Validate every record file against the schema and check referential
    integrity (citation keys, kernel names, unique ids).

    Returns a list of human-readable error strings; empty means valid.
    """
    base = Path(path) if path else dataset_dir()
    validator = Draft7Validator(_schema(base))
    errors: list[str] = []

    known_citations = {fp.stem for fp in (base / "citations").glob("*.json")}
    seen_ids = set()

    # Imported lazily to avoid a hard dependency cycle at module import time.
    from .export.reference import KERNELS

    for fp in sorted((base / "records").glob("*.json")):
        try:
            record = json.loads(fp.read_text())
        except json.JSONDecodeError as exc:
            errors.append(f"{fp.name}: invalid JSON ({exc})")
            continue

        for err in sorted(validator.iter_errors(record), key=lambda e: e.path):
            loc = "/".join(str(p) for p in err.path) or "<root>"
            errors.append(f"{fp.name}: {loc}: {err.message}")

        rid = record.get("id")
        if rid in seen_ids:
            errors.append(f"{fp.name}: duplicate id '{rid}'")
        seen_ids.add(rid)
        if rid and fp.stem != rid:
            errors.append(f"{fp.name}: filename does not match id '{rid}'")

        for key in _citation_keys(record):
            if key not in known_citations:
                errors.append(f"{fp.name}: unknown citation key '{key}'")

        kernel = record.get("kernel")
        if kernel is not None and kernel not in KERNELS:
            errors.append(f"{fp.name}: unknown kernel '{kernel}'")

        # Hypothesis-tier safety rule (spec §5, §10): immuno-oncology records ship
        # tier D, non-predictive — at the record level and every parameter.
        if record.get("subsystem") == "immuno_oncology":
            if record.get("tier") != "D":
                errors.append(
                    f"{fp.name}: immuno_oncology records must be tier D (hypothesis-tier), "
                    f"got '{record.get('tier')}'"
                )
            for p in record.get("parameters", []):
                if p.get("tier") != "D":
                    errors.append(
                        f"{fp.name}: immuno_oncology parameter '{p.get('symbol')}' must be tier D"
                    )

    # Evidence-based tier audit (spec §5, §9): a clinical TGI / survival record may
    # not claim a better tier than its recorded external validation supports.
    if not errors:
        from .audit import inflated_records
        from .load import load

        for f in inflated_records(load(str(base))):
            errors.append(
                f"{f.record_id}: tier '{f.assigned}' exceeds the evidence ceiling '{f.ceiling}' "
                f"(external_validation={f.has_external}, max_iiv_cv={f.max_iiv:.0f}, "
                f"validated_tumor_types={f.breadth})"
            )

    return errors


def _citation_keys(record: dict):
    if record.get("primary_citation"):
        yield record["primary_citation"]
    for p in record.get("parameters", []):
        if p.get("primary_citation"):
            yield p["primary_citation"]
    for pp in record.get("predictive_performance", []):
        if pp.get("citation"):
            yield pp["citation"]
    for fm in record.get("known_failure_modes", []):
        if fm.get("citation"):
            yield fm["citation"]
