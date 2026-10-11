"""Dataset validation: JSON Schema (structure), referential integrity, and semantic
integrity (every record must be internally consistent and actually computable)."""

from __future__ import annotations

import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
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

    known_citations: set[str] = set()
    for fp in sorted((base / "citations").glob("*.json")):
        try:
            citation = json.loads(fp.read_text())
        except json.JSONDecodeError as exc:
            errors.append(f"{fp.name}: invalid citation JSON ({exc})")
            continue
        if not isinstance(citation, dict):
            errors.append(f"{fp.name}: citation must be a JSON object")
            continue

        key = citation.get("key")
        if not isinstance(key, str) or not key:
            errors.append(f"{fp.name}: citation key must be a non-empty string")
            continue
        if fp.stem != key:
            errors.append(f"{fp.name}: filename does not match citation key '{key}'")
        if key in known_citations:
            errors.append(f"{fp.name}: duplicate citation key '{key}'")
        known_citations.add(key)

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

    # Semantic integrity: the schema proves a record is well-formed, not that it is
    # consistent or computable. Needs a loadable dataset, so it runs only once the
    # structural checks pass.
    if not errors:
        from .load import load

        errors.extend(_semantic_errors(load(str(base))))

    # Evidence-based tier audit (spec §5, §9): a clinical TGI / survival record may
    # not claim a better tier than its recorded external validation supports.
    if not errors:
        from .audit import inflated_records

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


def _semantic_errors(ds) -> list[str]:
    """Per-record and cross-record consistency checks the JSON Schema cannot express."""
    from .export.reference import KERNELS
    from .metrics import extract_tgi_metrics

    metric_keys = set(extract_tgi_metrics(np.array([0.0, 1.0]), np.array([1.0, 1.0]), 1.0))
    errors: list[str] = []
    default_links: Counter = Counter()
    baselines: Counter = Counter()

    for r in ds:
        def err(msg: str, rid: str = r.id) -> None:
            errors.append(f"{rid}: {msg}")

        symbols = [p.symbol for p in r.parameters]
        for sym, count in Counter(symbols).items():
            if count > 1:
                err(f"parameter '{sym}' is declared {count} times")
        for p in r.parameters:
            v = p.value
            if not v.units.strip():
                err(f"parameter '{p.symbol}' has empty units")
            # Every curated quantity (rate, size, capacity, potency, shape, scale, hazard
            # coefficient) is non-negative by the dataset's sign conventions.
            if v.central < 0:
                err(f"parameter '{p.symbol}' central value {v.central} is negative")
            if (v.low is not None and v.low > v.central) or (
                v.high is not None and v.high < v.central
            ):
                err(f"parameter '{p.symbol}' central {v.central} is outside [low, high] = "
                    f"[{v.low}, {v.high}]")
            if p.iiv_cv_percent is not None and p.iiv_cv_percent <= 0:
                err(f"parameter '{p.symbol}' iiv_cv_percent must be positive "
                    f"(got {p.iiv_cv_percent}); use null for 'not reported'")

        # A record's own derivation context must sit inside the envelope it declares
        # validated — otherwise it would be flagged out-of-context where it was derived.
        dc, tp = r.derivation_context, r.transportability
        if dc and tp:
            if dc.tumor_type and tp.validated_tumor_types and (
                dc.tumor_type not in tp.validated_tumor_types
            ):
                err(f"derivation tumor_type '{dc.tumor_type}' is not in its own "
                    f"validated_tumor_types {list(tp.validated_tumor_types)}")
            if dc.line_of_therapy and tp.validated_lines and (
                dc.line_of_therapy not in tp.validated_lines
            ):
                err(f"derivation line '{dc.line_of_therapy}' is not in its own "
                    f"validated_lines {list(tp.validated_lines)}")

        if r.kind == "context_baseline":
            # simulate() silently falls back to y0 = 100 when this is missing.
            if "baseline_tumor_size" not in symbols:
                err("context_baseline record has no 'baseline_tumor_size' parameter")
            elif r["baseline_tumor_size"].central <= 0:
                err("baseline_tumor_size must be positive")
            if dc:
                baselines[(dc.tumor_type, dc.line_of_therapy)] += 1
            continue

        if r.kernel not in KERNELS:
            err("model record has no bound kernel")
            continue
        spec = KERNELS[r.kernel]
        missing = [x for x in spec.record_symbols if x not in symbols]
        extra = [x for x in symbols if x not in spec.record_symbols]
        if missing:
            err(f"kernel '{spec.name}' requires parameter(s) {missing}")
        if extra:
            err(f"parameter(s) {extra} are not used by kernel '{spec.name}' "
                f"(expects {list(spec.record_symbols)})")

        if spec.kind == "survival":
            if dc and r.structure.get("default", True):
                default_links[
                    (dc.tumor_type, dc.line_of_therapy, r.structure.get("endpoint", "OS"))
                ] += 1
            metric = r.structure.get("link_metric")
            if metric is not None and metric not in metric_keys:
                err(f"link_metric '{metric}' is not a TGI metric (choose from "
                    f"{sorted(metric_keys)})")
            if spec.uses_baseline:
                table_errors = _baseline_table_errors(r.structure.get("baseline_survival"))
                for msg in table_errors:
                    err(msg)
                if table_errors:
                    continue

        if not missing and len(set(symbols)) == len(symbols):
            problem = _computation_problem(ds, r, spec)
            if problem:
                err(f"does not compute at its central values: {problem}")

    for (tumor_type, line, endpoint), count in default_links.items():
        if count > 1:
            errors.append(
                f"{count} default {endpoint} survival links for {tumor_type}/{line}: the default "
                f"must be unique (mark the alternatives structure.default = false)"
            )
    for (tumor_type, line), count in baselines.items():
        if count > 1:
            errors.append(f"{count} context_baseline records for {tumor_type}/{line}")
    return errors


def _baseline_table_errors(table) -> list[str]:
    """A tabulated (Cox) baseline must be a survival function: S(0) = 1, non-increasing,
    within [0, 1], on strictly increasing times from 0."""
    if not table:
        return ["kernel needs structure.baseline_survival (a tabulated baseline)"]
    times = np.asarray(table.get("times", []), dtype=float)
    surv = np.asarray(table.get("survival", []), dtype=float)
    if times.size != surv.size or times.size < 2:
        return ["baseline_survival times and survival must have equal length >= 2"]
    out = []
    if times[0] != 0 or np.any(np.diff(times) <= 0):
        out.append("baseline_survival times must start at 0 and strictly increase")
    if surv[0] != 1 or np.any(surv < 0) or np.any(np.diff(surv) > 0):
        out.append("baseline_survival must start at 1, stay >= 0, and never increase")
    return out


def _computation_problem(ds, record, spec) -> str | None:
    """Evaluate the record's kernel at its central values and return what is wrong with
    the output (``None`` if it is a valid trajectory / survival function / effect)."""
    from .export.reference import effect
    from .export.registry import kernel_values
    from .simulate import _survival_vals, simulate

    t = np.linspace(0.0, 104.0, 105)
    try:
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            if spec.kind == "ode":
                y = simulate(ds, record.id, context={"y0": 100.0}, t=t).tumor_size
                if not np.all(np.isfinite(y)) or np.any(y < -1e-6):
                    return "tumor size is not finite and non-negative"
            elif spec.kind == "survival":
                for x in (-0.5, 0.0, 0.5):
                    s = np.asarray(spec.analytic(t, _survival_vals(record, x)), dtype=float)
                    if (
                        not np.all(np.isfinite(s))
                        or abs(s[0] - 1.0) > 1e-9
                        or np.any(s < 0)
                        or np.any(np.diff(s) > 1e-12)
                    ):
                        return "survival curve is not a valid survival function"
            else:
                e = np.asarray(effect(spec, np.array([0.0, 1.0, 100.0]), kernel_values(record)))
                if not np.all(np.isfinite(e)) or np.any(e < 0):
                    return "effect is not finite and non-negative"
    except (ArithmeticError, ValueError, RuntimeError) as exc:
        return f"{type(exc).__name__}: {exc}"
    return None
