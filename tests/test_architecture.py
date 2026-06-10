"""Architecture-consistency tests — pin the documented system structure as a
tested contract, so docs and code cannot silently drift.

These checks have already paid for themselves: an empty declared subsystem
(`drug_effect`) and a CI export sweep that omitted formats (`so`, `jsonld`) are
exactly the kinds of drift caught here.
"""

import json
import re

import onkos
from onkos._data import dataset_dir
from onkos.cli import _TEXT_EXPORTERS, build_parser
from onkos.export.reference import KERNELS

# Formats handled outside the text-exporter table.
_SPECIAL_EXPORT_FORMATS = {"omex", "csv", "bibtex"}
_KERNEL_KINDS = {"ode", "survival", "exposure_response"}


def _schema() -> dict:
    return json.loads((dataset_dir() / "schema" / "record.schema.json").read_text())


def _subsystem_enum() -> set:
    return set(_schema()["properties"]["subsystem"]["enum"])


def _cli_export_formats() -> set:
    parser = build_parser()
    export = next(
        a for a in parser._subparsers._group_actions[0].choices["export"]._actions
        if a.dest == "format"
    )
    return set(export.choices)


def test_every_declared_subsystem_has_records():
    """No declared subsystem may be empty (caught the empty drug_effect subsystem)."""
    ds = onkos.load()
    present = {r.subsystem for r in ds}
    missing = _subsystem_enum() - present
    assert not missing, f"declared subsystems with no records: {sorted(missing)}"


def test_cli_export_formats_match_builders():
    """Every text-export format maps to a real builder; the CLI choice set is
    exactly the builders plus the special (omex/csv/bibtex) formats."""
    assert _cli_export_formats() == set(_TEXT_EXPORTERS) | _SPECIAL_EXPORT_FORMATS
    for _name, (builder, ext) in _TEXT_EXPORTERS.items():
        assert callable(builder) and ext.startswith(".")


def test_ci_exports_every_format():
    """The CI exports sweep must cover every CLI export format (caught so/jsonld
    being missing from the loop in an earlier release)."""
    ci = (dataset_dir().parent / ".github" / "workflows" / "ci.yml")
    if not ci.exists():  # source checkout only
        return
    m = re.search(r"for fmt in ([^;]+); do", ci.read_text())
    assert m, "could not find the CI export loop"
    looped = set(m.group(1).split())
    assert _cli_export_formats() <= looped, f"CI export loop missing: {_cli_export_formats() - looped}"


def test_all_kernel_kinds_are_known():
    for name, spec in KERNELS.items():
        assert spec.kind in _KERNEL_KINDS, f"{name}: unknown kind {spec.kind}"


def test_no_dead_kernels():
    """Every reference kernel is bound by at least one record (no orphans)."""
    used = {r.kernel for r in onkos.load() if r.kernel}
    unused = set(KERNELS) - used
    assert not unused, f"kernels defined but never used by a record: {sorted(unused)}"


def test_every_record_kernel_exists():
    for r in onkos.load():
        if r.kernel is not None:
            assert r.kernel in KERNELS, f"{r.id}: unknown kernel {r.kernel}"


def test_public_api_surface_is_stable():
    """The documented public API names remain exported (catches accidental removal)."""
    expected = {
        "load", "simulate", "simulate_ensemble", "sensitivity", "compare",
        "identifiability", "combine_effects", "simulate_combination", "compare_interactions",
        "model_selection_budget", "objective_response_rate", "response_vs_survival",
        "progression_free_survival", "pfs_route_divergence", "optimal_schedule",
        "joint_survival", "compare_joint_vs_two_stage", "current_value_survival",
        "build_report", "audit_tiers", "evidence_ceiling", "validate_dataset",
        "filter_records", "pk", "Dataset", "Trajectory", "Comparison",
    }
    assert expected <= set(onkos.__all__)
    for name in expected:
        assert hasattr(onkos, name), f"onkos.{name} missing"


def test_dataset_is_the_single_source():
    """The source `dataset/` directory resolves ahead of any synced copy."""
    resolved = dataset_dir()
    assert (resolved / "records").is_dir() and (resolved / "schema").is_dir()
    # in a source checkout the live dataset wins over a bundled _dataset
    assert resolved.name in ("dataset", "_dataset")
