"""Spec-completeness contract.

Guards the project invariant that a spec is never *written but unbuilt*. Every
forward-looking research spec under ``docs/specs/research/`` must be registered
here with (a) its landmark test file and (b) its shipped implementation
artifact — a module in the ``onkos`` package or a reference kernel. Adding a new
spec without an implementation, or deleting an implementation a spec depends on,
makes this test fail.

It also pins the two *enumerated* requirement sets in the v0.1 design spec: the
§2 growth-law family and the §7 export-format layer. Together these make "all
specs are complete" a checked condition rather than a claim.
"""

import importlib
from pathlib import Path

from onkos.export.reference import KERNELS

REPO = Path(__file__).resolve().parent.parent
SPEC_DIR = REPO / "docs" / "specs" / "research"
TESTS = REPO / "tests"

# research-spec stem -> (landmark test file, artifact kind, artifact name)
#   kind "module": ``onkos.<name>`` must import
#   kind "kernel": ``<name>`` must be a registered reference kernel
SPEC_INDEX = {
    "model-selection-uncertainty": ("test_combine.py", "module", "combine"),
    "practical-identifiability": ("test_identifiability.py", "module", "identify"),
    "combination-interaction": ("test_interaction.py", "module", "interaction"),
    "mechanistic-resistance": ("test_two_population.py", "kernel", "two_population_resistance"),
    "survival-metric-choice": ("test_survival_metric.py", "module", "simulate"),
    "model-selection-budget": ("test_budget.py", "module", "budget"),
    "recist-orr-surrogate": ("test_response.py", "module", "response"),
    "duration-of-response": ("test_duration.py", "module", "response"),
    "cross-context-generalization": ("test_context_library.py", "module", "compare"),
    "pfs-endpoint": ("test_pfs_routes.py", "module", "response"),
    "optimal-design": ("test_design.py", "module", "design"),
    "acquired-resistance": ("test_acquired.py", "kernel", "acquired_resistance"),
    "burden-auc-bridge-metric": ("test_burden_auc.py", "module", "metrics"),
    "joint-survival": ("test_joint.py", "module", "joint"),
    "loewe-additivity": ("test_loewe.py", "module", "interaction"),
    "exposure-response-extrapolation": ("test_dose_response.py", "module", "dose_response"),
    "early-surrogate-timing": ("test_early_surrogate.py", "module", "early_surrogate"),
    "model-discriminability": ("test_discriminability.py", "module", "discriminability"),
    "model-selection-atlas": ("test_atlas.py", "module", "atlas"),
    "von-bertalanffy-growth": ("test_landmarks.py", "kernel", "growth_von_bertalanffy"),
    "power-law-growth": ("test_landmarks.py", "kernel", "growth_power_law"),
}

# §2 — the declared in-scope unperturbed growth-law family.
GROWTH_LAW_KERNELS = {
    "growth_exponential",
    "growth_logistic",
    "growth_gompertz",
    "simeoni_exp_linear",
    "growth_von_bertalanffy",
    "growth_power_law",
}

# §7 — the declared export interop layer.
EXPORT_FORMATS = {
    "nonmem", "sbml", "pharmml", "so", "rxode2", "pumas", "vt-json",
    "jsonld", "omex", "csv", "bibtex",
}


def test_every_research_spec_is_indexed():
    """Every spec file on disk is accounted for, and vice versa (a bijection)."""
    on_disk = {p.stem for p in SPEC_DIR.glob("*.md")}
    indexed = set(SPEC_INDEX)
    assert not (on_disk - indexed), (
        f"research specs with no completeness entry (written but unindexed): "
        f"{sorted(on_disk - indexed)}"
    )
    assert not (indexed - on_disk), (
        f"indexed specs with no file on disk: {sorted(indexed - on_disk)}"
    )


def test_every_spec_has_its_landmark_test():
    for stem, (test_file, _, _) in SPEC_INDEX.items():
        assert (TESTS / test_file).exists(), f"{stem}: missing landmark test {test_file}"


def test_every_spec_has_a_shipped_artifact():
    for stem, (_, kind, name) in SPEC_INDEX.items():
        if kind == "module":
            importlib.import_module(f"onkos.{name}")  # ImportError if unbuilt
        elif kind == "kernel":
            assert name in KERNELS, f"{stem}: reference kernel '{name}' not registered"
        else:  # pragma: no cover - guards a typo in the index
            raise AssertionError(f"{stem}: unknown artifact kind {kind!r}")


def test_growth_law_family_complete():
    """Spec §2 enumerates the in-scope growth laws; all must be implemented."""
    missing = GROWTH_LAW_KERNELS - set(KERNELS)
    assert not missing, f"spec §2 growth laws not implemented: {sorted(missing)}"


def test_export_layer_complete():
    """Spec §7 enumerates the export interop layer; all formats must be available."""
    from onkos.cli import _TEXT_EXPORTERS

    available = set(_TEXT_EXPORTERS) | {"omex", "csv", "bibtex"}
    missing = EXPORT_FORMATS - available
    assert not missing, f"spec §7 export formats missing: {sorted(missing)}"
