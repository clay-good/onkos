"""Onkos — curated, citation-backed, tier-annotated TGI / survival dataset.

NOT a clinical decision tool. NOT a prognostic calculator. NOT a treatment
recommender. For drug-development methodology, simulation, and education only.
"""

from __future__ import annotations

from ._const import VERSION as __version__
from .audit import audit_tiers, evidence_ceiling
from .compare import Comparison, compare
from .filter import filter_records
from .load import DATASET_VERSION, Dataset, load
from .report import build_report
from .sensitivity import SensitivityResult, sensitivity
from .simulate import CLINICAL_USE, Trajectory, simulate
from .uncertainty import Ensemble, simulate_ensemble
from .validate import validate_dataset

__all__ = [
    "__version__",
    "DATASET_VERSION",
    "CLINICAL_USE",
    "Dataset",
    "Trajectory",
    "Comparison",
    "Ensemble",
    "SensitivityResult",
    "load",
    "simulate",
    "simulate_ensemble",
    "sensitivity",
    "compare",
    "filter_records",
    "validate_dataset",
    "build_report",
    "audit_tiers",
    "evidence_ceiling",
]
