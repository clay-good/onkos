"""Onkos — curated, citation-backed, tier-annotated TGI / survival dataset.

NOT a clinical decision tool. NOT a prognostic calculator. NOT a treatment
recommender. For drug-development methodology, simulation, and education only.
"""

from __future__ import annotations

from .compare import Comparison, compare
from .filter import filter_records
from .load import DATASET_VERSION, Dataset, load
from .simulate import CLINICAL_USE, Trajectory, simulate
from .validate import validate_dataset

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "DATASET_VERSION",
    "CLINICAL_USE",
    "Dataset",
    "Trajectory",
    "Comparison",
    "load",
    "simulate",
    "compare",
    "filter_records",
    "validate_dataset",
]
