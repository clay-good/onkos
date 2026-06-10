"""Onkos — curated, citation-backed, tier-annotated TGI / survival dataset.

NOT a clinical decision tool. NOT a prognostic calculator. NOT a treatment
recommender. For drug-development methodology, simulation, and education only.
"""

from __future__ import annotations

from . import pk
from ._const import VERSION as __version__
from .audit import audit_tiers, evidence_ceiling
from .budget import Budget, model_selection_budget
from .combine import ModelAverage
from .compare import Comparison, compare
from .design import OptimalDesign, optimal_schedule
from .filter import filter_records
from .identify import Identifiability, identifiability
from .interaction import (
    InteractionComparison,
    combine_effects,
    compare_interactions,
    simulate_combination,
)
from .load import DATASET_VERSION, Dataset, load
from .report import build_report
from .response import (
    PFSRouteDivergence,
    ProgressionFreeSurvival,
    ResponseRates,
    ResponseSurvival,
    objective_response_rate,
    pfs_route_divergence,
    progression_free_survival,
    response_vs_survival,
)
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
    "ModelAverage",
    "SensitivityResult",
    "Identifiability",
    "OptimalDesign",
    "optimal_schedule",
    "InteractionComparison",
    "Budget",
    "model_selection_budget",
    "ResponseRates",
    "ResponseSurvival",
    "objective_response_rate",
    "response_vs_survival",
    "ProgressionFreeSurvival",
    "PFSRouteDivergence",
    "progression_free_survival",
    "pfs_route_divergence",
    "load",
    "simulate",
    "simulate_ensemble",
    "sensitivity",
    "identifiability",
    "combine_effects",
    "simulate_combination",
    "compare_interactions",
    "compare",
    "filter_records",
    "validate_dataset",
    "build_report",
    "audit_tiers",
    "evidence_ceiling",
    "pk",
]
