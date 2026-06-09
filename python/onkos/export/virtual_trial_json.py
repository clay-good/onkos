"""virtual-trial JSON — parameters + tumor context + survival link for the
simulator/dashboard, carrying the mandatory NOT FOR CLINICAL USE flag."""

from __future__ import annotations

import json

from .._const import CLINICAL_USE
from .._const import VERSION as _V
from ..models import Record
from .annotate import PREDICTION_PROHIBITED, identifier_uris, is_hypothesis_tier
from .registry import get_kernel, kernel_values


def virtual_trial_dict(record: Record, *, tier=None, dataset_version: str = _V) -> dict:
    spec = get_kernel(record)
    dc = record.derivation_context
    tp = record.transportability
    return {
        # Minimal @context so the onkos:-prefixed keys below are valid JSON-LD.
        "@context": {
            "onkos": "https://onkos.dev/ns#",
            "dcterms": "http://purl.org/dc/terms/",
        },
        "@id": f"https://onkos.dev/record/{record.id}",
        "onkos:clinicalUse": CLINICAL_USE,
        "NOT_FOR_CLINICAL_USE": True,
        "DO_NOT_USE_FOR_PREDICTION": is_hypothesis_tier(record),
        "onkos:predictionStatus": PREDICTION_PROHIBITED if is_hypothesis_tier(record) else None,
        "datasetVersion": dataset_version,
        "id": record.id,
        "name": record.name,
        "kind": record.kind,
        "purpose": record.purpose,
        "subsystem": record.subsystem,
        "kernel": record.kernel,
        "tier": tier or record.tier,
        "review_status": record.review_status,
        "states": spec.states,
        "parameters": [
            {
                "symbol": p.symbol,
                "value": p.value.central,
                "units": p.value.units,
                "iiv_cv_percent": p.iiv_cv_percent,
                "tier": p.tier,
            }
            for p in record.parameters
        ],
        "kernel_values": kernel_values(record),
        "derivation_context": dc.__dict__ if dc else None,
        "transportability": {
            "validated_tumor_types": list(tp.validated_tumor_types),
            "validated_drug_classes": list(tp.validated_drug_classes),
            "validated_lines": list(tp.validated_lines),
            "out_of_context_action": tp.out_of_context_action,
        }
        if tp
        else None,
        "primary_citation": {
            "key": record.primary_citation.key if record.primary_citation else None,
            "doi": record.primary_citation.doi if record.primary_citation else None,
            "identifiers": identifier_uris(record),
        },
    }


def to_virtual_trial_json(record: Record, *, tier=None, dataset_version: str = _V) -> str:
    return json.dumps(
        virtual_trial_dict(record, tier=tier, dataset_version=dataset_version),
        indent=2,
        ensure_ascii=False,
    )
