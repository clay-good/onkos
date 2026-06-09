"""virtual-trial JSON — parameters + tumor context + survival link for the
simulator/dashboard, carrying the mandatory NOT FOR CLINICAL USE flag."""

from __future__ import annotations

import json

from .._const import CLINICAL_USE
from ..models import Record
from .annotate import identifier_uris
from .registry import get_kernel, kernel_values


def virtual_trial_dict(record: Record, *, tier=None, dataset_version: str = "0.1.0") -> dict:
    spec = get_kernel(record)
    dc = record.derivation_context
    tp = record.transportability
    return {
        "onkos:clinicalUse": CLINICAL_USE,
        "NOT_FOR_CLINICAL_USE": True,
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


def to_virtual_trial_json(record: Record, *, tier=None, dataset_version: str = "0.1.0") -> str:
    return json.dumps(
        virtual_trial_dict(record, tier=tier, dataset_version=dataset_version),
        indent=2,
        ensure_ascii=False,
    )
