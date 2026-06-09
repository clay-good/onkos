"""JSON-LD (linked-data) export.

Renders records as JSON-LD so the curation fields Onkos cares about — confidence
tier, clinical-use prohibition, derivation context, transportability, and the
``bqbiol:isDescribedBy`` DOI/PMID links — become real RDF triples a triple store
or reasoner can consume, not just JSON that happens to use ``onkos:`` keys.

The ``@context`` is the single one shipped in ``dataset/schema/context.jsonld``;
``tests/test_jsonld.py`` expands the output with rdflib and checks the expected
triples actually appear.
"""

from __future__ import annotations

import json
from functools import lru_cache

from .._const import CLINICAL_USE
from .._const import VERSION as _V
from .._data import dataset_dir
from ..load import Dataset
from ..models import Record
from .annotate import PREDICTION_PROHIBITED, identifier_uris, is_hypothesis_tier

_RECORD_IRI = "https://onkos.dev/record/"


@lru_cache(maxsize=1)
def load_context() -> dict:
    """The JSON-LD ``@context`` (from ``dataset/schema/context.jsonld``)."""
    path = dataset_dir() / "schema" / "context.jsonld"
    return json.loads(path.read_text())["@context"]


def record_node(record: Record, *, tier=None, dataset_version: str = _V) -> dict:
    """A single JSON-LD node (no ``@context``) for ``record``."""
    dc = record.derivation_context
    tp = record.transportability
    node = {
        "@id": _RECORD_IRI + record.id,
        "@type": "Model" if record.kind == "model" else "ContextBaseline",
        "recordId": record.id,
        "name": record.name,
        "kind": record.kind,
        "purpose": record.purpose,
        "subsystem": record.subsystem,
        "kernel": record.kernel,
        "confidenceTier": tier or record.tier,
        "reviewStatus": record.review_status,
        "clinicalUse": CLINICAL_USE,
        "datasetVersion": dataset_version,
    }
    if is_hypothesis_tier(record):
        node["predictionStatus"] = PREDICTION_PROHIBITED
    if dc:
        node["derivationContext"] = {
            k: v
            for k, v in {
                "drug": dc.drug,
                "drugClass": dc.drug_class,
                "tumorType": dc.tumor_type,
                "lineOfTherapy": dc.line_of_therapy,
            }.items()
            if v is not None
        }
    if tp:
        node["transportability"] = {
            "validatedTumorTypes": list(tp.validated_tumor_types),
            "validatedDrugClasses": list(tp.validated_drug_classes),
            "outOfContextAction": tp.out_of_context_action,
        }
    uris = identifier_uris(record)
    if uris:
        node["isDescribedBy"] = uris
    if record.primary_citation and record.primary_citation.doi:
        node["doi"] = record.primary_citation.doi
    node["hasParameter"] = [
        {
            "@id": f"{_RECORD_IRI}{record.id}#{p.symbol}",
            "symbol": p.symbol,
            "value": p.value.central,
            "units": p.value.units,
            "iivCvPercent": p.iiv_cv_percent,
            "confidenceTier": p.tier,
        }
        for p in record.parameters
    ]
    return node


def to_jsonld(record: Record, *, tier=None, dataset_version: str = _V) -> str:
    """A standalone JSON-LD document for one record."""
    doc = {"@context": load_context(), **record_node(record, tier=tier, dataset_version=dataset_version)}
    return json.dumps(doc, indent=2, ensure_ascii=False)


def dataset_jsonld(ds: Dataset) -> str:
    """The whole dataset as a single JSON-LD ``@graph``."""
    doc = {
        "@context": load_context(),
        "datasetVersion": ds.version,
        "@graph": [record_node(r, dataset_version=ds.version) for r in ds],
    }
    return json.dumps(doc, indent=2, ensure_ascii=False)
