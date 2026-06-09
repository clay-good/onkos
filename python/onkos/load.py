"""Load the dataset from disk into a typed :class:`Dataset`."""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from ._data import dataset_dir
from .models import (
    Citation,
    DerivationContext,
    Extraction,
    FailureMode,
    Parameter,
    PredictivePerformance,
    Record,
    Transportability,
    Value,
)

DATASET_VERSION = "0.1.0"


def _parameter(d: dict, citations: dict[str, Citation]) -> Parameter:
    v = d["value"]
    ex = d.get("extraction", {}) or {}
    return Parameter(
        symbol=d["symbol"],
        label=d["label"],
        value=Value(central=v["central"], units=v["units"], low=v.get("low"), high=v.get("high")),
        tier=d["tier"],
        iiv_cv_percent=d.get("iiv_cv_percent"),
        primary_citation=citations.get(d.get("primary_citation")),
        extraction=Extraction(
            review_status=ex.get("review_status", "unverified"),
            source_locator=ex.get("source_locator", ""),
            tier_rationale=ex.get("tier_rationale", ""),
        ),
    )


def _record(d: dict, citations: dict[str, Citation]) -> Record:
    dc = d.get("derivation_context")
    tp = d.get("transportability")
    return Record(
        id=d["id"],
        kind=d["kind"],
        purpose=d["purpose"],
        subsystem=d["subsystem"],
        tier=d["tier"],
        primary_citation=citations.get(d.get("primary_citation")),
        name=d.get("name", ""),
        description=d.get("description", ""),
        kernel=d.get("kernel"),
        structure=d.get("structure", {}),
        parameters=[_parameter(p, citations) for p in d.get("parameters", [])],
        derivation_context=DerivationContext(**dc) if dc else None,
        transportability=Transportability(
            validated_tumor_types=tuple(tp.get("validated_tumor_types", [])),
            validated_drug_classes=tuple(tp.get("validated_drug_classes", [])),
            validated_lines=tuple(tp.get("validated_lines", [])),
            out_of_context_action=tp.get("out_of_context_action", "warn"),
        )
        if tp
        else None,
        known_failure_modes=[FailureMode(**fm) for fm in d.get("known_failure_modes", [])],
        predictive_performance=[
            PredictivePerformance(**pp) for pp in d.get("predictive_performance", [])
        ],
        review_status=d.get("review_status", "unverified"),
        raw=d,
    )


class Dataset:
    """An in-memory, read-only view of the curated dataset."""

    def __init__(self, records: list[Record], citations: dict[str, Citation], version: str):
        self._records = {r.id: r for r in records}
        self.citations = citations
        self.version = version

    def __getitem__(self, record_id: str) -> Record:
        return self._records[record_id]

    def __contains__(self, record_id: str) -> bool:
        return record_id in self._records

    def __iter__(self) -> Iterator[Record]:
        return iter(self._records.values())

    def __len__(self) -> int:
        return len(self._records)

    @property
    def records(self) -> list[Record]:
        return list(self._records.values())

    def ids(self) -> list[str]:
        return sorted(self._records)

    def by_subsystem(self) -> Counter:
        return Counter(r.subsystem for r in self)

    def by_tier(self) -> Counter:
        return Counter(r.tier for r in self)

    def by_review_status(self) -> Counter:
        return Counter(r.review_status for r in self)


def load(path: str | None = None) -> Dataset:
    """Load all records and citations into a :class:`Dataset`."""
    base = Path(path) if path else dataset_dir()

    citations: dict[str, Citation] = {}
    cit_dir = base / "citations"
    if cit_dir.is_dir():
        for fp in sorted(cit_dir.glob("*.json")):
            d = json.loads(fp.read_text())
            citations[d["key"]] = Citation.from_dict(d)

    records: list[Record] = []
    for fp in sorted((base / "records").glob("*.json")):
        records.append(_record(json.loads(fp.read_text()), citations))

    return Dataset(records, citations, DATASET_VERSION)
