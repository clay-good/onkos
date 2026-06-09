"""Typed object model for Onkos records.

Records are structured objects, not scalars. A :class:`Record` is either a
``model`` (structure + parameters + derivation context) or a
``context_baseline`` (parameters with tumor-type/line provenance). Both share
one schema. Parameters are reachable by symbol via ``record["lambda"]``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

TIER_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3}


@dataclass(frozen=True)
class Citation:
    key: str
    title: str = ""
    authors: tuple = ()
    journal: str = ""
    year: int | None = None
    volume: str = ""
    issue: str = ""
    pages: str = ""
    doi: str | None = None
    pmid: str | None = None

    @classmethod
    def from_dict(cls, d: dict) -> Citation:
        return cls(
            key=d["key"],
            title=d.get("title", ""),
            authors=tuple(d.get("authors", [])),
            journal=d.get("journal", ""),
            year=d.get("year"),
            volume=str(d.get("volume", "")),
            issue=str(d.get("issue", "")),
            pages=str(d.get("pages", "")),
            doi=d.get("doi"),
            pmid=d.get("pmid"),
        )

    def bibtex(self) -> str:
        authors = " and ".join(self.authors)
        fields = [
            ("title", self.title),
            ("author", authors),
            ("journal", self.journal),
            ("year", self.year),
            ("volume", self.volume),
            ("number", self.issue),
            ("pages", self.pages),
            ("doi", self.doi),
        ]
        body = ",\n".join(
            f"  {k} = {{{v}}}" for k, v in fields if v not in (None, "", ())
        )
        return f"@article{{{self.key},\n{body}\n}}"


@dataclass(frozen=True)
class Value:
    central: float
    units: str
    low: float | None = None
    high: float | None = None


@dataclass(frozen=True)
class Extraction:
    review_status: str = "unverified"
    source_locator: str = ""
    tier_rationale: str = ""


@dataclass(frozen=True)
class Parameter:
    symbol: str
    label: str
    value: Value
    tier: str
    iiv_cv_percent: float | None = None
    primary_citation: Citation | None = None
    extraction: Extraction = field(default_factory=Extraction)

    @property
    def central(self) -> float:
        return self.value.central

    @property
    def units(self) -> str:
        return self.value.units


@dataclass(frozen=True)
class DerivationContext:
    drug: str | None = None
    drug_class: str | None = None
    tumor_type: str | None = None
    line_of_therapy: str | None = None
    trial: str | None = None
    n_patients: int | None = None
    measurement: str | None = None


@dataclass(frozen=True)
class Transportability:
    validated_tumor_types: tuple = ()
    validated_drug_classes: tuple = ()
    validated_lines: tuple = ()
    out_of_context_action: str = "warn"


@dataclass(frozen=True)
class FailureMode:
    condition: str
    behavior: str
    action: str
    citation: str | None = None


@dataclass(frozen=True)
class PredictivePerformance:
    metric: str
    value: float
    population: str | None = None
    citation: str | None = None


@dataclass
class Record:
    id: str
    kind: str
    purpose: str
    subsystem: str
    tier: str
    primary_citation: Citation | None
    name: str = ""
    description: str = ""
    kernel: str | None = None
    structure: dict = field(default_factory=dict)
    parameters: list[Parameter] = field(default_factory=list)
    derivation_context: DerivationContext | None = None
    transportability: Transportability | None = None
    known_failure_modes: list[FailureMode] = field(default_factory=list)
    predictive_performance: list[PredictivePerformance] = field(default_factory=list)
    review_status: str = "unverified"
    raw: dict = field(default_factory=dict)

    @property
    def _by_symbol(self) -> dict[str, Parameter]:
        return {p.symbol: p for p in self.parameters}

    def __getitem__(self, symbol: str) -> Parameter:
        return self._by_symbol[symbol]

    def __contains__(self, symbol: str) -> bool:
        return symbol in self._by_symbol

    def get(self, symbol: str, default=None):
        return self._by_symbol.get(symbol, default)

    def param_values(self) -> dict[str, float]:
        """Central values keyed by symbol, for kernel binding."""
        return {p.symbol: p.value.central for p in self.parameters}
