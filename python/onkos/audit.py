"""Evidence-based tier audit — make "tier assignment is partly numeric" real.

The spec (§5, §9) says a record's confidence tier is partly a *numeric* judgment
driven by external validation: did the TGI-metric → survival link hold in an
independent trial (external C-index), and how poorly identified are its kill /
resistance terms (IIV CV). Tiers are otherwise hand-set; this module derives the
*best tier the recorded evidence supports* (the "ceiling") and flags any record
whose assigned tier is **better** than its evidence — tier inflation, the
dangerous direction. The check runs in ``onkos validate`` so it cannot regress.

Only records whose tier reflects TGI→survival evidence are audited: clinical TGI
and survival-link models. Growth laws (structural forms), exposure-response, and
context baselines are tiered on other criteria and are left out. Preclinical and
hypothesis-tier immuno-oncology are governed by their own rules.
"""

from __future__ import annotations

from dataclasses import dataclass

from .load import Dataset
from .models import TIER_ORDER, Record

_EVIDENCE_TIERED_PURPOSES = {"tgi", "survival_link"}
_NON_CLINICAL_SUBSYSTEMS = {"preclinical_translation", "immuno_oncology"}
# A parameter known only to this CV or worse is "poorly identified" (spec §5's
# tier-C characteristic), which caps the achievable tier at C regardless of any
# external check — the resistance/kill term simply isn't pinned down.
_POORLY_IDENTIFIED_CV = 70.0


def is_evidence_tiered(record: Record) -> bool:
    return (
        record.purpose in _EVIDENCE_TIERED_PURPOSES
        and record.subsystem not in _NON_CLINICAL_SUBSYSTEMS
    )


def _has_external_validation(record: Record) -> bool:
    return any("external" in pp.metric.lower() for pp in record.predictive_performance)


def _max_iiv(record: Record) -> float:
    cvs = [p.iiv_cv_percent for p in record.parameters if p.iiv_cv_percent]
    return max(cvs) if cvs else 0.0


def _breadth(record: Record) -> int:
    tp = record.transportability
    return len(tp.validated_tumor_types) if tp else 0


def evidence_ceiling(record: Record) -> str:
    """The best tier the record's recorded evidence supports.

    A: external validation + broad context (>=2 validated tumor types).
    B: external validation (a single external check), parameters reasonably identified.
    C: no external validation, OR a poorly-identified (IIV CV >= 70%) kill/resistance
       term — spec §5's tier-C characteristic, which caps the ceiling at C."""
    has_external = _has_external_validation(record)
    breadth = _breadth(record)
    if has_external and breadth >= 2:
        ceiling = "A"
    elif has_external:
        ceiling = "B"
    else:
        ceiling = "C"
    if _max_iiv(record) >= _POORLY_IDENTIFIED_CV and TIER_ORDER[ceiling] < TIER_ORDER["C"]:
        ceiling = "C"
    return ceiling


@dataclass
class TierFinding:
    record_id: str
    assigned: str
    ceiling: str
    has_external: bool
    max_iiv: float
    breadth: int
    status: str  # "inflated" | "conservative" | "ok"


def audit_tiers(ds: Dataset) -> list[TierFinding]:
    findings = []
    for r in ds:
        if not is_evidence_tiered(r):
            continue
        ceiling = evidence_ceiling(r)
        if TIER_ORDER[r.tier] < TIER_ORDER[ceiling]:
            status = "inflated"  # assigned a better tier than the evidence supports
        elif TIER_ORDER[r.tier] > TIER_ORDER[ceiling]:
            status = "conservative"  # could be upgraded if the evidence is trusted
        else:
            status = "ok"
        findings.append(
            TierFinding(
                record_id=r.id,
                assigned=r.tier,
                ceiling=ceiling,
                has_external=_has_external_validation(r),
                max_iiv=_max_iiv(r),
                breadth=_breadth(r),
                status=status,
            )
        )
    return findings


def inflated_records(ds: Dataset) -> list[TierFinding]:
    return [f for f in audit_tiers(ds) if f.status == "inflated"]
