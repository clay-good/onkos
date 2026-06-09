"""Confidence-tier propagation and transportability enforcement.

The load-bearing rule, shared with Nidus and Hypnos: *worst input wins*. A
composed prediction inherits the worst tier among its components. Applying a
record outside its validated ``transportability`` envelope, or tripping a
``known_failure_mode``, forces a tier floor of ``D`` with an attached warning.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from .models import TIER_ORDER, Record


@dataclass
class TierResult:
    tier: str
    warnings: list[str]


def worst_tier(tiers: Iterable[str]) -> str:
    tiers = [t for t in tiers if t in TIER_ORDER]
    if not tiers:
        return "D"
    return max(tiers, key=lambda t: TIER_ORDER[t])


def transport_check(record: Record, *, tumor_type: str | None, line: str | None) -> list[str]:
    """Return warnings if applying ``record`` to this context leaves its
    validated envelope. Empty list means in-context."""
    warnings: list[str] = []
    tp = record.transportability
    if tp is None:
        return warnings

    if tumor_type and tp.validated_tumor_types and tumor_type not in tp.validated_tumor_types:
        warnings.append(
            f"{record.id}: tumor_type '{tumor_type}' is outside validated "
            f"{list(tp.validated_tumor_types)} -> {tp.out_of_context_action}"
        )
    if line and tp.validated_lines and line not in tp.validated_lines:
        warnings.append(
            f"{record.id}: line '{line}' is outside validated "
            f"{list(tp.validated_lines)} -> {tp.out_of_context_action}"
        )
    return warnings


def forces_tier_floor(record: Record) -> bool:
    tp = record.transportability
    return tp is not None and tp.out_of_context_action == "tier_down_to_D and warn"


def propagate(
    records: Iterable[Record], *, tumor_type: str | None = None, line: str | None = None
) -> TierResult:
    """Compose the tier for a simulation built from ``records`` in a context.

    Starts from the worst component tier, then applies transportability and
    failure-mode floors.
    """
    records = list(records)
    warnings: list[str] = []
    tier = worst_tier(r.tier for r in records)

    floored = False
    for r in records:
        w = transport_check(r, tumor_type=tumor_type, line=line)
        if w:
            warnings.extend(w)
            if forces_tier_floor(r):
                floored = True
        for fm in r.known_failure_modes:
            # Failure modes are condition-described; we surface them as advisories
            # so a user knows the regime in which this record becomes unstable.
            warnings.append(f"{r.id}: known failure mode — {fm.condition} ({fm.action})")

    if floored:
        tier = "D"
    return TierResult(tier=tier, warnings=warnings)
