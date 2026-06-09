"""Small filtering helpers over a :class:`~onkos.load.Dataset`."""

from __future__ import annotations

from collections.abc import Iterable

from .models import Record


def filter_records(
    records: Iterable[Record],
    *,
    subsystem: str | None = None,
    purpose: str | None = None,
    tier: str | None = None,
    tumor_type: str | None = None,
    review_status: str | None = None,
) -> list[Record]:
    out: list[Record] = []
    for r in records:
        if subsystem and r.subsystem != subsystem:
            continue
        if purpose and r.purpose != purpose:
            continue
        if tier and r.tier != tier:
            continue
        if review_status and r.review_status != review_status:
            continue
        if tumor_type:
            dc = r.derivation_context
            if not dc or dc.tumor_type != tumor_type:
                continue
        out.append(r)
    return out
