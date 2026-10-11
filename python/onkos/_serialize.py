"""Strict JSON serialization shared by every result object's ``to_json``.

Analyses legitimately produce non-finite numbers — a growth rate that never occurred
(``nan``), a parameter a design cannot estimate (``inf`` RSE), a median not reached.
Python's ``json`` would write those as the bare tokens ``NaN`` / ``Infinity``, which are
not JSON and are rejected by strict parsers (``jq``, JavaScript ``JSON.parse``). They are
written as ``null`` instead: "no finite value", in every consumer's language.
"""

from __future__ import annotations

import json
import math


def _finite_or_none(obj):
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _finite_or_none(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_finite_or_none(v) for v in obj]
    return obj


def to_strict_json(obj, *, indent: int = 2) -> str:
    """``json.dumps`` that maps non-finite floats to ``null`` (always valid JSON)."""
    return json.dumps(_finite_or_none(obj), indent=indent, allow_nan=False)
