"""Exposure-response model choice — the dose-extrapolation model-selection axis.

The TGI -> survival chain begins with one upstream modeling choice that the rest of
Onkos has so far taken as given: the **exposure-response (ER) model** that maps a drug
exposure to the effect magnitude driving the kill term. Emax, sigmoid-Emax, and power
are all standard, all fit comparably to the studied dose, and all *diverge* when you
extrapolate to a dose that trial did not study — which is exactly what a dose-selection
decision asks them to do.

This module makes that the project's core transportability thesis applied one layer
upstream. It **anchors** each ER shape to agree at a reference exposure ``(c_ref, e_ref)``
— the studied dose — so the curves are indistinguishable there, then reads off how far
their predicted effect (and the resulting population OS) diverge at other doses. The
headline: the divergence is ~0 at the studied dose and grows as you extrapolate,
especially *downward* (de-escalation), where the effect sits on the steep part of the
survival relationship. A dose-response model fit at one dose carries an unquantified
model-selection risk the moment it is used to pick another dose.

Population / trial level only. The reference point and the shapes are declared, the
analysis re-anchors the curated ER shapes (it does not refit them to data). NOT a dose
recommendation, NOT a per-patient prediction. A dose-response analysis never moves a
tier; it inherits the propagated tier of the chain it drives.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .export.registry import get_kernel, kernel_values
from .load import Dataset
from .simulate import simulate

__all__ = [
    "ER_SHAPE_RECORDS",
    "calibrated_er",
    "compare_er_extrapolation",
    "ExtrapolationComparison",
]

# The curated ER shape families used as the default model set: each is re-anchored to a
# common reference point, so the comparison isolates the SHAPE (saturating vs unbounded
# vs switch-like), not the curated scale.
ER_SHAPE_RECORDS = (
    "exposure_response.emax_generic",
    "exposure_response.power_generic",
    "exposure_response.sigmoid_emax_generic",
)


def calibrated_er(ds: Dataset, er_id: str, *, c_ref: float, e_ref: float) -> Callable[[float], float]:
    """Return the forward exposure->effect map of ``er_id``'s SHAPE, re-anchored so it
    passes through ``(c_ref, e_ref)``.

    The shape parameters (EC50, gamma, theta) are taken from the curated record; the
    single scale parameter (Emax or slope) is solved so the curve hits ``e_ref`` at
    ``c_ref``. So all shapes agree exactly at the reference dose and differ only in how
    they extrapolate away from it."""
    if c_ref <= 0 or e_ref <= 0:
        raise ValueError("c_ref and e_ref must be positive")
    rec = ds[er_id]
    name = get_kernel(rec).name
    v = kernel_values(rec)
    if name == "er_emax":
        ec50 = float(v["EC50"])
        emax = e_ref * (ec50 + c_ref) / c_ref

        def f(c):
            return emax * c / (ec50 + c)

        return f
    if name == "er_sigmoid_emax":
        ec50, g = float(v["EC50"]), float(v["gamma"])
        emax = e_ref * (ec50**g + c_ref**g) / c_ref**g

        def f(c):
            return emax * c**g / (ec50**g + c**g)

        return f
    if name == "er_power":
        theta = float(v["theta"])
        slope = e_ref / c_ref**theta

        def f(c):
            return slope * c**theta

        return f
    raise ValueError(f"no calibratable ER shape for kernel {get_kernel(rec).name!r} (record {er_id})")


@dataclass
class ExtrapolationComparison:
    """Effect and OS across ER shapes, over a dose grid, all anchored at ``(c_ref, e_ref)``."""

    record_id: str
    context: dict
    c_ref: float
    e_ref: float
    er_ids: list
    doses: np.ndarray
    rows: list = field(default_factory=list)  # per dose: {dose, effects, median_os, effect_divergence, os_divergence}
    tier: str = "C"
    warnings: list = field(default_factory=list)
    clinical_use: str = CLINICAL_USE

    def _row(self, dose: float) -> dict:
        return min(self.rows, key=lambda r: abs(r["dose"] - dose))

    def os_divergence_at(self, dose: float) -> float:
        return self._row(dose)["os_divergence"]

    def effect_divergence_at(self, dose: float) -> float:
        return self._row(dose)["effect_divergence"]

    @property
    def reference_os_divergence(self) -> float:
        """OS spread at the reference dose — ~0 by construction (the curves are anchored
        there), the control that proves the divergence elsewhere is the extrapolation."""
        return self._row(self.c_ref)["os_divergence"]

    @property
    def max_os_divergence(self) -> float:
        return max((r["os_divergence"] for r in self.rows), default=0.0)

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "record_id": self.record_id,
            "context": self.context,
            "c_ref": self.c_ref,
            "e_ref": self.e_ref,
            "er_ids": list(self.er_ids),
            "tier": self.tier,
            "reference_os_divergence": self.reference_os_divergence,
            "max_os_divergence": self.max_os_divergence,
            "rows": [
                {
                    "dose": r["dose"],
                    "effects": r["effects"],
                    "median_os_weeks": r["median_os"],
                    "effect_divergence": r["effect_divergence"],
                    "os_divergence": r["os_divergence"],
                }
                for r in self.rows
            ],
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def compare_er_extrapolation(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    er_ids: tuple | list | None = None,
    c_ref: float = 150.0,
    e_ref: float = 1.0,
    doses: np.ndarray | None = None,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> ExtrapolationComparison:
    """Anchor every ER shape in ``er_ids`` at ``(c_ref, e_ref)`` and, over ``doses``,
    report the effect each predicts and the population OS it drives through ``record_id``'s
    TGI -> survival chain — quantifying how much the ER-model choice matters at each dose.

    The spread is ~0 at ``c_ref`` (the anchor) and grows on extrapolation; the downward
    extrapolation (de-escalation) typically diverges most in OS."""
    context = context or {}
    er_ids = list(er_ids) if er_ids is not None else list(ER_SHAPE_RECORDS)
    if doses is None:
        doses = c_ref * np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    doses = np.asarray(doses, dtype=float)
    if t is None:
        t = np.linspace(0.0, 156.0, 313)

    curves = {er: calibrated_er(ds, er, c_ref=c_ref, e_ref=e_ref) for er in er_ids}
    cmp = ExtrapolationComparison(
        record_id=record_id, context=context, c_ref=float(c_ref), e_ref=float(e_ref),
        er_ids=er_ids, doses=doses,
    )
    tier = "C"
    for dose in doses:
        effects, medians = {}, {}
        for er, f in curves.items():
            e = max(float(f(dose)), 0.0)
            tr = simulate(ds, record_id, context=context, drug_effect=e, t=t,
                          survival_link=survival_link)
            tier = tr.tier
            effects[er] = e
            medians[er] = tr.median_os
        e_vals = list(effects.values())
        m_vals = [m for m in medians.values() if m is not None]
        cmp.rows.append({
            "dose": float(dose),
            "effects": effects,
            "median_os": medians,
            "effect_divergence": float(max(e_vals) - min(e_vals)) if e_vals else 0.0,
            "os_divergence": float(max(m_vals) - min(m_vals)) if len(m_vals) >= 2 else 0.0,
        })
    cmp.tier = tier
    return cmp
