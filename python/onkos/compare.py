"""Virtual-trial divergence view — the headline feature.

Across *every eligible* TGI model for a chosen tumor type / line / drug-effect,
overlay the simulated tumor-size and population OS curves, grey out the models
whose ``transportability`` envelope the context violates (with the reason), and
quantify how much the survival prediction depends on the model choice.

This makes model-selection risk in go/no-go decisions measurable. It is forward,
population-level simulation only — never a per-patient prognosis.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .load import Dataset
from .simulate import Trajectory, median_survival, simulate
from .tiers import forces_tier_floor, transport_check


@dataclass
class Comparison:
    context: dict
    drug_effect: float
    t: np.ndarray
    included: list[Trajectory] = field(default_factory=list)
    excluded: list[tuple[str, str]] = field(default_factory=list)

    def _divergence(self, endpoint: str) -> float:
        curves = [tr.survival.get(endpoint) for tr in self.included]
        curves = [c for c in curves if c is not None]
        if len(curves) < 2:
            return 0.0
        stacked = np.vstack(curves)
        return float(np.max(stacked.max(axis=0) - stacked.min(axis=0)))

    def _median_range(self, endpoint: str) -> tuple[float, float] | None:
        meds = [
            median_survival(tr.t, tr.survival.get(endpoint))
            for tr in self.included
            if tr.survival.get(endpoint) is not None
        ]
        meds = [m for m in meds if m is not None]
        return (min(meds), max(meds)) if meds else None

    @property
    def os_divergence(self) -> float:
        """Max pointwise spread across the population OS curves (0..1)."""
        return self._divergence("OS")

    @property
    def pfs_divergence(self) -> float:
        """Max pointwise spread across the population PFS curves (0..1)."""
        return self._divergence("PFS")

    @property
    def median_os_range(self) -> tuple[float, float] | None:
        return self._median_range("OS")

    @property
    def median_pfs_range(self) -> tuple[float, float] | None:
        return self._median_range("PFS")

    def to_dict(self, *, include_curves: bool = False) -> dict:
        """A JSON-serializable summary of the virtual-trial result.

        Per-model summaries (tier, median OS/PFS, key TGI metrics, warnings),
        the excluded models with reasons, and the OS/PFS divergence. With
        ``include_curves`` the tumor/OS/PFS arrays are included as plain lists for
        a dashboard or external simulator to ingest (spec §7)."""
        from ._const import CLINICAL_USE

        def model(tr: Trajectory) -> dict:
            d = {
                "id": tr.record_id,
                "tier": tr.tier,
                "median_os_weeks": tr.median_os,
                "median_pfs_weeks": tr.median_pfs,
                "week8_relative_change": tr.metrics.get("week8_relative_change"),
                "depth_of_response": tr.metrics.get("depth_of_response"),
                "tumor_growth_rate_kg": tr.metrics.get("tumor_growth_rate_kg"),
                "warnings": list(tr.warnings),
            }
            if include_curves:
                d["t"] = self.t.tolist()
                d["tumor_size"] = tr.tumor_size.tolist()
                d["survival"] = {k: v.tolist() for k, v in tr.survival.items()}
            return d

        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "drug_effect": self.drug_effect,
            "n_included": len(self.included),
            "os_divergence": self.os_divergence,
            "pfs_divergence": self.pfs_divergence,
            "median_os_range": list(self.median_os_range) if self.median_os_range else None,
            "median_pfs_range": list(self.median_pfs_range) if self.median_pfs_range else None,
            "included": [model(tr) for tr in self.included],
            "excluded": [{"id": rid, "reason": reason} for rid, reason in self.excluded],
        }

    def to_json(self, *, include_curves: bool = False, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(include_curves=include_curves), indent=indent)


def compare(
    ds: Dataset,
    *,
    purpose: str = "tgi",
    context: dict | None = None,
    drug_effect: float = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> Comparison:
    context = context or {}
    tumor_type = context.get("tumor_type")
    line = context.get("line") or context.get("line_of_therapy")
    if t is None:
        t = np.linspace(0.0, 104.0, 209)
    t = np.asarray(t, dtype=float)

    cmp = Comparison(context=context, drug_effect=drug_effect, t=t)

    # The clinical divergence view excludes non-clinical subsystems: preclinical
    # xenograft models (validated in mice) and the hypothesis-tier immuno-oncology
    # QSP models (not for prediction).
    _excluded_subsystems = {"preclinical_translation", "immuno_oncology"}
    candidates = [
        r
        for r in ds
        if r.purpose == purpose
        and r.kind == "model"
        and r.kernel
        and r.subsystem not in _excluded_subsystems
    ]
    for r in sorted(candidates, key=lambda x: x.id):
        warns = transport_check(r, tumor_type=tumor_type, line=line)
        if warns and forces_tier_floor(r):
            cmp.excluded.append((r.id, "; ".join(warns)))
            continue
        traj = simulate(
            ds, r.id, context=context, drug_effect=drug_effect, exposure=exposure,
            exposure_response=exposure_response, t=t, survival_link=survival_link,
        )
        cmp.included.append(traj)

    return cmp
