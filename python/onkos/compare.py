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
