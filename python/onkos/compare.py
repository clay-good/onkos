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

    @property
    def os_divergence(self) -> float:
        """Max pointwise spread across the population OS curves (0..1)."""
        curves = [tr.os_curve for tr in self.included if tr.os_curve is not None]
        if len(curves) < 2:
            return 0.0
        stacked = np.vstack(curves)
        return float(np.max(stacked.max(axis=0) - stacked.min(axis=0)))

    @property
    def median_os_range(self) -> tuple[float, float] | None:
        meds = [
            median_survival(tr.t, tr.os_curve)
            for tr in self.included
            if tr.os_curve is not None
        ]
        meds = [m for m in meds if m is not None]
        if not meds:
            return None
        return (min(meds), max(meds))


def compare(
    ds: Dataset,
    *,
    purpose: str = "tgi",
    context: dict | None = None,
    drug_effect: float = 1.0,
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

    candidates = [
        r for r in ds if r.purpose == purpose and r.kind == "model" and r.kernel
    ]
    for r in sorted(candidates, key=lambda x: x.id):
        warns = transport_check(r, tumor_type=tumor_type, line=line)
        if warns and forces_tier_floor(r):
            cmp.excluded.append((r.id, "; ".join(warns)))
            continue
        traj = simulate(
            ds, r.id, context=context, drug_effect=drug_effect, t=t, survival_link=survival_link
        )
        cmp.included.append(traj)

    return cmp
