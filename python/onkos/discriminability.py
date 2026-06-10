"""Model discriminability — can a trial even tell the competing models apart?

The model-selection arc (v0.21-v0.37) quantified how much the survival forecast depends on
modeling choices: which TGI model, which resistance mechanism and origin, which bridge
metric, which survival structure. This module asks the question that closes the loop: given
two models' population OS curves, **what trial would it take to distinguish them?** — and
conversely, are the divergences we have been quantifying even *detectable* in a realistic
trial?

The answer reframes the project's load-bearing idea. A model-selection divergence that needs
tens of thousands of events to detect is **practically unidentifiable from the trial**: the
model choice cannot be resolved by the data, only by assumption. That is exactly what makes
the silent-transport risk silent — and here it is quantified, in events.

The discriminability of two survival curves is a logrank power calculation. The required
number of events to distinguish them at power ``1-beta`` and two-sided level ``alpha`` is the
Schoenfeld formula ``d = 4 (z_{1-alpha/2} + z_{1-beta})^2 / (ln HR)^2`` for 1:1 allocation,
where ``HR`` is the follow-up-horizon hazard ratio between the curves (the ratio of cumulative
hazards at the trial horizon — exact under proportional hazards, an average otherwise).

Design / trial level only — this is a power calculation over published model structures, the
same family as :mod:`onkos.identify` (parameter identifiability) and :mod:`onkos.design`
(optimal schedule), one level up at the **model**-identifiability question. NOT a real trial
design, NOT a recommendation, NOT an individual quantity. A discriminability analysis never
moves a tier.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .compare import compare
from .load import Dataset
from .simulate import simulate

__all__ = [
    "required_events",
    "horizon_hazard_ratio",
    "discriminating_events",
    "model_discriminability",
    "ModelDiscriminability",
]

# Illustrative phase-3-scale bounds for classifying a comparison's feasibility (declared, not
# a clinical threshold): a large oncology OS trial reads out a few hundred events; a few
# thousand is at the edge of feasibility; tens of thousands is practically impossible.
_FEASIBLE_EVENTS = 500.0
_INFEASIBLE_EVENTS = 3000.0


def _ppf(p: float) -> float:
    """Standard-normal quantile (inverse CDF). scipy is already a project dependency."""
    from scipy.special import ndtri

    return float(ndtri(p))


def required_events(hazard_ratio: float, *, power: float = 0.8, alpha: float = 0.05) -> float:
    """Events to distinguish two survival curves by a logrank test at ``power`` and two-sided
    ``alpha`` (Schoenfeld, 1:1 allocation): ``d = 4 (z_{1-a/2}+z_{1-b})^2 / (ln HR)^2``.

    ``inf`` when ``HR = 1`` (identical curves — never distinguishable). Symmetric in
    ``HR`` and ``1/HR`` (the events depend only on ``|ln HR|``)."""
    hr = float(hazard_ratio)
    if not np.isfinite(hr) or hr <= 0.0 or hr == 1.0:
        return float("inf")
    z = _ppf(1.0 - alpha / 2.0) + _ppf(power)
    return 4.0 * z * z / (np.log(hr) ** 2)


def horizon_hazard_ratio(curve_a, curve_b, *, floor: float = 1e-6) -> float:
    """The follow-up-horizon hazard ratio between two survival curves: the ratio of their
    cumulative hazards ``H = -ln S`` at the trial horizon (the last point). Exact under
    proportional hazards; a horizon-average otherwise. ``inf`` if the reference curve has no
    cumulative hazard (flat at 1)."""
    sa = float(np.clip(np.asarray(curve_a, dtype=float)[-1], floor, 1.0))
    sb = float(np.clip(np.asarray(curve_b, dtype=float)[-1], floor, 1.0))
    ha, hb = -np.log(sa), -np.log(sb)
    if hb <= 0.0:
        return float("inf")
    return ha / hb


def discriminating_events(
    ds: Dataset,
    record_a: str,
    record_b: str,
    *,
    context: dict,
    survival_link: str | None = None,
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
    power: float = 0.8,
    alpha: float = 0.05,
) -> dict:
    """The horizon hazard ratio and required events to distinguish two models' OS curves in
    ``context`` under ``survival_link`` (the default link if None)."""
    if t is None:
        t = np.linspace(0.0, 312.0, 625)
    t = np.asarray(t, dtype=float)
    ca = simulate(ds, record_a, context=context, drug_effect=drug_effect, t=t,
                  survival_link=survival_link).os_curve
    cb = simulate(ds, record_b, context=context, drug_effect=drug_effect, t=t,
                  survival_link=survival_link).os_curve
    if ca is None or cb is None:
        return {"record_a": record_a, "record_b": record_b, "hazard_ratio": float("nan"),
                "required_events": float("inf")}
    hr = horizon_hazard_ratio(ca, cb)
    return {
        "record_a": record_a,
        "record_b": record_b,
        "hazard_ratio": hr,
        "required_events": required_events(hr, power=power, alpha=alpha),
    }


@dataclass
class ModelDiscriminability:
    """Pairwise discriminability of a context's TGI models under one survival link."""

    context: dict
    survival_link: str | None
    power: float
    alpha: float
    models: list
    pairs: list = field(default_factory=list)  # {record_a, record_b, hazard_ratio, required_events}
    tier: str = "C"
    clinical_use: str = CLINICAL_USE

    @property
    def indistinguishable_pairs(self) -> list:
        """Pairs needing more than the infeasible-events bound — the model choice cannot be
        resolved by a realistic trial (practically unidentifiable from OS)."""
        return [p for p in self.pairs if p["required_events"] > _INFEASIBLE_EVENTS]

    @property
    def feasible_pairs(self) -> list:
        return [p for p in self.pairs if p["required_events"] < _FEASIBLE_EVENTS]

    @property
    def n_indistinguishable(self) -> int:
        return len(self.indistinguishable_pairs)

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "survival_link": self.survival_link,
            "power": self.power,
            "alpha": self.alpha,
            "tier": self.tier,
            "n_models": len(self.models),
            "n_indistinguishable": self.n_indistinguishable,
            "feasible_events_bound": _FEASIBLE_EVENTS,
            "infeasible_events_bound": _INFEASIBLE_EVENTS,
            "pairs": [
                {**p, "required_events": (None if not np.isfinite(p["required_events"])
                                          else round(p["required_events"], 1))}
                for p in self.pairs
            ],
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def model_discriminability(
    ds: Dataset,
    *,
    context: dict,
    survival_link: str | None = None,
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
    power: float = 0.8,
    alpha: float = 0.05,
) -> ModelDiscriminability:
    """For every pair of a context's eligible TGI models, the required events to distinguish
    their OS curves under ``survival_link``. Flags pairs that need an infeasible trial —
    model choices a realistic study cannot resolve, only assume.

    Under the default week-8 link the resistance-model pairs (which diverge only in the
    regrowth tail) are typically practically indistinguishable, while pairs that differ in
    early shrinkage are easily distinguished — the v0.24 / v0.32 finding made quantitative."""
    if t is None:
        t = np.linspace(0.0, 312.0, 625)
    t = np.asarray(t, dtype=float)
    cmp = compare(ds, purpose="tgi", context=context, drug_effect=drug_effect, t=t)
    models = [tr.record_id for tr in cmp.included]
    tier = cmp.included[0].tier if cmp.included else "C"
    curves = {
        tr.record_id: simulate(ds, tr.record_id, context=context, drug_effect=drug_effect, t=t,
                               survival_link=survival_link).os_curve
        for tr in cmp.included
    }
    md = ModelDiscriminability(
        context=context, survival_link=survival_link, power=power, alpha=alpha,
        models=models, tier=tier,
    )
    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            a, b = models[i], models[j]
            ca, cb = curves[a], curves[b]
            if ca is None or cb is None:
                hr, d = float("nan"), float("inf")
            else:
                hr = horizon_hazard_ratio(ca, cb)
                d = required_events(hr, power=power, alpha=alpha)
            md.pairs.append({"record_a": a, "record_b": b, "hazard_ratio": hr, "required_events": d})
    return md
