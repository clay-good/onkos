"""Early-surrogate readout timing — *when* you read the surrogate is a model-selection axis.

The metric-choice work (v0.25 week-8 vs k_g, v0.33 integrated burden) asked *which*
on-treatment quantity predicts survival. This module asks the orthogonal question the
ctDNA era has made urgent: *when* do you read it? The field is pushing the surrogate
readout ever earlier — circulating-tumor-DNA (ctDNA) "molecular response" at week 2-4,
well before a RECIST tumor-size change is reliable at week 8 — on the premise that an
earlier signal is an earlier go/no-go.

Onkos models ctDNA molecular response as proportional to tumor burden (the standard
first-order shedding assumption), so the modeled distinction between a ctDNA readout and
a RECIST-size readout is purely the **landmark time** (genomic / mutational ctDNA content
is out of scope, spec §2). That isolates the timing question cleanly: ``landmark_response``
is the relative tumor-burden change at an arbitrary landmark week, and
``surrogate_timing_fidelity`` measures how well the early-landmark ranking of a context's
models agrees with a tail-aware *durable-benefit* reference (the k_g survival link).

The finding: **earliness trades against fidelity.** An earlier landmark sits closer to the
nadir, before any resistant regrowth, so it over-rewards deep-but-doomed early responders;
the discordance with the durable-benefit ranking falls monotonically as the landmark moves
later. The ctDNA-driven push to week-2-4 readouts maximizes that surrogate-timing bias.

Population / trial level only. NOT an individual molecular-response prediction, NOT a
go/no-go recommendation. A timing analysis never moves a tier; it inherits the propagated
tier of the trajectories it ranks.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .compare import compare
from .load import Dataset
from .simulate import simulate

__all__ = [
    "landmark_response",
    "discordant_pairs",
    "surrogate_timing_fidelity",
    "SurrogateTiming",
]


def landmark_response(t: np.ndarray, v: np.ndarray, week: float) -> float:
    """Relative tumor-burden change at ``week`` (the early-surrogate readout):
    ``(v(week) - v0) / v0``, ``v0 = v[0]``.

    Generalizes the fixed week-8 surrogate to an arbitrary landmark; under the
    proportional-shedding assumption this is also the modeled ctDNA molecular response at
    ``week``. More negative = deeper early response."""
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    y0 = float(v[0])
    if y0 == 0.0:
        return float("nan")
    return float((np.interp(week, t, v) - y0) / y0)


def discordant_pairs(rank_a: list, rank_b: list) -> int:
    """Number of item pairs the two rankings order oppositely (a Kendall-style count).
    Both rankings must contain the same items."""
    pa = {x: i for i, x in enumerate(rank_a)}
    pb = {x: i for i, x in enumerate(rank_b)}
    items = list(rank_a)
    n = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i], items[j]
            if (pa[a] - pa[b]) * (pb[a] - pb[b]) < 0:
                n += 1
    return n


@dataclass
class SurrogateTiming:
    """How an early-surrogate model ranking's fidelity to a durable-benefit reference
    varies with the readout landmark week."""

    context: dict
    reference_link: str
    landmark_weeks: list
    reference_ranking: list  # model ids, best (longest durable OS) first
    rows: list = field(default_factory=list)  # per week: {week, ranking, discordant_pairs}
    total_pairs: int = 0
    tier: str = "C"
    clinical_use: str = CLINICAL_USE

    def _row(self, week: float) -> dict:
        return min(self.rows, key=lambda r: abs(r["week"] - week))

    def discordance_at(self, week: float) -> int:
        return self._row(week)["discordant_pairs"]

    @property
    def earliest_discordance(self) -> int:
        return self.rows[0]["discordant_pairs"] if self.rows else 0

    @property
    def latest_discordance(self) -> int:
        return self.rows[-1]["discordant_pairs"] if self.rows else 0

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "reference_link": self.reference_link,
            "reference_ranking": list(self.reference_ranking),
            "total_pairs": self.total_pairs,
            "tier": self.tier,
            "rows": [
                {"week": r["week"], "discordant_pairs": r["discordant_pairs"], "ranking": list(r["ranking"])}
                for r in self.rows
            ],
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def _median_os(ds, rid, context, link, t):
    return simulate(ds, rid, context=context, drug_effect=1.0, t=t, survival_link=link).median_os


def surrogate_timing_fidelity(
    ds: Dataset,
    *,
    context: dict,
    landmark_weeks=(2.0, 4.0, 8.0, 12.0, 16.0, 24.0, 36.0, 52.0),
    reference_link: str | None = None,
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
) -> SurrogateTiming:
    """For each landmark week, rank a context's eligible TGI models by their early-surrogate
    response (most shrinkage = best) and count how many model pairs that ranking orders
    oppositely to the tail-aware **durable-benefit** reference ranking (median OS under the
    k_g survival link, which rewards slow regrowth).

    The earlier the landmark, the more the surrogate over-rewards deep-but-doomed early
    responders — so the discordance falls as the landmark moves later."""
    if t is None:
        t = np.linspace(0.0, 260.0, 521)
    t = np.asarray(t, dtype=float)
    tumor_type = context.get("tumor_type")

    if reference_link is None:
        # The tail-aware durable-benefit reference: the context's k_g (growth-rate) OS link.
        reference_link = f"survival_link.{str(tumor_type).lower()}_os_growth_rate"

    cmp = compare(ds, purpose="tgi", context=context, drug_effect=drug_effect, t=t)
    models = [tr.record_id for tr in cmp.included]
    trajectories = {tr.record_id: (tr.t, tr.tumor_size) for tr in cmp.included}
    tier = cmp.included[0].tier if cmp.included else "C"

    # Durable-benefit reference ranking: longer k_g-link OS = better.
    ref_os = {rid: (_median_os(ds, rid, context, reference_link, t) or -1.0) for rid in models}
    reference_ranking = sorted(models, key=lambda r: -ref_os[r])

    st = SurrogateTiming(
        context=context, reference_link=reference_link, landmark_weeks=list(landmark_weeks),
        reference_ranking=reference_ranking, total_pairs=len(models) * (len(models) - 1) // 2,
        tier=tier,
    )
    for w in landmark_weeks:
        resp = {rid: landmark_response(*trajectories[rid], w) for rid in models}
        # most shrinkage (most negative) = the surrogate's "best" model
        ranking = sorted(models, key=lambda r: resp[r])
        st.rows.append({
            "week": float(w),
            "ranking": ranking,
            "discordant_pairs": discordant_pairs(ranking, reference_ranking),
        })
    return st
