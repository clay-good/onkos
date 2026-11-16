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
from .tiers import worst_tier

__all__ = [
    "landmark_response",
    "discordant_pairs",
    "readout_se",
    "ranking_flip_probability",
    "surrogate_timing_fidelity",
    "SurrogateTiming",
]


def landmark_response(t: np.ndarray, v: np.ndarray, week: float) -> float:
    """Relative tumor-burden change at ``week`` (the early-surrogate readout):
    ``(v(week) - v0) / v0``, ``v0 = v[0]``.

    Generalizes the fixed week-8 surrogate to an arbitrary landmark; under the
    proportional-shedding assumption this is also the modeled ctDNA molecular response at
    ``week``. More negative = deeper early response. ``nan`` if ``week`` lies off the grid."""
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    y0 = float(v[0])
    # A landmark off the grid is unobserved: np.interp would clamp to the nearest end and
    # report, say, the week-4 response as the week-8 readout.
    if y0 == 0.0 or not (t[0] <= week <= t[-1]):
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


_LOG_FLOOR = 1e-3  # detection-limit floor on relative size before taking logs (as in metrics.py)


def readout_se(noise_cv: float, n_patients: int) -> float:
    """Standard error of the difference between two arms' mean log-change readouts.

    Each patient's readout is ``log(v_week / v_0)``, both measurements carrying independent
    multiplicative noise with coefficient of variation ``noise_cv`` (log-SD
    ``s = sqrt(log(1 + cv^2))``), so one patient's log-change has variance ``2 s^2``, an
    arm mean ``2 s^2 / n`` and the difference of two arms ``4 s^2 / n``: ``SE = 2 s / sqrt(n)``."""
    if not np.isfinite(noise_cv) or noise_cv < 0:
        raise ValueError(f"noise_cv={noise_cv} must be a finite number >= 0")
    if int(n_patients) != n_patients or n_patients < 1:
        raise ValueError(f"n_patients={n_patients} must be a positive integer")
    return 2.0 * float(np.sqrt(np.log1p(noise_cv * noise_cv))) / float(np.sqrt(n_patients))


def ranking_flip_probability(delta: float, se: float) -> float:
    """Probability that a noisy readout orders two models opposite to their true order:
    ``Phi(-|delta| / se)`` for a true difference ``delta`` read with standard error ``se``.
    ``0.5`` for a tie (a coin flip); ``0`` for a non-zero difference read without noise."""
    from scipy.special import ndtr

    d = abs(float(delta))
    if se <= 0:
        return 0.5 if d == 0 else 0.0
    return float(ndtr(-d / se))


def _log_change(t, v, week) -> float:
    r = landmark_response(t, v, week)
    return float(np.log(max(1.0 + r, _LOG_FLOOR))) if np.isfinite(r) else float("nan")


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
    noise_cv: float | None = None  # per-measurement CV of the noisy-readout model, if any
    n_patients: int | None = None
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
        noisy = self.noise_cv is not None
        extra = ("expected_flips", "expected_discordant_pairs") if noisy else ()
        d = {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "reference_link": self.reference_link,
            "reference_ranking": list(self.reference_ranking),
            "total_pairs": self.total_pairs,
            "tier": self.tier,
            "rows": [
                {"week": r["week"], "discordant_pairs": r["discordant_pairs"],
                 "ranking": list(r["ranking"]), **{k: r[k] for k in extra}}
                for r in self.rows
            ],
        }
        if noisy:
            d["noise_cv"] = self.noise_cv
            d["n_patients"] = self.n_patients
        return d

    def to_json(self, *, indent: int = 2) -> str:
        from ._serialize import to_strict_json

        return to_strict_json(self.to_dict(), indent=indent)


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
    noise_cv: float | None = None,
    n_patients: int = 50,
) -> SurrogateTiming:
    """For each landmark week, rank a context's eligible TGI models by their early-surrogate
    response (most shrinkage = best) and count how many model pairs that ranking orders
    oppositely to the tail-aware **durable-benefit** reference ranking (median OS under the
    k_g survival link, which rewards slow regrowth).

    The earlier the landmark, the more the surrogate over-rewards deep-but-doomed early
    responders — so the discordance falls as the landmark moves later.

    ``noise_cv`` (opt-in) adds a measurement-noise model: each arm's readout is the mean of
    ``n_patients`` noisy log-changes (:func:`readout_se`). Each row then also reports
    ``expected_flips`` — the expected number of model pairs a trial's noisy readout orders
    opposite to its own noise-free ranking (precision) — and
    ``expected_discordant_pairs`` against the durable-benefit reference (fidelity under
    noise). Left ``None``, the analysis is unchanged."""
    se = None if noise_cv is None else readout_se(noise_cv, n_patients)
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
    tier = worst_tier([tr.tier for tr in cmp.included]) if cmp.included else "C"

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
        row = {
            "week": float(w),
            "ranking": ranking,
            "discordant_pairs": discordant_pairs(ranking, reference_ranking),
        }
        if se is not None:
            lc = {rid: _log_change(*trajectories[rid], w) for rid in models}
            flips = expected = 0.0
            for i, a in enumerate(reference_ranking):
                for b in reference_ranking[i + 1:]:
                    p_flip = ranking_flip_probability(lc[a] - lc[b], se)
                    flips += p_flip
                    # The reference prefers a; the readout agrees when a shrank more.
                    agrees = lc[a] < lc[b]
                    tie = lc[a] == lc[b]
                    expected += 0.5 if tie else (p_flip if agrees else 1.0 - p_flip)
            row["expected_flips"] = flips
            row["expected_discordant_pairs"] = expected
        st.rows.append(row)
    if se is not None:
        st.noise_cv, st.n_patients = float(noise_cv), int(n_patients)
    return st
