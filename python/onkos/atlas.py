"""Model-selection atlas — one index and one per-context survey of every axis.

Across eighteen versions Onkos turned each silent modeling choice in the TGI → survival
chain into an explicit, quantified **model-selection axis**: which TGI model, which
resistance mechanism and origin, which bridge metric, which survival structure, which
exposure-response shape, which readout timing — and, finally, whether a trial could even
tell the models apart. Each lives in its own module with its own headline. This module is
the synthesis layer: a single declarative **registry** of the axes (the source of truth for
"what model-selection risk does Onkos quantify") and a single **per-context survey** that
runs each applicable axis and returns its native headline.

It is deliberately a *navigational survey*, not a variance decomposition. The axes are not
orthogonal, their headlines are in different units (weeks of median-OS spread, discordant
model pairs, required trial events), and the magnitudes depend on the operating point — so
the atlas reports each axis in its **own** unit, flags ``comparable = False``, and points to
:func:`onkos.model_selection_budget` for the rigorous, orthogonal two-factor decomposition.
The value is discoverability and a one-call "where is the model-selection risk for this
context?" map, with each axis's own module and CLI for the deep dive.

Population / trial level only. The atlas inherits the propagated tier of the context's
models and never moves it; it emits no individual quantity and no recommendation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .load import Dataset

__all__ = [
    "Axis",
    "AXES",
    "AtlasEntry",
    "Atlas",
    "model_selection_atlas",
]


@dataclass(frozen=True)
class Axis:
    """A model-selection axis: a silent modeling choice Onkos has made explicit."""

    key: str
    label: str
    varies: str  # the choice it varies
    finding: str  # the canonical one-line finding
    module: str  # onkos.<module>
    cli: str  # the CLI command
    scope: str  # "single-agent" | "combination"
    version: str


# The registry — the single source of truth for the model-selection axes. The README
# cheat-sheet mirrors this; the landmark suite checks it is complete and well-formed.
AXES: tuple = (
    Axis("tgi_model", "TGI model", "which tumor-growth model (compare().included)",
         "matched in-context models can imply median OS from ~54 to ~94 wk",
         "onkos.compare", "onkos simulate --compare", "single-agent", "v0.13"),
    Axis("survival_link", "survival link / bridge metric",
         "which on-treatment metric drives the hazard (week-8 / k_g / integrated burden)",
         "the metric choice inverts the resistance-model ranking and re-ranks a complete responder",
         "onkos.simulate", "onkos simulate (survival_link=)", "single-agent", "v0.25"),
    Axis("survival_structure", "survival structure (PH vs joint)",
         "two-stage proportional hazards vs the joint current-value link",
         "the joint link's hazard ratio rises 10x-255x as a clone regrows (non-PH); it inverts week-8",
         "onkos.joint", "onkos joint", "single-agent", "v0.34"),
    Axis("exposure_response", "exposure-response shape",
         "which ER shape (Emax / power / sigmoid) maps dose to effect",
         "invisible at the studied dose, a ~19 wk OS swing on de-escalation",
         "onkos.dose_response", "onkos dose-response", "single-agent", "v0.36"),
    Axis("readout_timing", "early-surrogate readout timing",
         "when the surrogate is read (the ctDNA push to week 2-4 vs RECIST week 8)",
         "earlier readout trades fidelity to durable benefit (9/10 -> 3/10 discordant)",
         "onkos.early_surrogate", "onkos early-surrogate", "single-agent", "v0.37"),
    Axis("model_discriminability", "model discriminability (meta-axis)",
         "whether a trial can tell the competing models apart at all",
         "the resistance mechanism/origin needs 1e4-1e5 events: practically unidentifiable",
         "onkos.discriminability", "onkos discriminability", "single-agent", "v0.38"),
    Axis("additivity_reference", "combination additivity reference",
         "which 'no-interaction' null (HSA / Bliss / Loewe) for a drug combination",
         "the reference moves combined OS; only Loewe passes the sham-combination test",
         "onkos.interaction", "onkos loewe", "combination", "v0.35"),
)

# The rigorous, orthogonal companion the atlas is NOT (deliberately): the variance budget.
BUDGET_REFERENCE = (
    "onkos.model_selection_budget decomposes Var(median OS) orthogonally over the TGI-model "
    "and survival-link factors (+ parameter noise) — the rigorous two-factor partition; the "
    "atlas is the broader, non-orthogonal survey across all axes."
)


@dataclass
class AtlasEntry:
    key: str
    label: str
    headline: float | int | None
    unit: str
    detail: str

    def to_dict(self) -> dict:
        return {"key": self.key, "label": self.label, "headline": self.headline,
                "unit": self.unit, "detail": self.detail}


@dataclass
class Atlas:
    """A per-context survey of the single-agent model-selection axes."""

    context: dict
    entries: list = field(default_factory=list)
    tier: str = "C"
    comparable: bool = False
    note: str = BUDGET_REFERENCE
    clinical_use: str = CLINICAL_USE

    def get(self, key: str) -> AtlasEntry | None:
        return next((e for e in self.entries if e.key == key), None)

    @property
    def os_swing_axes(self) -> list:
        """The axes whose headline is a median-OS spread in weeks — loosely comparable as
        'weeks of OS riding on this one choice' (NOT orthogonal; see the budget)."""
        return [e for e in self.entries if e.unit == "weeks (median-OS spread)"]

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "tier": self.tier,
            "comparable": self.comparable,
            "note": self.note,
            "entries": [e.to_dict() for e in self.entries],
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def _os_range(curves) -> float | None:
    meds = [m for m in curves if m is not None]
    return float(max(meds) - min(meds)) if len(meds) >= 2 else None


def model_selection_atlas(
    ds: Dataset,
    *,
    context: dict,
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
) -> Atlas:
    """Run each applicable single-agent model-selection axis for ``context`` and return its
    native headline — a one-call survey of where the model-selection risk lies.

    Each axis is reported in its own unit (the atlas is a survey, not a decomposition;
    ``comparable = False``). Combination axes (the additivity reference) need a regimen and
    are catalogued in :data:`AXES` but not run here."""
    from .budget import eligible_survival_links
    from .compare import compare
    from .discriminability import model_discriminability
    from .dose_response import compare_er_extrapolation
    from .early_surrogate import surrogate_timing_fidelity
    from .joint import joint_survival
    from .simulate import simulate

    if t is None:
        t = np.linspace(0.0, 312.0, 625)
    t = np.asarray(t, dtype=float)

    cmp = compare(ds, purpose="tgi", context=context, drug_effect=drug_effect, t=t)
    models = [tr.record_id for tr in cmp.included]
    tier = cmp.included[0].tier if cmp.included else "C"
    atlas = Atlas(context=context, tier=tier)

    # 1) TGI-model axis — OS divergence across the eligible models.
    rng = cmp.median_os_range
    atlas.entries.append(AtlasEntry(
        "tgi_model", "TGI model",
        round(rng[1] - rng[0], 1) if rng else None, "weeks (median-OS spread)",
        f"{len(models)} eligible models" + (f", median OS {rng[0]:.0f}-{rng[1]:.0f} wk" if rng else ""),
    ))

    # 2) survival-link axis — the worst-case median-OS spread across eligible OS links, over
    #    all models (max-over-models, so a complete responder with n/r links cannot mask it).
    links = eligible_survival_links(ds, context, "OS")
    if models and len(links) >= 2:
        per_model = []
        for m in models:
            meds = [simulate(ds, m, context=context, drug_effect=drug_effect, t=t,
                             survival_link=lk).median_os for lk in links]
            r = _os_range(meds)
            if r is not None:
                per_model.append(r)
        atlas.entries.append(AtlasEntry(
            "survival_link", "survival link / bridge metric",
            round(max(per_model), 1) if per_model else None, "weeks (median-OS spread)",
            f"{len(links)} eligible OS links (week-8 / k_g / …); worst-case model swing",
        ))

    # 3) survival-structure axis — the largest joint-vs-two-stage median gap over the models,
    #    and the peak non-proportional hazard ratio.
    gaps, peak_phv = [], 0.0
    for m in models:
        js = joint_survival(ds, m, context=context, drug_effect=drug_effect, t=t, alpha=1.0)
        if np.isfinite(js.ph_violation):
            peak_phv = max(peak_phv, js.ph_violation)
        if js.median_os is not None and js.two_stage_median_os is not None:
            gaps.append(abs(js.median_os - js.two_stage_median_os))
    if models:
        atlas.entries.append(AtlasEntry(
            "survival_structure", "survival structure (PH vs joint)",
            round(max(gaps), 1) if gaps else None, "weeks (median-OS spread)",
            f"joint hazard ratio peaks at {peak_phv:.0f}x the week-8 value (non-proportional)",
        ))

    # 4) exposure-response axis — worst-case OS divergence across ER shapes on dose
    #    extrapolation, over the models.
    if models:
        er_div = max(compare_er_extrapolation(ds, m, context=context, c_ref=150.0,
                                              e_ref=1.0).max_os_divergence for m in models)
        atlas.entries.append(AtlasEntry(
            "exposure_response", "exposure-response shape", round(er_div, 1),
            "weeks (median-OS spread)", "max over the dose grid (0 at the studied dose)",
        ))

    # 5) readout-timing axis — early-vs-late surrogate discordance with durable benefit.
    st = surrogate_timing_fidelity(ds, context=context, drug_effect=drug_effect, t=t)
    atlas.entries.append(AtlasEntry(
        "readout_timing", "early-surrogate readout timing", st.earliest_discordance,
        f"discordant model pairs / {st.total_pairs}",
        f"earliest readout {st.earliest_discordance}/{st.total_pairs} vs latest {st.latest_discordance}/{st.total_pairs}",
    ))

    # 6) discriminability meta-axis — model pairs a realistic trial cannot resolve.
    md = model_discriminability(ds, context=context, drug_effect=drug_effect, t=t)
    atlas.entries.append(AtlasEntry(
        "model_discriminability", "model discriminability (meta-axis)", md.n_indistinguishable,
        f"indistinguishable model pairs / {len(md.pairs)}",
        "pairs needing an infeasible trial to tell apart (week-8 OS)",
    ))

    return atlas
