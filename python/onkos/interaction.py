"""Drug-combination interaction models — the interaction assumption as a model-
selection axis.

Oncology is overwhelmingly combination therapy, yet a composed survival forecast for
a combination silently depends on one unmeasured choice: *how do the two drugs'
effects combine?* Highest-single-agent, additive (the Bliss/Loewe null), or
synergistic — these give different predicted benefits from the *same* single-agent
activity, and the difference is routinely assumed rather than measured.

This module makes the interaction model a first-class, quantified model-selection
axis. It combines two single-agent effect magnitudes ``E_A, E_B`` (the kernel's
drug-effect scalar) under each declared interaction rule, feeds the result through the
*existing* TGI -> survival chain unchanged, and reports how much the predicted outcome
depends on which interaction you assumed. It is the combination-therapy analog of the
virtual-trial divergence view (across TGI models) and the kill-mechanism axis (across
kill mechanisms): the same "make the silent modeling choice visible" move one layer up.

**Synergy is an assumption, not a finding.** Onkos does not estimate the interaction
parameter from data — distinguishing synergy from additivity requires a combination
trial designed for it. ``psi`` is a *declared* input (default 0, the additive null);
a non-zero value carries a warning and never presents as a curated quantity. The
headline is the *spread* across interaction assumptions, not a single synergistic
answer.

Population / regimen level only. NOT a per-patient quantity, NOT a therapy ranking,
NOT a recommendation. Combination-vs-monotherapy curves are descriptive simulation of
what an interaction model predicts, never a choice between treatments. The underlying
TGI model's propagated tier governs and cannot be raised by any interaction assumption.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .load import Dataset
from .simulate import Trajectory, simulate

__all__ = [
    "INTERACTION_MODELS",
    "SYNERGY_IS_AN_ASSUMPTION",
    "combine_effects",
    "bliss_fraction",
    "simulate_combination",
    "compare_interactions",
    "InteractionComparison",
]

INTERACTION_MODELS = ("hsa", "additive", "greco")

# Printed wherever a non-zero interaction parameter appears, so the distinction the
# project depends on — assumption vs measurement — cannot be lost.
SYNERGY_IS_AN_ASSUMPTION = (
    "interaction parameter is a DECLARED assumption (not estimated from data); "
    "distinguishing synergy from additivity needs a combination trial designed for it"
)


# --------------------------------------------------------------------------- #
# Pure interaction math — the landmark-tested combination rules.              #
# Combine two non-negative single-agent effect magnitudes into one.           #
# --------------------------------------------------------------------------- #


def combine_effects(e_a: float, e_b: float, *, model: str = "additive", psi: float = 0.0) -> float:
    """Combine two single-agent effect magnitudes ``e_a, e_b >= 0`` into the combined
    effect ``E_AB`` under an interaction ``model``:

    * ``hsa``      — highest single agent: ``max(e_a, e_b)`` (the conservative null);
    * ``additive`` — Bliss-independence / effect-additive null: ``e_a + e_b``
      (for log-linear kill, Bliss independence *is* additive rates; see
      :func:`bliss_fraction`);
    * ``greco``    — interaction index: ``e_a + e_b + psi*sqrt(e_a*e_b)``, ``psi=0``
      additive, ``psi>0`` synergy, ``psi<0`` antagonism (clamped at 0).

    ``psi`` is a declared assumption, never fitted from the dataset.
    """
    a = max(float(e_a), 0.0)
    b = max(float(e_b), 0.0)
    if model == "hsa":
        return max(a, b)
    if model == "additive":
        return a + b
    if model == "greco":
        return max(0.0, a + b + float(psi) * np.sqrt(a * b))
    raise ValueError(f"unknown interaction model {model!r}; choose from {INTERACTION_MODELS}")


def bliss_fraction(e_a: float, e_b: float) -> float:
    """The Bliss-independence combined *fractional* effect over the two log-linear
    survival fractions: ``1 - (1-f_a)(1-f_b)`` with ``f = 1 - exp(-E)``.

    Equals ``1 - exp(-(e_a + e_b))`` — i.e. Bliss independence corresponds to *additive*
    effect magnitudes for a log-linear (exponential) kill process. Exposed so the
    identity is checkable, not asserted (landmark suite)."""
    fa = 1.0 - np.exp(-max(float(e_a), 0.0))
    fb = 1.0 - np.exp(-max(float(e_b), 0.0))
    return float(1.0 - (1.0 - fa) * (1.0 - fb))


# --------------------------------------------------------------------------- #
# Simulation bridge — feed the combined effect through the existing chain.    #
# --------------------------------------------------------------------------- #


def simulate_combination(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    effect_a: float,
    effect_b: float,
    interaction: str = "additive",
    psi: float = 0.0,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> Trajectory:
    """Simulate one combination regimen: combine ``effect_a`` and ``effect_b`` under
    ``interaction`` into a single drug effect and run it through ``record_id``'s
    existing TGI -> survival chain. Returns a standard :class:`~onkos.simulate.Trajectory`
    (tier and warnings from the underlying model, with the synergy-is-an-assumption note
    appended for a non-zero ``psi``).

    The combined effect is a *modeling assumption*, not a curated parameter — exactly
    the choice this module exists to make visible."""
    e_ab = combine_effects(effect_a, effect_b, model=interaction, psi=psi)
    tr = simulate(ds, record_id, context=context, drug_effect=e_ab, t=t,
                  survival_link=survival_link)
    if interaction == "greco" and psi != 0.0:
        tr.warnings = list(tr.warnings) + [f"synergy_is_an_assumption: {SYNERGY_IS_AN_ASSUMPTION}"]
    return tr


# --------------------------------------------------------------------------- #
# The divergence view — the interaction model as a model-selection axis.      #
# --------------------------------------------------------------------------- #


@dataclass
class InteractionComparison:
    record_id: str
    context: dict
    effect_a: float
    effect_b: float
    psi: float
    t: np.ndarray
    trajectories: dict = field(default_factory=dict)  # label -> Trajectory
    combined_effects: dict = field(default_factory=dict)  # label -> E_AB
    tier: str = "C"
    warnings: list = field(default_factory=list)

    def _curves(self, endpoint: str):
        return [tr.survival.get(endpoint) for tr in self.trajectories.values()
                if tr.survival.get(endpoint) is not None]

    def _divergence(self, endpoint: str) -> float:
        curves = self._curves(endpoint)
        if len(curves) < 2:
            return 0.0
        stacked = np.vstack(curves)
        return float(np.max(stacked.max(axis=0) - stacked.min(axis=0)))

    @property
    def os_divergence(self) -> float:
        """Max pointwise OS spread across the interaction models — how much the
        predicted survival depends on the interaction assumption alone."""
        return self._divergence("OS")

    @property
    def pfs_divergence(self) -> float:
        return self._divergence("PFS")

    @property
    def median_os(self) -> dict:
        return {k: tr.median_os for k, tr in self.trajectories.items()}

    @property
    def median_os_range(self) -> tuple | None:
        meds = [m for m in self.median_os.values() if m is not None]
        return (min(meds), max(meds)) if meds else None

    def to_dict(self) -> dict:
        """JSON-serializable result. Carries the clinical-use prohibition and the
        synergy-is-an-assumption note; the combined effects are recorded as the
        assumptions they are."""
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "record_id": self.record_id,
            "context": self.context,
            "effect_a": self.effect_a,
            "effect_b": self.effect_b,
            "psi": self.psi,
            "synergy_note": SYNERGY_IS_AN_ASSUMPTION,
            "tier": self.tier,
            "combined_effects": self.combined_effects,
            "median_os_weeks": self.median_os,
            "os_divergence": self.os_divergence,
            "pfs_divergence": self.pfs_divergence,
            "median_os_range": list(self.median_os_range) if self.median_os_range else None,
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def compare_interactions(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    effect_a: float,
    effect_b: float,
    psi: float = 0.5,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> InteractionComparison:
    """Simulate one combination (the same ``effect_a, effect_b``) under every
    interaction model — ``hsa``, ``additive``, and Greco synergy/antagonism at
    ``±psi`` — and quantify how much the predicted OS depends on the interaction
    assumption. The interaction-model analog of the virtual-trial divergence view.

    ``psi`` sets the magnitude of the synergy / antagonism bracket (a declared
    assumption); the additive case is always included as the null."""
    if t is None:
        t = np.linspace(0.0, 156.0, 313)
    t = np.asarray(t, dtype=float)

    specs = [("hsa", "hsa", 0.0), ("additive", "additive", 0.0)]
    if psi != 0.0:
        specs += [(f"greco+{abs(psi):g}", "greco", abs(psi)),
                  (f"greco-{abs(psi):g}", "greco", -abs(psi))]

    cmp = InteractionComparison(
        record_id=record_id, context=context or {}, effect_a=effect_a, effect_b=effect_b,
        psi=psi, t=t,
    )
    for label, model, p in specs:
        tr = simulate_combination(
            ds, record_id, context=context, effect_a=effect_a, effect_b=effect_b,
            interaction=model, psi=p, t=t, survival_link=survival_link,
        )
        cmp.trajectories[label] = tr
        cmp.combined_effects[label] = combine_effects(effect_a, effect_b, model=model, psi=p)
    # Tier is shared across interaction models (same underlying chain + context).
    cmp.tier = next(iter(cmp.trajectories.values())).tier
    cmp.warnings = [SYNERGY_IS_AN_ASSUMPTION] if psi != 0.0 else []
    return cmp
