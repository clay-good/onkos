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

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .export.registry import get_kernel, kernel_values
from .load import Dataset
from .simulate import Trajectory, simulate

__all__ = [
    "INTERACTION_MODELS",
    "ADDITIVITY_REFERENCES",
    "SYNERGY_IS_AN_ASSUMPTION",
    "combine_effects",
    "bliss_fraction",
    "simulate_combination",
    "compare_interactions",
    "InteractionComparison",
    "ERCurve",
    "er_curve",
    "loewe_effect",
    "combine_doses",
    "simulate_dose_combination",
    "compare_additivity_references",
    "AdditivityComparison",
]

INTERACTION_MODELS = ("hsa", "additive", "greco")

# The dose-level "no-interaction" references (the additivity null is itself a choice):
#   hsa   — highest single agent (the conservative null);
#   bliss — Bliss independence / effect-additive (add the two effect magnitudes);
#   loewe — Loewe dose-additivity (the isobole over the dose-response curves), the only
#           one that satisfies the sham-combination identity (a drug with itself is
#           exactly additive). Requires the exposure-response curves, not just the effects.
ADDITIVITY_REFERENCES = ("hsa", "bliss", "loewe")

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


# --------------------------------------------------------------------------- #
# Dose-level Loewe additivity — the additivity REFERENCE as a model choice.   #
#                                                                             #
# The effect-level nulls above (hsa/additive) combine the two effect          #
# MAGNITUDES. Loewe additivity instead combines the two DOSES via each drug's #
# dose-response (exposure-response) curve: the combination (d_a, d_b) produces #
# effect E iff  d_a/D_a(E) + d_b/D_b(E) = 1,  where D_x(E) is the dose of      #
# drug x ALONE producing effect E (the inverse ER curve). Loewe is the only   #
# null that satisfies the SHAM-COMBINATION identity — a drug combined with     #
# itself is exactly additive: d_a/D(E)+d_b/D(E)=1 -> D(E)=d_a+d_b ->           #
# E=f(d_a+d_b). Bliss/effect-additivity fails this for any saturating curve.  #
# --------------------------------------------------------------------------- #


@dataclass
class ERCurve:
    """An exposure-response curve as a forward map (dose -> effect), its analytic
    inverse (effect -> dose, for effect < emax), and the effect ceiling ``emax``
    (``inf`` for an unbounded power curve). Built from a record by :func:`er_curve`;
    constructible from plain callables so the Loewe core is landmark-testable alone."""

    forward: Callable[[float], float]
    inverse: Callable[[float], float]
    emax: float


def er_curve(ds: Dataset, er_id: str) -> ERCurve:
    """Build an :class:`ERCurve` from an exposure-response record, with the analytic
    inverse for the three curated ER kernels (Emax, sigmoid-Emax, power)."""
    rec = ds[er_id]
    spec = get_kernel(rec)
    v = kernel_values(rec)
    name = spec.name
    if name == "er_emax":
        emax, ec50 = float(v["Emax"]), float(v["EC50"])

        def fwd(d):
            return emax * d / (ec50 + d)

        def inv(e):  # valid 0 <= e < emax
            return ec50 * e / (emax - e)

        return ERCurve(forward=fwd, inverse=inv, emax=emax)
    if name == "er_sigmoid_emax":
        emax, ec50, g = float(v["Emax"]), float(v["EC50"]), float(v["gamma"])

        def fwd(d):
            return emax * d**g / (ec50**g + d**g)

        def inv(e):
            return ec50 * (e / (emax - e)) ** (1.0 / g)

        return ERCurve(forward=fwd, inverse=inv, emax=emax)
    if name == "er_power":
        slope, theta = float(v["slope"]), float(v["theta"])

        def fwd(d):
            return slope * d**theta

        def inv(e):  # unbounded: no finite emax
            return (e / slope) ** (1.0 / theta)

        return ERCurve(forward=fwd, inverse=inv, emax=float("inf"))
    raise ValueError(f"no analytic inverse for ER kernel {name!r} (record {er_id})")


def loewe_effect(dose_a: float, dose_b: float, *, curve_a: ERCurve, curve_b: ERCurve) -> float:
    """The Loewe dose-additive combined effect of ``dose_a`` of drug A and ``dose_b`` of
    drug B, solving the isobole ``d_a/D_a(E) + d_b/D_b(E) = 1`` for ``E``.

    Single-agent limits (one dose 0) return the other drug's effect exactly. If the
    combined effect would exceed the shared effect ceiling it is clamped there (the curves
    cannot jointly express more than ``min(emax_a, emax_b)``)."""
    da, db = max(float(dose_a), 0.0), max(float(dose_b), 0.0)
    if db <= 0.0:
        return float(curve_a.forward(da))
    if da <= 0.0:
        return float(curve_b.forward(db))
    ceiling = min(curve_a.emax, curve_b.emax)
    hi = ceiling - 1e-9 if np.isfinite(ceiling) else max(curve_a.forward(da), curve_b.forward(db)) * 10.0

    def isobole(e: float) -> float:
        # d_a/D_a(E) + d_b/D_b(E) - 1; strictly decreasing in E (D_x increases with E).
        return da / curve_a.inverse(e) + db / curve_b.inverse(e) - 1.0

    lo = 1e-9
    if isobole(lo) <= 0.0:  # even an infinitesimal effect already over-satisfies (degenerate)
        return 0.0
    if isobole(hi) >= 0.0:  # combination demands more effect than the shared ceiling allows
        return float(hi)
    # Bisection (monotone isobole) — no scipy dependency, deterministic.
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if isobole(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def combine_doses(
    ds: Dataset, dose_a: float, dose_b: float, *, er_a: str, er_b: str, reference: str = "loewe"
) -> float:
    """Combined effect magnitude from two DOSES under a dose-level additivity
    ``reference`` (``hsa`` / ``bliss`` / ``loewe``), using the ER curves ``er_a``/``er_b``.

    ``hsa`` and ``bliss`` reduce to the effect-level nulls of :func:`combine_effects`
    evaluated at the single-agent effects ``f_a(d_a)``, ``f_b(d_b)``; ``loewe`` is the
    dose-additive isobole. The reference is a declared modeling choice, never fitted."""
    ca, cb = er_curve(ds, er_a), er_curve(ds, er_b)
    e_a, e_b = float(ca.forward(dose_a)), float(cb.forward(dose_b))
    if reference == "hsa":
        return max(e_a, e_b)
    if reference == "bliss":
        return e_a + e_b
    if reference == "loewe":
        return loewe_effect(dose_a, dose_b, curve_a=ca, curve_b=cb)
    raise ValueError(f"unknown additivity reference {reference!r}; choose from {ADDITIVITY_REFERENCES}")


def simulate_dose_combination(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    dose_a: float,
    dose_b: float,
    er_a: str,
    er_b: str,
    reference: str = "loewe",
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> Trajectory:
    """Combine two doses under a dose-level additivity ``reference`` into one drug effect
    and run it through ``record_id``'s existing TGI -> survival chain."""
    e_ab = combine_doses(ds, dose_a, dose_b, er_a=er_a, er_b=er_b, reference=reference)
    return simulate(ds, record_id, context=context, drug_effect=e_ab, t=t, survival_link=survival_link)


@dataclass
class AdditivityComparison:
    """OS across the dose-level additivity references for one dose pair."""

    record_id: str
    context: dict
    dose_a: float
    dose_b: float
    er_a: str
    er_b: str
    t: np.ndarray
    trajectories: dict = field(default_factory=dict)  # reference -> Trajectory
    combined_effects: dict = field(default_factory=dict)  # reference -> E_AB
    tier: str = "C"
    warnings: list = field(default_factory=list)

    @property
    def os_divergence(self) -> float:
        curves = [tr.survival.get("OS") for tr in self.trajectories.values()
                  if tr.survival.get("OS") is not None]
        if len(curves) < 2:
            return 0.0
        stacked = np.vstack(curves)
        return float(np.max(stacked.max(axis=0) - stacked.min(axis=0)))

    @property
    def median_os(self) -> dict:
        return {k: tr.median_os for k, tr in self.trajectories.items()}

    @property
    def median_os_range(self) -> tuple | None:
        meds = [m for m in self.median_os.values() if m is not None]
        return (min(meds), max(meds)) if meds else None

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "record_id": self.record_id,
            "context": self.context,
            "dose_a": self.dose_a,
            "dose_b": self.dose_b,
            "er_a": self.er_a,
            "er_b": self.er_b,
            "tier": self.tier,
            "combined_effects": self.combined_effects,
            "median_os_weeks": self.median_os,
            "os_divergence": self.os_divergence,
            "median_os_range": list(self.median_os_range) if self.median_os_range else None,
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def compare_additivity_references(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    dose_a: float,
    dose_b: float,
    er_a: str,
    er_b: str,
    t: np.ndarray | None = None,
    survival_link: str | None = None,
) -> AdditivityComparison:
    """Simulate the SAME dose pair under every additivity reference — ``hsa``, ``bliss``
    (effect-additive), and ``loewe`` (dose-additive) — and quantify how much the predicted
    OS depends on which "no-interaction" reference you assume. The dose-level analog of
    :func:`compare_interactions`; the headline is the spread, with Loewe the self-consistent
    reference (it alone passes the sham-combination test)."""
    if t is None:
        t = np.linspace(0.0, 156.0, 313)
    t = np.asarray(t, dtype=float)
    cmp = AdditivityComparison(
        record_id=record_id, context=context or {}, dose_a=dose_a, dose_b=dose_b,
        er_a=er_a, er_b=er_b, t=t,
    )
    for ref in ADDITIVITY_REFERENCES:
        tr = simulate_dose_combination(
            ds, record_id, context=context, dose_a=dose_a, dose_b=dose_b,
            er_a=er_a, er_b=er_b, reference=ref, t=t, survival_link=survival_link,
        )
        cmp.trajectories[ref] = tr
        cmp.combined_effects[ref] = combine_doses(
            ds, dose_a, dose_b, er_a=er_a, er_b=er_b, reference=ref
        )
    cmp.tier = next(iter(cmp.trajectories.values())).tier
    return cmp
