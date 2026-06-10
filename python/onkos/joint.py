"""Joint longitudinal–survival modeling — the current-value link.

Every survival link Onkos has shipped is **two-stage**: collapse the simulated tumor
trajectory to one scalar covariate — the week-8 change (v0.12), the growth-rate constant
``k_g`` (v0.25), the integrated burden (v0.33) — then apply a parametric/Cox baseline with
that *static* covariate. A static covariate means a **proportional hazard**: the hazard
ratio between two tumors is constant over time. That is exactly the assumption the joint
longitudinal–survival model relaxes.

This module adds the **current-value** association structure (the canonical joint-model
link): the instantaneous hazard tracks the *current* tumor size,

    λ(t) = λ₀(t) · exp(α · log(v(t) / y0)),     S(t) = exp(-∫₀ᵗ λ(u) du)

with a Weibull baseline hazard λ₀ and an association coefficient ``α`` between log tumor
size and log hazard. The cumulative hazard is integrated as a Stieltjes sum against the
**analytic** baseline cumulative hazard ``H₀(t) = (t/scale)^shape``, so the implementation
is exact in the two limits that make it a strict generalization:

- a tumor held at baseline (``v ≡ y0``) recovers the Weibull **baseline** survival exactly;
- a tumor held at a constant ``c·y0`` recovers the **two-stage** Weibull-PH curve exactly,
  with covariate ``x = log c`` and ``β = α`` — i.e. the v0.33 burden link is the
  constant-trajectory special case of this model.

The payload: for a shrink-then-regrow tumor the hazard ratio is **time-varying** — low
while the tumor is small, then rising as a resistant clone regrows — a **non-proportional
hazard** none of the two-stage links can represent. So "two-stage vs joint" is a
model-selection axis at the survival-link layer, and the disagreement is concentrated in
exactly the resistance models whose tail the week-8 surrogate is blind to.

Population / trial level only. ``α`` and the baseline are **declared, illustrative**
parameters, never fitted here (a real joint model estimates them jointly from data; Onkos
simulates the structure forward). NOT an individual prognosis, NOT a therapy ranking. A
joint analysis never moves a tier; it inherits the propagated tier of the trajectory it
summarizes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .load import Dataset
from .simulate import median_survival, simulate

__all__ = [
    "current_value_survival",
    "JointSurvival",
    "joint_survival",
    "JointComparison",
    "compare_joint_vs_two_stage",
]

_FLOOR = 1e-3  # detection-limit floor on relative size (a complete response), as in metrics.py


# --------------------------------------------------------------------------- #
# Pure core — population survival under a time-varying (current-value) hazard  #
# ratio. Landmark-tested in isolation (no record needed).                      #
# --------------------------------------------------------------------------- #


def current_value_survival(
    t: np.ndarray, hazard_ratio: np.ndarray, *, shape: float, scale: float
) -> np.ndarray:
    """Population survival ``S(t) = exp(-∫₀ᵗ λ₀(u)·HR(u) du)`` for a Weibull baseline
    hazard ``λ₀`` and a **time-varying** hazard ratio ``HR(u)``.

    The integral is a Stieltjes sum against the analytic baseline cumulative hazard
    ``H₀(t) = (t/scale)^shape`` (``∫ λ₀ du`` in closed form), with the trapezoid rule on
    ``HR``. A constant ``HR`` therefore telescopes to ``HR·H₀(t)`` exactly — the two-stage
    Weibull-PH curve — so this is a strict generalization, not an approximation, of the
    proportional-hazards link."""
    t = np.asarray(t, dtype=float)
    hr = np.asarray(hazard_ratio, dtype=float)
    h0 = (t / scale) ** shape  # analytic baseline cumulative hazard
    incr = (hr[1:] + hr[:-1]) * 0.5 * np.diff(h0)
    cum = np.concatenate([[0.0], np.cumsum(incr)])
    return np.exp(-cum)


# --------------------------------------------------------------------------- #
# Binding to a record + context                                               #
# --------------------------------------------------------------------------- #


def _default_weibull_os_link(ds: Dataset, tumor_type, line):
    """The context's default parametric (Weibull) OS link — the baseline hazard the joint
    model builds on. Non-default links (Cox/k_g/burden) are not used as the baseline."""
    for r in ds:
        if r.purpose != "survival_link" or r.kernel != "survival_weibull_ph":
            continue
        dc = r.derivation_context
        if (
            dc
            and dc.tumor_type == tumor_type
            and (line is None or dc.line_of_therapy == line)
            and r.structure.get("default", True)
            and r.structure.get("endpoint", "OS") == "OS"
        ):
            return r
    return None


@dataclass
class JointSurvival:
    """A current-value joint-model survival curve beside its two-stage counterpart."""

    record_id: str
    t: np.ndarray
    alpha: float
    os_curve: np.ndarray  # current-value (joint) survival
    two_stage_curve: np.ndarray  # the default week-8 two-stage survival, for comparison
    hazard_ratio: np.ndarray  # the time-varying HR(t) = exp(alpha·log(v/y0))
    tier: str
    warnings: list[str] = field(default_factory=list)
    clinical_use: str = CLINICAL_USE

    @property
    def median_os(self) -> float | None:
        return median_survival(self.t, self.os_curve)

    @property
    def two_stage_median_os(self) -> float | None:
        return median_survival(self.t, self.two_stage_curve)

    def hazard_ratio_at(self, week: float) -> float:
        """The current-value hazard ratio at ``week`` (1.0 ⇒ tumor at baseline size)."""
        return float(np.interp(week, self.t, self.hazard_ratio))

    @property
    def ph_violation(self) -> float:
        """How far the hazard ratio departs from proportionality over the horizon:
        ``HR(end) / HR(8 wk)``. Exactly 1 for a constant-size (proportional) trajectory;
        ``>> 1`` for a regrowing tumor (the hazard rises in the tail), ``< 1`` for a
        deepening responder."""
        hr8 = float(np.interp(8.0, self.t, self.hazard_ratio))
        return float(self.hazard_ratio[-1] / hr8) if hr8 > 0 else float("inf")


def joint_survival(
    ds: Dataset,
    record_id: str,
    *,
    context: dict,
    drug_effect: float | None = 1.0,
    exposure=None,
    t: np.ndarray | None = None,
    alpha: float = 1.0,
    baseline_link: str | None = None,
    floor: float = _FLOOR,
) -> JointSurvival:
    """Current-value joint-model OS for ``record_id`` in ``context``.

    The tumor trajectory and its two-stage (default week-8) OS come from
    :func:`onkos.simulate`; the Weibull baseline hazard (``shape``/``scale``) is read from
    the context's default Weibull OS link (or ``baseline_link``); the instantaneous hazard
    ratio is ``exp(alpha·log(max(v/y0, floor)))``. The propagated tier and transport
    warnings ride through unchanged from the simulation."""
    if t is None:
        t = np.linspace(0.0, 260.0, 521)
    t = np.asarray(t, dtype=float)

    tr = simulate(ds, record_id, context=context, drug_effect=drug_effect, exposure=exposure, t=t)
    tumor_type = context.get("tumor_type")
    line = context.get("line") or context.get("line_of_therapy")

    link = ds[baseline_link] if baseline_link else _default_weibull_os_link(ds, tumor_type, line)
    if link is None:
        raise ValueError(
            f"no default Weibull OS link for context {context!r}; the joint model needs a "
            "parametric baseline hazard (pass baseline_link=...)"
        )
    shape = float(link["weibull_shape"].central)
    scale = float(link["weibull_scale"].central)

    y0 = float(tr.tumor_size[0])
    rel = np.maximum(tr.tumor_size / y0, floor)
    hr = np.exp(alpha * np.log(rel))
    joint_curve = current_value_survival(t, hr, shape=shape, scale=scale)

    return JointSurvival(
        record_id=record_id,
        t=t,
        alpha=float(alpha),
        os_curve=joint_curve,
        two_stage_curve=tr.os_curve if tr.os_curve is not None else np.full_like(t, np.nan),
        hazard_ratio=hr,
        tier=tr.tier,
        warnings=list(tr.warnings),
    )


@dataclass
class JointComparison:
    """Two-stage vs joint OS across a context's eligible TGI models."""

    context: dict
    alpha: float
    rows: list[dict]  # one per model: record_id, two_stage_median, joint_median, ph_violation
    clinical_use: str = CLINICAL_USE

    @property
    def rank_discordant_pairs(self) -> int:
        """Model pairs the two link structures order oppositely (a non-responder-robust
        count: only pairs both links resolve to a finite median, ordered strictly)."""
        finite = [
            r for r in self.rows
            if r["two_stage_median"] is not None and r["joint_median"] is not None
        ]
        n = 0
        for i in range(len(finite)):
            for j in range(i + 1, len(finite)):
                a, b = finite[i], finite[j]
                ts = a["two_stage_median"] - b["two_stage_median"]
                jt = a["joint_median"] - b["joint_median"]
                if ts * jt < 0:
                    n += 1
        return n

    @property
    def max_ph_violation(self) -> float:
        vals = [r["ph_violation"] for r in self.rows if np.isfinite(r["ph_violation"])]
        return max(vals) if vals else float("nan")


def compare_joint_vs_two_stage(
    ds: Dataset,
    *,
    context: dict,
    purpose: str = "tgi",
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
    alpha: float = 1.0,
) -> JointComparison:
    """For every eligible TGI model in ``context``, the two-stage (week-8) median OS, the
    joint (current-value) median OS, and the per-model PH-violation ``HR(end)/HR(8wk)``.

    The headline is twofold: the joint link **re-ranks** the models relative to the
    two-stage surrogate (because it weights the regrowth tail), and it exposes the
    **non-proportional hazard** (``ph_violation ≫ 1``) for the resistance models the
    two-stage surrogate cannot represent."""
    from .compare import compare

    if t is None:
        t = np.linspace(0.0, 260.0, 521)
    cmp = compare(ds, purpose=purpose, context=context, drug_effect=drug_effect, t=t)
    rows: list[dict] = []
    for tr in cmp.included:
        js = joint_survival(
            ds, tr.record_id, context=context, drug_effect=drug_effect, t=t, alpha=alpha
        )
        rows.append(
            {
                "record_id": tr.record_id,
                "two_stage_median": js.two_stage_median_os,
                "joint_median": js.median_os,
                "ph_violation": js.ph_violation,
                "tier": js.tier,
            }
        )
    return JointComparison(context=context, alpha=float(alpha), rows=rows)
