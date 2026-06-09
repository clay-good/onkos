"""Model-selection uncertainty & model averaging — the third uncertainty axis.

The virtual-trial divergence view (:mod:`onkos.compare`) shows *that* the
eligible published TGI models disagree and by how much. This module is the
inferential completion of that view: it splits the total predictive uncertainty
of a composed survival forecast into

* **WITHIN** — parameter (IIV) noise, the Axis-1 uncertainty that a bigger trial
  could shrink, obtained per model from :func:`onkos.uncertainty.ensemble_samples`;
* **BETWEEN** — model-selection disagreement, the irreducible Axis-3 risk that
  more data on any *one* model cannot resolve,

via the law of total variance, and combines the eligible models into a
model-averaged trajectory whose weights are *declared* and whose disagreement
travels with it. The single headline number is::

    model_selection_fraction = BETWEEN / (WITHIN + BETWEEN)  in [0, 1]

Guardrails (spec §5, §9 of the research spec):

* Averaging **cannot raise a tier** — the averaged tier is the worst included.
* Only ``Comparison.included`` (in-context) models are averaged; transported /
  failure-mode models stay excluded under every scheme.
* The point estimate is structurally inseparable from its
  ``model_selection_fraction`` — a result without it is a test failure.
* Weights are **forecast-combination weights** (Bates–Granger), explicitly *not*
  posterior model probabilities: the Onkos models are fit to different trials, so
  a posterior model probability is not identifiable and would be invented.

Population / trial level only. No individual prognosis, no therapy ranking.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .models import TIER_ORDER
from .tiers import worst_tier
from .uncertainty import ensemble_samples

__all__ = [
    "ModelAverage",
    "SCHEMES",
    "WEIGHTS_ARE_COMBINATION_NOT_POSTERIOR",
    "tier_scores",
    "compute_weights",
    "decompose",
    "average_curve",
]

# The weights are combination weights, not model posteriors — printed wherever
# weights appear so the distinction the project depends on cannot be lost.
WEIGHTS_ARE_COMBINATION_NOT_POSTERIOR = (
    "forecast-combination weights (Bates-Granger), NOT posterior model probabilities"
)

SCHEMES = ("equal", "tier", "evidence")

# Tier score: A:B:C:D = 8:4:2:1, i.e. A:B:C = 4:2:1 (research spec §3). A better-
# validated model speaks louder, but only by a *declared* factor, never a fitted one.
_TIER_SCORE = {t: float(2 ** (3 - rank)) for t, rank in TIER_ORDER.items()}


def tier_scores() -> dict:
    """The declared tier->weight scores (a convention, not a probability)."""
    return dict(_TIER_SCORE)


# --------------------------------------------------------------------------- #
# Pure combination math (no simulation) — the landmark-tested estimator core.  #
# --------------------------------------------------------------------------- #


def _normalize(w: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    s = w.sum()
    if not np.isfinite(s) or s <= 0:
        return np.full(w.size, 1.0 / w.size)
    return w / s


def compute_weights(scheme: str, tiers, evidence) -> tuple:
    """Return ``(weights, notes)`` for a scheme over the eligible models.

    ``tiers`` is a per-model tier string; ``evidence`` is a per-model external
    C-index (or ``None``). Weights are non-negative and sum to 1.
    """
    tiers = list(tiers)
    m = len(tiers)
    notes: list[str] = []
    if m == 0:
        return np.array([]), notes
    if scheme == "equal":
        return _normalize(np.ones(m)), notes
    if scheme == "tier":
        return _normalize(np.array([_TIER_SCORE.get(t, 1.0) for t in tiers])), notes
    if scheme == "evidence":
        raw = np.array([max(0.0, e - 0.5) if e is not None else 0.0 for e in evidence])
        if not np.any(raw > 0):
            notes.append(
                "evidence weighting unavailable (no eligible model has an external "
                "C-index > 0.5); fell back to equal weights"
            )
            return _normalize(np.ones(m)), notes
        if any(e is None for e in evidence):
            notes.append(
                "models without a recorded external C-index get zero data-weight under "
                "the evidence scheme"
            )
        return _normalize(raw), notes
    raise ValueError(f"unknown weighting scheme {scheme!r}; choose from {SCHEMES}")


def decompose(means, variances, weights) -> tuple:
    """Law of total variance over a mixture of models.

    Returns ``(point, within, between, fraction)`` where
    ``point = Σ w_m·mean_m``, ``within = Σ w_m·var_m`` (parameter noise),
    ``between = Σ w_m·(mean_m − point)²`` (model-selection), and
    ``fraction = between / (within + between)``.
    """
    means = np.asarray(means, dtype=float)
    variances = np.asarray(variances, dtype=float)
    w = _normalize(weights)
    point = float(np.sum(w * means))
    within = float(np.sum(w * variances))
    between = float(np.sum(w * (means - point) ** 2))
    total = within + between
    fraction = float(between / total) if total > 0 else 0.0
    return point, within, between, fraction


def average_curve(curves, weights, within_curves=None) -> tuple:
    """Pointwise convex combination of per-model survival curves.

    ``curves`` is an ``(M, T)`` array of per-model mean survival functions
    ``S_m(t)``. Returns ``(S_bar, between_t, within_t)`` — the averaged curve, the
    pointwise between-model variance ``Σ w_m (S_m − S_bar)²``, and (if
    ``within_curves`` given) the pointwise within-model variance.
    """
    curves = np.asarray(curves, dtype=float)
    w = _normalize(weights)[:, None]
    s_bar = np.sum(w * curves, axis=0)
    between_t = np.sum(w * (curves - s_bar) ** 2, axis=0)
    within_t = None
    if within_curves is not None:
        within_t = np.sum(w * np.asarray(within_curves, dtype=float), axis=0)
    return s_bar, between_t, within_t


# --------------------------------------------------------------------------- #
# The result object.                                                          #
# --------------------------------------------------------------------------- #


@dataclass
class ModelAverage:
    target: str
    endpoint: str
    scheme: str
    n: int
    t: np.ndarray
    point: float
    curve: np.ndarray  # S_bar(t)
    within_var: float
    between_var: float
    model_selection_fraction: float
    weights: dict  # record_id -> w_m
    per_model_mean: dict  # record_id -> E[Q | m]
    tier: str
    weight_sensitivity: float
    between_curve: np.ndarray  # pointwise between-model variance of S_bar(t)
    within_curve: np.ndarray  # pointwise within-model variance
    scheme_points: dict = field(default_factory=dict)  # scheme -> point
    warnings: list = field(default_factory=list)

    @property
    def between_band(self) -> np.ndarray:
        """Pointwise between-model standard deviation around ``curve``."""
        return np.sqrt(self.between_curve)

    def to_dict(self, *, include_curves: bool = False) -> dict:
        """JSON-serializable result. The clinical-use prohibition and the
        ``model_selection_fraction`` are always present — the point estimate is
        never reported without the disagreement that qualifies it."""
        d = {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "onkos:modelSelectionUncertainty": {
                "fraction": self.model_selection_fraction,
                "scheme": self.scheme,
            },
            "target": self.target,
            "endpoint": self.endpoint,
            "scheme": self.scheme,
            "weights_meaning": WEIGHTS_ARE_COMBINATION_NOT_POSTERIOR,
            "n": self.n,
            "point": self.point,
            "within_var": self.within_var,
            "between_var": self.between_var,
            "model_selection_fraction": self.model_selection_fraction,
            "weights": self.weights,
            "per_model_mean": self.per_model_mean,
            "tier": self.tier,
            "weight_sensitivity": self.weight_sensitivity,
            "scheme_points": self.scheme_points,
            "warnings": list(self.warnings),
        }
        if include_curves:
            d["t"] = self.t.tolist()
            d["curve"] = self.curve.tolist()
            d["between_curve"] = self.between_curve.tolist()
            d["within_curve"] = self.within_curve.tolist()
        return d

    def to_json(self, *, include_curves: bool = False, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(include_curves=include_curves), indent=indent)


# --------------------------------------------------------------------------- #
# Binding the math to a Comparison.                                           #
# --------------------------------------------------------------------------- #


def _evidence_for(record, endpoint: str):
    """The external C-index for ``endpoint`` from a record's predictive_performance."""
    want = f"{endpoint}_C_index_external"
    best = None
    for pp in record.predictive_performance:
        if pp.metric == want:
            return pp.value
        if "C_index" in pp.metric:  # any C-index as a fallback
            best = pp.value
    return best


def _target_samples(samples, target: str, endpoint: str) -> np.ndarray:
    if target in ("median_os_weeks", "median_pfs_weeks"):
        ep = "OS" if target == "median_os_weeks" else "PFS"
        return samples.median.get(ep, np.full(samples.n, np.nan))
    if target in samples.metrics:
        return samples.metrics[target]
    raise KeyError(
        f"unknown target {target!r}; choose median_os_weeks, median_pfs_weeks, "
        f"or a metric key {sorted(samples.metrics)}"
    )


def model_average(
    comparison,
    *,
    target: str = "median_os_weeks",
    endpoint: str = "OS",
    weights: str = "equal",
    n: int = 200,
    seed: int = 0,
) -> ModelAverage:
    """Combine ``comparison.included`` into a model-averaged forecast with the
    law-of-total-variance decomposition. See :class:`ModelAverage`."""
    included = comparison.included
    warnings: list[str] = []
    if not included:
        raise ValueError("no eligible (in-context) models to average")

    # Per-model ensembles (Axis-1 noise) reusing the shared sampling core.
    per = [
        ensemble_samples(
            comparison.ds, tr.record_id, context=comparison.context,
            drug_effect=comparison.drug_effect, exposure=comparison.exposure,
            exposure_response=comparison.exposure_response, t=comparison.t,
            survival_link=comparison.survival_link, n=n, seed=seed,
        )
        for tr in included
    ]

    ids = [tr.record_id for tr in included]
    tiers = [tr.tier for tr in included]
    evidence = [_evidence_for(comparison.ds[tr.record_id], endpoint) for tr in included]

    # Scalar target: per-model mean / variance over the finite draws.
    means, variances, scalar_idx = [], [], []
    for i, s in enumerate(per):
        arr = _target_samples(s, target, endpoint)
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            means.append(np.nan)
            variances.append(np.nan)
        else:
            means.append(float(finite.mean()))
            variances.append(float(finite.var()))
            scalar_idx.append(i)
    if not scalar_idx:
        raise ValueError(f"no included model produced a finite {target!r} over the horizon")
    if len(scalar_idx) < len(per):
        dropped = [ids[i] for i in range(len(per)) if i not in scalar_idx]
        warnings.append(
            f"{len(dropped)} model(s) produced no finite {target} over the horizon and were "
            f"dropped from the point estimate: {dropped}"
        )

    sub_means = np.array([means[i] for i in scalar_idx])
    sub_vars = np.array([variances[i] for i in scalar_idx])
    sub_tiers = [tiers[i] for i in scalar_idx]
    sub_evidence = [evidence[i] for i in scalar_idx]
    sub_ids = [ids[i] for i in scalar_idx]

    # Curves (pointwise S_bar) over the requested endpoint, across the models that
    # have it. The curve set can differ from the scalar set (e.g. a non-crossing
    # median still contributes a full curve).
    curve_idx = [i for i, s in enumerate(per) if endpoint in s.survival]
    s_bar = between_curve = within_curve = None
    if curve_idx:
        curve_means = np.array([per[i].survival[endpoint].mean(axis=0) for i in curve_idx])
        curve_withins = np.array([per[i].survival[endpoint].var(axis=0) for i in curve_idx])
        c_tiers = [tiers[i] for i in curve_idx]
        c_evidence = [evidence[i] for i in curve_idx]
        cw, _ = compute_weights(weights, c_tiers, c_evidence)
        s_bar, between_curve, within_curve = average_curve(curve_means, cw, curve_withins)
    else:
        t = comparison.t
        s_bar = np.full(t.size, np.nan)
        between_curve = np.zeros(t.size)
        within_curve = np.zeros(t.size)

    # Headline scheme decomposition + cross-scheme sensitivity.
    scheme_points: dict = {}
    chosen_w = None
    for scheme in SCHEMES:
        w, notes = compute_weights(scheme, sub_tiers, sub_evidence)
        point, *_ = decompose(sub_means, sub_vars, w)
        scheme_points[scheme] = point
        if scheme == weights:
            chosen_w = w
            warnings.extend(notes)
    if chosen_w is None:  # weights not in SCHEMES
        chosen_w, notes = compute_weights(weights, sub_tiers, sub_evidence)
        warnings.extend(notes)
        scheme_points[weights] = decompose(sub_means, sub_vars, chosen_w)[0]

    point, within_var, between_var, frac = decompose(sub_means, sub_vars, chosen_w)

    pts = list(scheme_points.values())
    weight_sensitivity = float(max(pts) - min(pts))
    if point and weight_sensitivity > 0.1 * abs(point):
        warnings.append(
            "central estimate depends materially on the weighting choice "
            f"(swing {weight_sensitivity:.1f} across schemes vs point {point:.1f})"
        )

    averaged_tier = worst_tier(tiers)
    if len(included) == 1:
        warnings.append(
            "single_eligible_model: only one in-context model; model_selection_fraction=0 "
            "is an absence of cross-checks, not a clean bill of health"
        )

    return ModelAverage(
        target=target,
        endpoint=endpoint,
        scheme=weights,
        n=n,
        t=comparison.t,
        point=point,
        curve=s_bar,
        within_var=within_var,
        between_var=between_var,
        model_selection_fraction=frac,
        weights=dict(zip(sub_ids, (float(x) for x in chosen_w))),
        per_model_mean=dict(zip(sub_ids, (float(x) for x in sub_means))),
        tier=averaged_tier,
        weight_sensitivity=weight_sensitivity,
        between_curve=between_curve,
        within_curve=within_curve,
        scheme_points=scheme_points,
        warnings=warnings,
    )


def uncertainty_decomposition(comparison, *, target: str = "median_os_weeks",
                              endpoint: str = "OS", n: int = 200, seed: int = 0) -> dict:
    """Per-scheme decomposition table for a context: each weighting scheme's
    point estimate, within/between variance, and model-selection fraction."""
    table = {}
    for scheme in SCHEMES:
        ma = model_average(
            comparison, target=target, endpoint=endpoint, weights=scheme, n=n, seed=seed,
        )
        table[scheme] = {
            "point": ma.point,
            "within_var": ma.within_var,
            "between_var": ma.between_var,
            "model_selection_fraction": ma.model_selection_fraction,
            "weights": ma.weights,
        }
    return table
