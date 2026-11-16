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
from .tiers import worst_tier

__all__ = [
    "required_events",
    "horizon_hazard_ratio",
    "discriminating_events",
    "model_discriminability",
    "ModelDiscriminability",
    "logrank_statistic",
    "sample_event_times",
    "simulated_power",
    "simulated_required_events",
    "simulated_discriminability",
    "SimulatedDiscriminability",
    "feasibility_class",
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
    if not (0.0 < power < 1.0 and 0.0 < alpha < 1.0):
        raise ValueError(f"power and alpha must lie in (0, 1) (got power={power}, alpha={alpha})")
    hr = float(hazard_ratio)
    if not np.isfinite(hr) or hr <= 0.0 or hr == 1.0:
        return float("inf")
    z = _ppf(1.0 - alpha / 2.0) + _ppf(power)
    return 4.0 * z * z / (np.log(hr) ** 2)


def horizon_hazard_ratio(curve_a, curve_b, *, floor: float = 1e-6) -> float:
    """The follow-up-horizon hazard ratio between two survival curves: the ratio of their
    cumulative hazards ``H = -ln S`` at the trial horizon. Exact under proportional hazards;
    a horizon-average otherwise. ``inf`` if the reference curve has no cumulative hazard
    (flat at 1).

    The horizon is the last time at which *both* curves are still above ``floor``: once both
    arms are effectively all dead a trial observes nothing more, and clamping two vanished
    curves to the floor would read ``HR = 1`` for any pair (declaring clearly different
    models indistinguishable). Under proportional hazards the ratio is the same at every
    time, so this only removes the artifact."""
    a = np.asarray(curve_a, dtype=float)
    b = np.asarray(curve_b, dtype=float)
    alive = np.flatnonzero((a > floor) & (b > floor))
    k = int(alive[-1]) if alive.size else 0
    sa = float(np.clip(a[k], floor, 1.0))
    sb = float(np.clip(b[k], floor, 1.0))
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
        from ._serialize import to_strict_json

        return to_strict_json(self.to_dict(), indent=indent)


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
    tier = worst_tier([tr.tier for tr in cmp.included]) if cmp.included else "C"
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


# --------------------------------------------------------------------------- #
# Simulation-based logrank power — the check on the closed-form Schoenfeld core. #
# Schoenfeld is exact only asymptotically and under proportional hazards; this   #
# simulates the trial itself (event times drawn from the two population curves,  #
# uniform accrual, administrative censoring at the horizon) and runs the logrank #
# test, so it holds for crossing / non-proportional curves and small trials too. #
# --------------------------------------------------------------------------- #


def logrank_statistic(time, event, group) -> float:
    """The two-sample logrank Z statistic (group ``True`` vs ``False``), with the standard
    hypergeometric variance and exact handling of tied event times. ``0`` when no events
    carry information (e.g. no events at all)."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=bool)
    group = np.asarray(group, dtype=bool)
    if not (time.shape == event.shape == group.shape) or time.ndim != 1:
        raise ValueError("time, event and group must be aligned 1-D arrays")
    if np.any(np.isnan(time)):
        raise ValueError("time must not contain NaN")
    all_sorted = np.sort(time)
    g_sorted = np.sort(time[group])
    tj, dj = np.unique(time[event], return_counts=True)
    if tj.size == 0:
        return 0.0
    t1, d1_counts = np.unique(time[event & group], return_counts=True)
    d1 = np.zeros(tj.size)
    d1[np.searchsorted(tj, t1)] = d1_counts
    nj = all_sorted.size - np.searchsorted(all_sorted, tj, "left")  # at risk, both arms
    n1 = g_sorted.size - np.searchsorted(g_sorted, tj, "left")  # at risk, group arm
    frac = n1 / nj
    expected = float(np.sum(dj * frac))
    var_terms = dj * frac * (1.0 - frac) * (nj - dj) / np.maximum(nj - 1, 1)
    variance = float(np.sum(np.where(nj > 1, var_terms, 0.0)))
    if variance <= 0.0:
        return 0.0
    return (float(d1.sum()) - expected) / np.sqrt(variance)


def sample_event_times(t, survival, n: int, rng) -> np.ndarray:
    """Draw ``n`` event times from a population survival curve by inverse transform
    (linear interpolation between grid points). A draw the curve never reaches within the
    grid (``S(t_end) > u``) is ``inf``: no event by the horizon."""
    t = np.asarray(t, dtype=float)
    s = np.minimum.accumulate(np.clip(np.asarray(survival, dtype=float), 0.0, 1.0))
    if s.shape != t.shape:
        raise ValueError("survival must be aligned to t")
    u = rng.uniform(size=int(n))
    idx = np.searchsorted(-s, -u, "left")  # first grid index with S <= u
    out = np.full(u.size, np.inf)
    hit = idx < t.size
    i1 = idx[hit]
    i0 = np.maximum(i1 - 1, 0)
    s0, s1 = s[i0], s[i1]
    drop = s0 - s1
    w = np.where(drop > 0, (s0 - u[hit]) / np.where(drop > 0, drop, 1.0), 0.0)
    out[hit] = t[i0] + w * (t[i1] - t[i0])
    return out


def _check_power_inputs(power: float | None, alpha: float, reps: int, accrual: float, t) -> None:
    if power is not None and not 0.0 < power < 1.0:
        raise ValueError(f"power={power} must lie in (0, 1)")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha={alpha} must lie in (0, 1)")
    if reps < 10:
        raise ValueError(f"reps={reps}: need at least 10 simulated trials")
    if not 0.0 <= accrual < float(t[-1]):
        raise ValueError(f"accrual={accrual} must lie in [0, horizon)")


def simulated_power(
    curve_a, curve_b, t, *, n_per_arm: int, alpha: float = 0.05, reps: int = 400,
    accrual: float = 0.0, seed: int = 0,
) -> dict:
    """Monte-Carlo power of a two-sided logrank test to distinguish two population survival
    curves with ``n_per_arm`` patients per arm. Patients enter uniformly over ``[0, accrual]``
    weeks and are censored at the horizon ``t[-1]``. Returns the power, its Monte-Carlo
    standard error, and the mean number of observed events."""
    t = np.asarray(t, dtype=float)
    _check_power_inputs(None, alpha, reps, accrual, t)
    if n_per_arm < 1:
        raise ValueError(f"n_per_arm={n_per_arm} must be >= 1")
    n = int(n_per_arm)
    rng = np.random.default_rng(seed)
    z_crit = _ppf(1.0 - alpha / 2.0)
    group = np.r_[np.ones(n, dtype=bool), np.zeros(n, dtype=bool)]
    hits = 0
    events = 0
    for _ in range(int(reps)):
        times = np.concatenate([sample_event_times(t, curve_a, n, rng),
                                sample_event_times(t, curve_b, n, rng)])
        censor = t[-1] - rng.uniform(0.0, accrual, 2 * n)
        observed = times <= censor
        events += int(observed.sum())
        z = logrank_statistic(np.where(observed, times, censor), observed, group)
        hits += abs(z) > z_crit
    pw = hits / reps
    return {
        "n_per_arm": n,
        "power": pw,
        "mc_se": float(np.sqrt(max(pw * (1.0 - pw), 1e-12) / reps)),
        "mean_events": events / reps,
        "reps": int(reps),
    }


def simulated_required_events(
    curve_a, curve_b, t, *, power: float = 0.8, alpha: float = 0.05, reps: int = 400,
    accrual: float = 0.0, seed: int = 0, max_events: float = _INFEASIBLE_EVENTS,
) -> dict:
    """The smallest trial (patients per arm, and the events it observes) whose simulated
    logrank power reaches ``power`` — the simulation counterpart of :func:`required_events`.
    Searched by bisection on patients per arm with common random numbers. ``events`` is
    ``inf`` when even a trial observing about ``max_events`` events falls short."""
    t = np.asarray(t, dtype=float)
    _check_power_inputs(power, alpha, reps, accrual, t)
    a = np.asarray(curve_a, dtype=float)
    b = np.asarray(curve_b, dtype=float)
    # Event probability per patient with the shortest follow-up bounds the trial size.
    k = int(np.searchsorted(t, t[-1] - accrual, "right")) - 1
    per_patient = max(2.0 - a[k] - b[k], 1e-9) / 2.0

    def run(n):
        return simulated_power(a, b, t, n_per_arm=n, alpha=alpha, reps=reps,
                               accrual=accrual, seed=seed)

    n_hi = max(2, int(np.ceil(max_events / (2.0 * per_patient))))
    top = run(n_hi)
    if top["power"] < power:
        return {"events": float("inf"), "n_per_arm": None, "power": top["power"],
                "feasible": False}
    n_lo, best = 1, top
    while n_hi - n_lo > 1:
        mid = (n_lo + n_hi) // 2
        r = run(mid)
        if r["power"] >= power:
            n_hi, best = mid, r
        else:
            n_lo = mid
    return {"events": best["mean_events"], "n_per_arm": best["n_per_arm"],
            "power": best["power"], "feasible": True}


def feasibility_class(events: float) -> str:
    """``feasible`` (< 500 events), ``borderline`` (500-3000) or ``infeasible`` (> 3000):
    the declared, illustrative phase-3-scale bounds used throughout this module."""
    if events < _FEASIBLE_EVENTS:
        return "feasible"
    if events <= _INFEASIBLE_EVENTS:
        return "borderline"
    return "infeasible"


@dataclass
class SimulatedDiscriminability:
    """Closed-form (Schoenfeld) vs simulated logrank discriminability for a context."""

    context: dict
    survival_link: str | None
    joint_alpha: float | None
    power: float
    alpha: float
    reps: int
    accrual: float
    models: list
    pairs: list = field(default_factory=list)
    tier: str = "C"
    clinical_use: str = CLINICAL_USE

    @property
    def misclassified_pairs(self) -> list:
        """Pairs the closed form puts in a different feasibility class than the simulated
        trial does — where trusting Schoenfeld would change the conclusion."""
        return [p for p in self.pairs if p["schoenfeld_misclassifies"]]

    def to_dict(self) -> dict:
        def num(x):
            return round(float(x), 3) if x is not None and np.isfinite(x) else None

        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "survival_link": self.survival_link,
            "joint_alpha": self.joint_alpha,
            "power": self.power,
            "alpha": self.alpha,
            "reps": self.reps,
            "accrual_weeks": self.accrual,
            "tier": self.tier,
            "n_models": len(self.models),
            "n_misclassified": len(self.misclassified_pairs),
            "pairs": [{k: (num(v) if isinstance(v, float) else v) for k, v in p.items()}
                      for p in self.pairs],
        }

    def to_json(self, *, indent: int = 2) -> str:
        from ._serialize import to_strict_json

        return to_strict_json(self.to_dict(), indent=indent)


def simulated_discriminability(
    ds: Dataset,
    *,
    context: dict,
    survival_link: str | None = None,
    joint_alpha: float | None = None,
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
    power: float = 0.8,
    alpha: float = 0.05,
    reps: int = 400,
    accrual: float = 0.0,
    seed: int = 0,
) -> SimulatedDiscriminability:
    """For every pair of a context's eligible TGI models: the Schoenfeld events (closed form,
    from the horizon hazard ratio) beside the events a simulated logrank trial actually needs.

    ``joint_alpha`` uses the current-value joint link (:mod:`onkos.joint`) instead of a
    two-stage link — the non-proportional case where the horizon hazard ratio is only an
    average and the closed form can misjudge a pair. Pairs needing more than the infeasible
    bound are reported as infeasible rather than simulated at enormous size."""
    if t is None:
        t = np.linspace(0.0, 312.0, 625)
    t = np.asarray(t, dtype=float)
    _check_power_inputs(power, alpha, reps, accrual, t)
    if joint_alpha is not None and survival_link is not None:
        raise ValueError("pass either survival_link or joint_alpha, not both")
    if joint_alpha is not None and not np.isfinite(joint_alpha):
        raise ValueError(f"joint_alpha={joint_alpha} must be finite")

    cmp = compare(ds, purpose="tgi", context=context, drug_effect=drug_effect, t=t)
    models = [tr.record_id for tr in cmp.included]
    tier = worst_tier([tr.tier for tr in cmp.included]) if cmp.included else "C"
    if joint_alpha is not None:
        from .joint import joint_survival

        curves = {m: joint_survival(ds, m, context=context, drug_effect=drug_effect, t=t,
                                    alpha=joint_alpha).os_curve for m in models}
    else:
        curves = {m: simulate(ds, m, context=context, drug_effect=drug_effect, t=t,
                              survival_link=survival_link).os_curve for m in models}

    out = SimulatedDiscriminability(
        context=context, survival_link=survival_link,
        joint_alpha=None if joint_alpha is None else float(joint_alpha),
        power=power, alpha=alpha, reps=int(reps), accrual=float(accrual), models=models,
        tier=tier,
    )
    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            a, b = models[i], models[j]
            ca, cb = curves[a], curves[b]
            if ca is None or cb is None:
                continue
            hr = horizon_hazard_ratio(ca, cb)
            d_closed = required_events(hr, power=power, alpha=alpha)
            sim = simulated_required_events(ca, cb, t, power=power, alpha=alpha, reps=reps,
                                            accrual=accrual, seed=seed)
            d_sim = sim["events"]
            closed_cls, sim_cls = feasibility_class(d_closed), feasibility_class(d_sim)
            ratio = d_closed / d_sim if np.isfinite(d_closed) and np.isfinite(d_sim) else None
            out.pairs.append({
                "record_a": a,
                "record_b": b,
                "hazard_ratio": float(hr),
                "schoenfeld_events": float(d_closed),
                "simulated_events": float(d_sim),
                "events_ratio": None if ratio is None else float(ratio),
                "simulated_n_per_arm": sim["n_per_arm"],
                "schoenfeld_class": closed_cls,
                "simulated_class": sim_cls,
                "schoenfeld_misclassifies": closed_cls != sim_cls,
            })
    return out

