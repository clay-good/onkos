"""The model-selection budget — structural variance decomposition (the capstone).

Onkos has, one axis at a time, made every structural choice in a composed survival
forecast first-class: which TGI model (:mod:`onkos.combine`), whether the parameters
are estimable (:mod:`onkos.identify`), which interaction model
(:mod:`onkos.interaction`), which resistance mechanism, and which on-treatment metric /
survival model drives the hazard (the ``link_metric`` / ``survival_link`` axis). This
module accounts for them *together*: it splits the total uncertainty of a forecast into
**parameter noise** and the **between-axis variance contributed by each structural
choice**, and names which assumption is the biggest single driver.

It is variance-based sensitivity analysis (a balanced two-way ANOVA / first-order Sobol
decomposition) applied to the *discrete structural choices* rather than to the continuous
parameters — the structural analog of :mod:`onkos.sensitivity`'s parameter tornado, and
the honest synthesis of the whole model-selection arc:

    Var(Q) = WITHIN(parameter) + V_model(TGI choice) + V_link(survival choice) + V_inter

Collapsing the survival-link factor to one level recovers exactly the v0.21 within/between
split — the budget is a strict generalization.

Population / trial level only. NOT an individual prediction, NOT a therapy ranking, NOT a
model recommendation: it attributes variance across *assumptions*. The structural share is
irreducible by more data (the honest opposite of false precision), and the budget cannot
raise a tier.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .compare import compare
from .load import Dataset
from .tiers import worst_tier
from .uncertainty import ensemble_samples

__all__ = [
    "Budget",
    "variance_components",
    "eligible_survival_links",
    "model_selection_budget",
]

_AXIS_LABELS = {
    "parameter": "parameter noise (within-model IIV)",
    "tgi_model": "TGI-model choice",
    "survival_link": "survival-link choice (metric / structure)",
    "interaction": "model x link interaction",
}


# --------------------------------------------------------------------------- #
# Pure variance-component algebra — the landmark-tested estimator core.        #
# Balanced two-way layout: rows = TGI model (factor A), cols = survival link    #
# (factor B); ``cell_means`` and ``cell_within`` are (M, L) arrays.             #
# --------------------------------------------------------------------------- #


def variance_components(cell_means, cell_within) -> dict:
    """Balanced two-way variance-component decomposition of a target ``Q``.

    ``cell_means[a, b] = E[Q | model a, link b]`` (over parameter IIV) and
    ``cell_within[a, b] = Var[Q | model a, link b]`` (parameter noise in the cell).
    Returns ``within``, ``v_model``, ``v_link``, ``v_inter``, ``between``, ``total``,
    the grand mean, and the simplex of fractions. The components are non-negative and
    sum exactly to ``total`` (the ANOVA sum-of-squares identity, balanced design)."""
    means = np.atleast_2d(np.asarray(cell_means, dtype=float))
    within_cells = np.atleast_2d(np.asarray(cell_within, dtype=float))
    grand = float(means.mean())
    within = float(within_cells.mean())

    row_means = means.mean(axis=1)  # per TGI model, averaged over links
    col_means = means.mean(axis=0)  # per survival link, averaged over models
    v_model = float(np.mean((row_means - grand) ** 2))   # main effect A
    v_link = float(np.mean((col_means - grand) ** 2))    # main effect B
    v_cells = float(np.mean((means - grand) ** 2))       # variance of all cell means
    v_inter = max(0.0, v_cells - v_model - v_link)       # residual = interaction (>= 0)

    between = v_cells
    total = within + between
    if total > 0:
        fractions = {
            "parameter": within / total,
            "tgi_model": v_model / total,
            "survival_link": v_link / total,
            "interaction": v_inter / total,
        }
    else:
        fractions = {"parameter": 1.0, "tgi_model": 0.0, "survival_link": 0.0,
                     "interaction": 0.0}
    return {
        "grand_mean": grand,
        "within": within,
        "v_model": v_model,
        "v_link": v_link,
        "v_inter": v_inter,
        "between": between,
        "total": total,
        "fractions": fractions,
    }


# --------------------------------------------------------------------------- #
# The result object.                                                          #
# --------------------------------------------------------------------------- #


@dataclass
class Budget:
    context: dict
    endpoint: str
    target: str
    n: int
    models: list  # factor-A levels (TGI model record ids)
    links: list  # factor-B levels (survival-link record ids)
    grand_mean: float
    within: float
    v_model: float
    v_link: float
    v_inter: float
    total: float
    fractions: dict  # parameter / tgi_model / survival_link / interaction -> share
    tier: str
    cell_means: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    warnings: list = field(default_factory=list)

    @property
    def between(self) -> float:
        return self.v_model + self.v_link + self.v_inter

    @property
    def structural_fraction(self) -> float:
        """Share of the total that is irreducible structural-choice risk (everything
        except the parameter/within component) — not shrinkable by a bigger trial."""
        return self.fractions["tgi_model"] + self.fractions["survival_link"] + \
            self.fractions["interaction"]

    @property
    def dominant(self) -> str:
        """The axis with the largest share — standardize / validate this first."""
        return max(self.fractions, key=self.fractions.get)

    def to_dict(self) -> dict:
        """JSON-serializable budget. Carries the clinical-use prohibition; the
        structural share and dominant axis travel with the point variances."""
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "endpoint": self.endpoint,
            "target": self.target,
            "n": self.n,
            "models": list(self.models),
            "links": list(self.links),
            "grand_mean": self.grand_mean,
            "components": {
                "within": self.within,
                "v_model": self.v_model,
                "v_link": self.v_link,
                "v_inter": self.v_inter,
            },
            "total": self.total,
            "fractions": self.fractions,
            "structural_fraction": self.structural_fraction,
            "dominant": self.dominant,
            "tier": self.tier,
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


# --------------------------------------------------------------------------- #
# Binding the algebra to a forecast.                                          #
# --------------------------------------------------------------------------- #


def eligible_survival_links(ds: Dataset, context: dict, endpoint: str = "OS") -> list:
    """Every curated survival link for ``context`` and ``endpoint`` — default *and*
    non-default (week-8 Weibull, Cox, k_g, …) — since each is an alternative survival
    model a user could pick. Sorted by id for a deterministic factor order."""
    tumor_type = context.get("tumor_type")
    line = context.get("line") or context.get("line_of_therapy")
    out = []
    for r in ds:
        if r.purpose != "survival_link":
            continue
        dc = r.derivation_context
        if (
            dc
            and dc.tumor_type == tumor_type
            and (line is None or dc.line_of_therapy == line)
            and r.structure.get("endpoint", "OS") == endpoint
        ):
            out.append(r.id)
    return sorted(out)


def _target_median(samples, endpoint: str) -> np.ndarray:
    return samples.median.get(endpoint, np.full(samples.n, np.nan))


def model_selection_budget(
    ds: Dataset,
    *,
    context: dict,
    endpoint: str = "OS",
    target: str = "median_os_weeks",
    drug_effect: float = 1.0,
    t: np.ndarray | None = None,
    n: int = 200,
    seed: int = 0,
) -> Budget:
    """Decompose the total variance of ``target`` for a context across the structural
    factors — TGI model (``compare().included``) and survival link (all eligible) —
    plus parameter noise, via the balanced two-way variance-component identity. See
    :class:`Budget`.

    Cells that produce no finite target over the horizon drop their row/column so the
    design stays balanced; the drop is recorded in ``warnings``."""
    if t is None:
        t = np.linspace(0.0, 312.0, 625)  # long horizon so every cell reaches its median
    t = np.asarray(t, dtype=float)
    if target != "median_os_weeks":
        # The grid target is the endpoint median; other targets are a clean extension.
        raise ValueError("model_selection_budget currently supports target='median_os_weeks'")

    cmp = compare(ds, purpose="tgi", context=context, drug_effect=drug_effect, t=t)
    models = [tr.record_id for tr in cmp.included]
    links = eligible_survival_links(ds, context, endpoint)
    warnings: list[str] = []
    if not models:
        raise ValueError("no eligible (in-context) TGI models for this context")
    if not links:
        raise ValueError(f"no survival links for this context and endpoint {endpoint!r}")

    means = np.full((len(models), len(links)), np.nan)
    withins = np.zeros((len(models), len(links)))
    for a, m in enumerate(models):
        for b, lk in enumerate(links):
            s = ensemble_samples(
                ds, m, context=context, drug_effect=drug_effect, t=t,
                survival_link=lk, n=n, seed=seed,
            )
            arr = _target_median(s, endpoint)
            finite = arr[np.isfinite(arr)]
            if finite.size:
                means[a, b] = float(finite.mean())
                withins[a, b] = float(finite.var())

    # Keep the design balanced: drop any model (row) or link (col) with a missing cell.
    good_rows = ~np.isnan(means).any(axis=1)
    good_cols = ~np.isnan(means).any(axis=0)
    if not good_rows.all() or not good_cols.all():
        dropped_m = [models[i] for i in range(len(models)) if not good_rows[i]]
        dropped_l = [links[j] for j in range(len(links)) if not good_cols[j]]
        if dropped_m:
            warnings.append(f"dropped {len(dropped_m)} model(s) with no finite {target}: {dropped_m}")
        if dropped_l:
            warnings.append(f"dropped {len(dropped_l)} link(s) with no finite {target}: {dropped_l}")
    models = [m for i, m in enumerate(models) if good_rows[i]]
    links = [lk for j, lk in enumerate(links) if good_cols[j]]
    means = means[np.ix_(good_rows, good_cols)]
    withins = withins[np.ix_(good_rows, good_cols)]
    if means.size == 0:
        raise ValueError(f"no complete (model x link) cell produced a finite {target}")

    comp = variance_components(means, withins)

    if len(models) == 1:
        warnings.append(
            "single_tgi_model: only one in-context TGI model; v_model=0 is an absence of "
            "cross-checks, not a clean bill of health"
        )
    if len(links) == 1:
        warnings.append(
            "single_survival_link: only one survival link for this endpoint; v_link=0 is an "
            "absence of cross-checks, not a clean bill of health"
        )

    tier = worst_tier([tr.tier for tr in cmp.included]) if cmp.included else "C"

    return Budget(
        context=context,
        endpoint=endpoint,
        target=target,
        n=n,
        models=models,
        links=links,
        grand_mean=comp["grand_mean"],
        within=comp["within"],
        v_model=comp["v_model"],
        v_link=comp["v_link"],
        v_inter=comp["v_inter"],
        total=comp["total"],
        fractions=comp["fractions"],
        tier=tier,
        cell_means=means,
        warnings=warnings,
    )
