# Onkos — research spec: the integrated tumor burden — a third TGI→OS bridge metric

**Status:** implemented in v0.33.0 (`log_burden_auc` metric + `survival_link.nsclc_os_burden_auc`
record). This is the design-of-record; written in the v0.1 house style. Every value is illustrative
and `unverified` by design — the infrastructure is the contribution.

**The two-stage TGI→OS surrogate has a free choice of *which* on-treatment number drives the
hazard, and that choice changes the answer.** v0.25 made that choice an explicit, swappable field
(`structure.link_metric`) and shipped two values: the default **week-8 relative change** (the early-
shrinkage surrogate, blind to the tail) and the **growth-rate constant `k_g`** (the post-nadir
regrowth slope, blind to depth). It showed the metric choice *inverts* the model ranking. This spec
adds the natural third option — the **integrated tumor burden**, the time-averaged log relative
tumor size over the horizon (the AUC of the log-size curve) — which is the only one of the three
that is sensitive to **both** depth of response and the regrowth tail. The finding: a "comprehensive"
integrated metric does **not** dissolve the metric-choice axis. It produces a *third, distinct*
ranking, and in doing so it exposes a specific pathology of the pure-tail metric.

> Why integrate at all? The week-8 change reads one early timepoint; `k_g` reads one terminal slope.
> Both throw away most of the trajectory. The integrated burden keeps it: it is the log geometric-
> mean relative tumor size a patient population carries across the observation window. A tumor that
> shrinks deeply and stays down accrues little burden; one that shrinks deeply then regrows fast, or
> one that never shrinks, accrues a lot. It is the single covariate that "sees" the whole curve —
> and it is a standard summary, since TGI models live in log space (sizes are log-normal).

---

## 1. The problem this extends

A survival link consumes one TGI metric as its hazard covariate `x` (spec §6; `simulate.py`):

```
S(t) = exp(-(t/scale)^shape * exp(beta * x)),   x = metrics[link.structure.link_metric]
```

v0.25 made `link_metric` a declared, swappable field. The two shipped covariates summarize opposite
ends of the trajectory:

| metric | what it reads | what it is blind to |
| --- | --- | --- |
| `week8_relative_change` (default) | the tumor size at one early landmark (week 8) | the entire post-week-8 tail (regrowth, durability) |
| `tumor_growth_rate_kg` | the log-linear slope of the late regrowth phase | the depth of response (whether the tumor ever shrank at all) |

The new covariate fills the gap between them:

```
log_burden_auc = (1/T) ∫₀ᵀ log( max(v(t)/y0, floor) ) dt        # time-averaged log relative burden
```

— the AUC of the log relative tumor-size curve, time-averaged over the horizon `T`. It reads the
*whole* trajectory: depth lowers it (a deep nadir contributes large-negative log-burden), and a
regrowth tail raises it (sustained large size accrues positive log-burden). Eradication is floored at
a detection limit (`v/y0 = 1e-3`, a complete response — clinically indistinguishable from zero) so the
integral is a stable summary, not a `−∞`-dominated one. It is a pure post-processing metric over the
*existing* kernels — exactly the v0.25 move (one metric + one record, no kernel, no schema change).

**Why this is the right deepening (and the right scope).** It (1) is the reserve list's "more bridge
metrics for the v0.25 `link_metric` field"; (2) is pure post-processing — a metric added to the
model-agnostic Stein/Bruno panel and one non-default survival-link record, with **zero** change to
the default divergence view (the new link is `default=false`, reached only via `survival_link=`); (3)
completes the depth/tail pair with the one summary that integrates both, so the metric-choice axis is
now sampled at its two extremes *and* its center; (4) is *safe by construction* — a population
trajectory summary over published model structures, no patient data; (5) reinforces the load-bearing
message: even the "use everything" metric is still a *choice*, and it disagrees with both extremes.

---

## 2. The result — a third ranking, and the pathology it exposes

For NSCLC first line at unit drug effect, the five eligible TGI models, ranked by median OS under
each of the three bridge metrics (horizon 260 wk; illustrative):

| Model | nadir / y₀ | week-8 mOS | **k_g** mOS | **burden** mOS |
| --- | --- | --- | --- | --- |
| Norton-Simon (complete responder) | 0.00 | 58 | **102** | **151** |
| Claret (phenom. resistance) | 0.03 | 91 | 39 | 81 |
| two-population (mechanistic) | 0.03 | 94 | 32 | 66 |
| acquired resistance | 0.10 | 92 | 32 | 56 |
| Wang biexponential (minimal responder) | 0.75 | 54 | **44** | **43** |

Read the three orderings (best → worst OS):

- **week-8** (depth-only, tail-blind): two-pop > acquired > Claret > **Norton-Simon** > Wang.
  The complete responder ranks **4th of 5** — its early kill is gradual, so its week-8 shrinkage is
  shallow (−0.29) and the surrogate undervalues it. The deep *early* shrinkers (two-pop, acquired) top
  the list even though they are doomed in the tail.
- **k_g** (tail-only, depth-blind): Norton-Simon > **Wang** > Claret > two-pop > acquired.
  The complete responder jumps to 1st (no regrowth). **But Wang — a tumor that barely responds (nadir
  is 75% of baseline) — ranks 2nd**, because its regrowth *slope* happens to be slow. `k_g` cannot see
  that Wang never shrank.
- **burden** (depth *and* tail): Norton-Simon > Claret > two-pop > acquired > **Wang**.
  The complete responder stays 1st, the deep-but-doomed resistance models demote to the middle, and
  **Wang drops to last** — because the integrated burden penalizes a tumor that spent the whole
  horizon large, regardless of how slowly it grew.

Two findings, both robust to the horizon (156–260 wk) and the eradication floor:

1. **The integrated metric is a third distinct total order** — it agrees with neither week-8 nor `k_g`.
   "Which bridge metric" remains a live model-selection axis even when you reach for the metric that
   uses the whole trajectory. There is no neutral choice.
2. **The integrated metric repairs a specific pathology of `k_g`.** The pure-tail metric ranks the
   *minimal responder* (Wang) second-best, because a slow regrowth slope looks favorable in isolation;
   it is blind to the fact the tumor never got small. The integrated burden, which weighs depth, ranks
   Wang last — where a clinician would. Tail-sensitivity without depth-sensitivity is its own failure
   mode, and an integrated metric is the honest fix.

**The honest framing.** Onkos does not declare a winning bridge metric — all three are illustrative,
tier-C, `unverified`, and each is defensible from the literature (early shrinkage: Wang/Claret; `k_g`:
Stein/Bruno; integrated longitudinal burden: registrational longitudinal-size modeling). It makes the
metric a *declared, swappable field* and shows that the choice carries an unquantified model-selection
risk: a two-stage surrogate that silently picks one metric inherits a ranking the other metrics
contradict.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** The burden-AUC link carries its propagated tier (C); out-of-context
  transport floors to D with a warning, exactly as for the week-8 and `k_g` links.
- **Default view untouched.** The new link is `structure.default = false`, so the default
  virtual-trial divergence view, every default export, and the v0.21–v0.32 default-link numbers are
  **byte-identical**. The metric is added to the panel (every trajectory now reports it) but no default
  consumer reads it.
- **Population / trial level only.** The metric summarizes a population trajectory; it is never an
  individual prediction and never ranks a therapy — it ranks *models under a context*.
- **A metric, not a verdict.** `log_burden_auc` is a descriptive summary of a model's simulated curve,
  not a claim about a drug. The eradication floor is a stated convention, not a measurement.

---

## 4. Validation landmarks

The metric is a closed-form summary, so the landmarks (`tests/test_burden_auc.py`) are exact:

| Landmark | Condition |
| --- | --- |
| **Baseline ⇒ zero burden** | a trajectory held at baseline (`v ≡ y0`) has `log_burden_auc = 0` — the same zero-point convention as `week8_relative_change`, so a no-effect tumor maps to the baseline hazard. |
| **Constant shrink ⇒ log of the ratio** | a trajectory held at `c·y0` (constant `c`) has `log_burden_auc = log(c)` exactly (the time-average of a constant). |
| **Monotone in burden** | for two flat trajectories, the larger constant size has the larger (less negative) burden — the covariate orders by integrated size. |
| **Eradication is floored, finite** | a trajectory reaching `v = 0` yields a finite burden (bounded below by `log(floor)`), not `−∞` — the integral is stable under complete response. |
| **Tail-sensitive where week-8 is blind** | matched on the week-8 change, a deep-then-regrow trajectory has a larger burden than a deep-and-stays-down one (the two-population vs complete-responder contrast) — the metric sees the tail. |
| **Depth-sensitive where k_g is blind** | the minimal responder (shallow nadir, slow regrowth) has a larger burden than a deep responder with the *same* `k_g` — the metric sees the depth. |
| **Re-ranks models** | the burden-AUC OS ranking of the NSCLC model set differs from both the week-8 and the `k_g` ranking (a third order). |
| **Horizon monotonicity** | extending the horizon over a regrowing tail does not decrease the burden (more time large ⇒ more burden) — the documented cumulative-summary property. |
| **Tier & guardrails** | the propagated tier rides through; out-of-context transport floors to D; the default view is unchanged. |

---

## 5. API, CLI, and surface

The burden-AUC link is an ordinary non-default survival link — every existing surface works on it via
`survival_link=`, exactly like the `k_g` link:

```python
# The integrated-burden bridge metric drives the OS hazard when selected explicitly.
b = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0,
                   survival_link="survival_link.nsclc_os_burden_auc")
b.metrics["log_burden_auc"]      # the covariate (now in every trajectory's metric panel)
b.median_os                      # OS under the integrated-burden link

# It is an eligible survival link, so it joins the model-selection-budget V_link factor for free.
from onkos.budget import eligible_survival_links
eligible_survival_links(ds, ctx, "OS")   # now 4 for NSCLC/first: week8, cox, k_g, burden_auc
```

**No new module or CLI command** — the contribution is a metric + a record + the finding, surfaced
through the existing `simulate` / `compare` / `budget` paths and a figure. The metric appears in
`extract_tgi_metrics`'s panel, so `onkos`'s existing metric-aware surfaces report it automatically.

---

## 6. Source anchors (methodological; values illustrative)

- **Longitudinal tumor-size → OS modeling.** Wang, Y. et al. (2009), *Elucidation of relationship
  between tumor size and survival in non-small-cell lung cancer patients can aid early decision making
  in clinical drug development*, Clin. Pharmacol. Ther. 86:167–174 (DOI 10.1038/clpt.2009.64) — the
  registrational longitudinal-tumor-size modeling family the integrated-burden covariate belongs to
  (the primary citation, shared with the NSCLC week-8 and Cox links).
- **Growth-rate-constant framework.** Stein, W.D. et al. (2008) (DOI 10.1634/theoncologist.2008-0075)
  — the `k_g` metric this one is contrasted against (the pure-tail extreme).
- **Tumor-dynamics / TGI-metric review.** Bruno et al. (2020) (DOI 10.1158/1078-0432.CCR-19-0287) —
  surveys the spectrum of TGI metrics linked to survival, of which an integrated/time-averaged burden
  is one.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not a winning metric.** The spec adds a third *option*, never a recommendation that the integrated
  burden is the right covariate; the point is the choice exists and disagrees.
- **Not fitted.** `beta`/`scale`/`shape` are declared, illustrative link parameters, never estimated
  from data here; the eradication floor is a stated convention.
- **Not cross-context (yet).** The introduction adds the metric (dataset-wide, since it is in the
  panel) and one NSCLC link record (mirroring v0.25's single-context introduction of the `k_g` link,
  broadened later). Adding burden-AUC links to breast/CRC/HCC/melanoma is a clean breadth follow-on
  that would make the third-ranking finding dataset-wide.
- **Not a new endpoint.** This is a bridge metric for the *existing* OS link, not a new survival
  endpoint; PFS routes (v0.30) and ORR/DoR (v0.27–v0.28) are untouched.

---

## 8. Safety & scope (unchanged hard line)

- **Population / trial level only.** A trajectory summary and population OS over published model
  structures; nothing is an individual prediction.
- **No therapy ranking.** The metric ranks *models under a context*, never treatments.
- **Cannot raise a tier.** The burden-AUC link is tier C; out-of-context transport floors to D.
- **The line, restated.** Any feature that takes a real patient's tumor measurement and returns a
  prognosis or a therapy choice **does not get built.** Adding a third way to summarize a simulated
  trajectory changes none of this.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Metric** | `log_burden_auc` added to `extract_tgi_metrics` (floored, version-agnostic trapezoid, time-averaged). | The closed-form landmarks (baseline⇒0, constant⇒log c, floored eradication) pass. |
| **2 — Record** | `survival_link.nsclc_os_burden_auc` (`default=false`, `link_metric=log_burden_auc`, calibrated `beta`/`scale`), Wang primary citation, tier C. | `onkos validate` passes; it appears in `eligible_survival_links` for NSCLC/first (now 4). |
| **3 — The third-ranking finding** | the burden-AUC OS ranking, its distinctness from week-8 and `k_g`, and the `k_g`-repairs-Wang pathology, via the existing simulate/budget machinery. | The dynamics landmarks and the §2 table hold; default view byte-identical. |
| **4 — Surfaces** | a three-metric figure + a CI-executed notebook; README section; report/figures refreshed for the new record. | The three-way metric axis is visualized and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, tested integrated-burden bridge metric beside
the early-shrinkage and growth-rate metrics. Step 3 is the payload: showing that the metric that uses
the whole trajectory produces a *third* model ranking — and that the pure-tail metric's apparent
comprehensiveness hides a depth-blind pathology an integrated metric repairs — is the quantitative core
of this spec, shipped as a tested artifact.
