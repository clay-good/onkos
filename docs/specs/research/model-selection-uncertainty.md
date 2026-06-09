# Onkos — research spec: model-selection uncertainty & model averaging

**Status:** implemented in v0.21.0 (`onkos.combine`; steps 1–4 of §10). This
remains the design-of-record; the methodological source anchors of §8 are
documented but their Crossref-verified citation curation is still pending, honest
by design. Written in the v0.1 house style; every parameter and weight is
illustrative and `unverified` by design — the infrastructure is the contribution.

**The third uncertainty axis.** Make the disagreement between published TGI models —
the project's named load-bearing risk — a *quantified variance component* and a
*principled, honestly-weighted central forecast*, rather than the descriptive spread
the divergence view reports today. Decompose total predictive uncertainty into
within-model (parameter) and between-model (model-selection) variance via the law of
total variance; combine the eligible models into a model-averaged trajectory whose
weights are declared and whose disagreement travels with it. Never let averaging
manufacture confidence, raise a tier, or cross the population→individual line.

> The virtual-trial divergence view already greys out out-of-context models and shows
> how far the survivors disagree. This spec is the inferential completion of that view:
> from *"the models disagree by this much"* to *"here is the central tendency of the
> published models, here is how much of the remaining uncertainty is irreducible
> model-choice risk versus estimable parameter noise, and here is how fragile that
> answer is to how we weighted them."* That last clause is the honesty guardrail — the
> average is never reported without the disagreement that qualifies it.

---

## 1. The problem this extends

Onkos quantifies two of the three uncertainties in a composed survival forecast:

| Axis | Question | Today |
| --- | --- | --- |
| **Parameter (within-model)** | How much does the forecast move under the reported IIV CV of *this* model's terms? | ✅ `onkos.uncertainty.simulate_ensemble` — lognormal Monte-Carlo bands. |
| **Transportability / tier** | Is *this* model even allowed in *this* context, and how trustworthy is it? | ✅ `transportability` + tier propagation; out-of-context → tier-D + exclusion. |
| **Model-selection (between-model)** | Of the models that *are* eligible, how much does the answer depend on *which one* you pick — and what is the defensible combined answer? | ⚠️ **Descriptive only.** `Comparison.os_divergence` = max pointwise spread; `median_os_range` = min/max. No central estimate, no variance decomposition, no weighting. |

The essay names model-and-context-selection risk as *the* load-bearing idea. The
transportability machinery handles the *context* half rigorously. The *model-selection*
half is currently surfaced but not measured: a user sees five OS curves and a spread
number, and is left to eyeball the consequence. The deepening is to give the spread an
inferential structure that is standard in forecasting and in model-informed drug
development, while staying inside Onkos's no-false-precision discipline.

**Why this is the right deepening (and the right scope).** It (1) advances the project's
own stated thesis rather than adding breadth; (2) is pure post-processing of the existing
`compare` ensemble — no new dataset subsystem, near-zero schema change; (3) has direct
regulatory-science precedent (model averaging is established in dose-finding, §8); (4) is
*safe by construction* — it produces a trial-level central tendency and its disagreement,
never an individual prediction or a therapy ranking; and (5) sharpens, rather than
softens, the honest message, because the headline output is a *fraction of uncertainty
that no amount of better estimation can remove*.

---

## 2. The statistical framework

Fix a context (tumor type, line), a drug-effect (or exposure profile), and a scalar or
functional target `Q`. The eligible set is exactly what `compare()` already returns in
`Comparison.included`: models `m = 1..M` that passed the transportability check. Out-of-
context models stay in `excluded` and are **not** averaged — averaging never rehabilitates
a transported model.

Let `w_m ≥ 0`, `Σ_m w_m = 1` be model weights (§3). Treat the composed prediction as a
**mixture over models**, where within each model the parameters vary per their stored IIV.

**Model-averaged point estimate** (target `Q` = median OS in weeks, week-8 metric, `S(t)`
at a landmark time, …):

```
Q̄  =  Σ_m  w_m · E[Q | model m]          # E over within-model parameter IIV
```

**Law of total variance — the decomposition that is the whole point:**

```
Var(Q)  =  Σ_m w_m · Var(Q | m)            +   Σ_m w_m · ( E[Q | m] − Q̄ )²
        =  WITHIN  (parameter, Axis 1)      +   BETWEEN (model-selection, Axis 3)
```

and the single new headline number:

```
model_selection_fraction  =  BETWEEN / (WITHIN + BETWEEN)   ∈ [0, 1]
```

This fraction answers the question a go/no-go committee actually has: *of everything I am
uncertain about in this forecast, how much would shrink if I ran a bigger trial and nailed
the parameters (within), versus how much is structural disagreement between equally-
published models that more data on any one of them will not resolve (between)?* A high
between-fraction is precisely the signal that sends programs into doomed phase-3 trials,
and it has never had a number.

**Functional target (the survival curve itself).** Apply the same algebra pointwise. The
model-averaged survival function

```
S̄(t)  =  Σ_m w_m · S_m(t)
```

is a valid survival function (a convex combination of monotone-decreasing curves in
`[0,1]` is monotone-decreasing in `[0,1]`), and the pointwise between-model variance is
`Σ_m w_m (S_m(t) − S̄(t))²`. This yields a model-averaged OS/PFS curve with a between-model
band — the natural successor to the current overlay-of-curves figure.

**Composition with Axis 1.** `WITHIN` is obtained for free by running the existing
`simulate_ensemble` per included model and reusing its sample variance of `Q`; `BETWEEN`
needs only the per-model means `E[Q|m]`. So the integrated run is: ensemble-per-model →
combine by the identity above. No new numerical machinery, only book-keeping.

---

## 3. Weighting schemes — and an honesty boundary we will not cross

The weights decide how much each published model speaks. Onkos ships a small, **declared**
registry of schemes and *always reports the headline target under more than one*, because
weight choice is itself an uncertainty (§4).

| Scheme | `w_m ∝` | When it is the honest default | Caveat surfaced |
| --- | --- | --- | --- |
| `equal` | `1` | No basis to prefer one published model over another. **The default.** | None; maximally agnostic. |
| `tier` | declared tier score (e.g. A:B:C = 4:2:1) | A better-validated model should speak louder, but only by a *declared*, not fitted, factor. | Scores are a convention, not a probability; printed alongside. |
| `evidence` | `max(0, C_index_external − 0.5)` from `predictive_performance` | An external-validation metric exists for the eligible models. The most defensible data-driven scheme. | Only available where external validation was recorded; models without it fall back or get zero data-weight (flagged). |

**What we deliberately do NOT do — and why it matters.** Classical Bayesian model
averaging assigns `w_m = P(model m | data)`, and stacking/pseudo-BMA optimize predictive
weights — *both require the candidate models to have been fit to a common dataset*. Onkos
models are fit to **different trials, drugs, and tumor types**; there is no shared
likelihood, so a posterior model probability is **not identifiable** and would be a
fabricated quantity. Onkos therefore frames its weights as **forecast-combination weights**
(in the Bates–Granger sense), explicitly *not* posterior model probabilities, and prints
that distinction wherever weights appear. This is the same refusal that makes the rest of
the project trustworthy: we do not invent a number the data cannot support.

---

## 4. Weight-scheme sensitivity — the meta-uncertainty

A model average is only as honest as it is robust to its weights. The combiner computes the
headline target (e.g. model-averaged median OS) under **all** applicable schemes and reports
the swing:

```
weight_sensitivity  =  max_scheme Q̄(scheme)  −  min_scheme Q̄(scheme)
```

If `weight_sensitivity` is large relative to the median, the model-averaged point estimate
is fragile and is reported with an explicit *"central estimate depends materially on the
weighting choice"* warning — a fourth, smaller honesty signal nested inside the third axis.
A converged answer across schemes is itself informative: it means the combination is robust
and the divergence, while present, is weight-insensitive.

---

## 5. Tier & guardrail propagation (worst-input-wins still governs)

Model averaging **cannot raise confidence**. The rules:

- **Averaged tier = worst tier among the included models.** You cannot average a C-tier and
  a B-tier model into a B-tier forecast. This preserves the project invariant and is
  enforced (a landmark test, §6) so it cannot regress.
- **Excluded models stay excluded.** The averaged set is exactly `Comparison.included`;
  transported / failure-mode-tripping models never enter the weights, under any scheme.
- **The average never ships without its disagreement.** A `ModelAverage` result that omits
  `model_selection_fraction` is a schema/test failure. The point estimate and the between-
  model variance are a single inseparable object in every export and figure (§7).
- **A near-degenerate ensemble is flagged, not hidden.** `M = 1` (only one eligible model)
  yields `model_selection_fraction = 0` *and* a `single_eligible_model` warning — a zero is
  not a clean bill of health, it is an absence of cross-checks.

---

## 6. Reference kernel & validation landmarks

No new ODE kernel — the combiner is pure post-processing over existing trajectories. It
gets its own **landmark suite** in the spirit of `tests/test_landmarks.py`: closed-form
properties of the *combination math itself*, so the implementation is provably the
estimator it claims to be, not merely self-consistent.

| Landmark | Closed-form condition |
| --- | --- |
| **Degenerate set** | `M = 1` → `S̄ ≡ S_1`, `BETWEEN = 0`, `model_selection_fraction = 0`, `single_eligible_model` warning present. |
| **Equal-weight identity** | `Q̄` equals the arithmetic mean of `{E[Q\|m]}`. |
| **Identical components** | `M` copies of one model → `BETWEEN = 0` for *any* weights. |
| **Law of total variance** | `Total = WITHIN + BETWEEN` to ≤1e-9 on a constructed mixture of lognormals with known moments. |
| **Weight normalization & inertness** | weights sum to 1; a zero-weight model leaves `Q̄`, `S̄`, and both variances unchanged. |
| **Convex-hull bound** | `Q̄ ∈ [min_m E[Q\|m], max_m E[Q\|m]]` and `S̄(t) ∈ [min_m S_m(t), max_m S_m(t)]` — the average never extrapolates beyond the model set. |
| **Survival-function validity** | `S̄(t)` monotone non-increasing, `S̄(0)=1`, `S̄ ∈ [0,1]`. |
| **Monotone re-weighting** | raising one model's weight moves `Q̄` monotonically toward `E[Q\|that model]`. |
| **Tier floor** | averaged tier equals the worst included tier (string-min over A<B<C<D). |

This mirrors the existing two-axis validation strategy: round-trip proves exports match the
kernel; landmarks prove the kernel *is* the model it names. Here the "kernel" is the
combination estimator, and the landmarks prove it *is* the law of total variance and a
convex forecast combination — not an unconstrained curve fit.

---

## 7. API, CLI, and export surface

**Python.** A new `onkos.combine` module and a `ModelAverage` dataclass, reached through the
existing `Comparison`:

```python
cmp = onkos.compare(ds, purpose="tgi", context=dict(tumor_type="NSCLC", line="first"),
                    drug_effect=1.0)

ma = cmp.model_average(target="median_os_weeks", endpoint="OS",
                       weights="equal", n=200)   # n -> within-model ensemble depth

ma.point                       # model-averaged median OS (weeks)
ma.curve                       # S̄(t): the averaged OS survival function
ma.within_var, ma.between_var  # the law-of-total-variance components
ma.model_selection_fraction    # BETWEEN / TOTAL  — the headline number
ma.weights                     # {record_id: w_m}, and the scheme name
ma.weight_sensitivity          # swing of `point` across applicable schemes
ma.tier                        # worst included tier (cannot be raised)
ma.warnings                    # single_eligible_model / weight-fragility / etc.

dec = cmp.uncertainty_decomposition(target="median_os_weeks")   # per-scheme table
```

`ModelAverage.to_dict()/.to_json()` carry the universal `onkos:clinicalUse =
"PROHIBITED …"` flag and `NOT_FOR_CLINICAL_USE: true`, exactly as `Comparison.to_dict()`
does today.

**CLI.**

```bash
onkos compare --average --weights equal     --decompose      # adds the Axis-3 block
onkos compare --average --weights evidence   --json          # machine-readable result
```

**Exports.** The virtual-trial JSON gains an optional `model_average` block (averaged
curve, `within_var`/`between_var`, `model_selection_fraction`, weights + scheme,
`weight_sensitivity`, the worst-tier, warnings, and the clinical-use prohibition). JSON-LD /
RDF gains one predicate, `onkos:modelSelectionUncertainty`, carrying the fraction and the
scheme so the disagreement is a resolvable triple alongside `onkos:confidenceTier`. No new
NONMEM/SBML/PharmML surface — model averaging is an analysis over models, not a model.

**Schema impact: essentially none.** Weights derive from existing fields (`tier`,
`predictive_performance`); scheme definitions live in code constants, not in the dataset;
the eligible set is computed by `compare`. The only optional dataset-side addition is a
declared per-context **model-set note** (a short rationale string on which models *should*
be considered comparable for a context) — and even that is deferrable; the computed
eligible set is the default.

---

## 8. Source anchors (methodological; DOIs added at curation time)

These are well-established methods, not Onkos parameters; they anchor the framework and are
added to `dataset/citations/` through the normal Crossref/PubMed-verified pipeline, with
`review_status` honest until a human confirms each.

- **Forecast combination.** Bates, J.M. & Granger, C.W.J. (1969), *The Combination of
  Forecasts*, Operational Research Quarterly — the origin of weighted forecast averaging
  and the reason Onkos calls its weights *combination* weights, not model posteriors.
- **Bayesian model averaging (and why it does not directly apply here).** Hoeting, Madigan,
  Raftery & Volinsky (1999), *Bayesian Model Averaging: A Tutorial*, Statistical Science —
  the canonical BMA reference; cited to mark the boundary (common-data requirement) Onkos
  cannot meet across trials.
- **Stacking / pseudo-BMA caveats.** Yao, Vehtari, Simpson & Gelman (2018), *Using Stacking
  to Average Bayesian Predictive Distributions*, Bayesian Analysis — why predictive-weight
  optimization needs a shared predictive task.
- **Model averaging in drug development (the regulatory precedent).** Bretz, Pinheiro &
  Branson (2005), *Combining Multiple Comparisons and Modeling Techniques in Dose-Response
  Studies* (MCP-Mod), Biometrics — model-averaged dose-response is established, FDA/EMA-
  qualified methodology; Onkos's TGI-survival combiner is the same idea one layer up.
- **Model averaging in pharmacometrics.** Buatois, Ueckert, Frey, Retout & Mentré (2018),
  *Comparison of Model Averaging and Model Selection in Dose-Finding Trials Analyzed by
  Nonlinear Mixed-Effect Models*, AAPS Journal — direct precedent for averaging NLME models
  of the kind Onkos curates.
- **Variance decomposition.** The law of total variance (standard); used here exactly as in
  variance-based sensitivity analysis, consistent with `onkos.sensitivity`.

---

## 9. Safety & scope (unchanged hard line)

Model averaging is the most tempting place in the whole project to overreach, because a
single combined curve *looks* like an answer. The guardrails are therefore explicit:

- **Population / trial level only.** `S̄(t)` is the central tendency of *published models for
  a trial-level context*. It is **not** an estimate of any person's survival, **not** a
  prognosis, **not** a probability for an individual. No individual-level output is added.
- **No therapy ranking.** The combiner averages over *models for one context*; it never
  compares *drugs or regimens* to rank them. Cross-regimen comparison remains out of scope.
- **The average never stands alone.** It is structurally inseparable from its
  `model_selection_fraction` and worst-tier; an export or figure showing the point estimate
  without its disagreement is a test failure, not a styling choice.
- **No invented confidence.** Weights are declared forecast-combination weights, not fitted
  posterior probabilities (§3); weight-scheme fragility is reported, not smoothed (§4);
  averaging cannot raise a tier (§5).
- **The line, restated.** Any feature that takes a real patient's tumor measurement and
  returns a combined survival estimate or a therapy choice **does not get built.** Model
  averaging changes none of this; it makes the trial-level disagreement legible, and stops
  there.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Decomposition core** | `onkos.combine` with the law-of-total-variance split over `Comparison.included`, reusing `simulate_ensemble` for `WITHIN`; `model_selection_fraction`; equal-weight scheme. Landmark suite (§6). | The headline fraction computes and every landmark passes; tier floor enforced. |
| **2 — Weighting & sensitivity** | `tier` and `evidence` schemes; `weight_sensitivity` and the fragility warning; the *combination-not-posterior* labeling everywhere weights appear. | All three schemes report side by side; fragility is surfaced, not hidden. |
| **3 — Curve averaging & presentation** | Pointwise `S̄(t)` with between-model band; `ModelAverage.to_dict/json` with the clinical-use flag; `onkos compare --average`; vt-JSON `model_average` block; `onkos:modelSelectionUncertainty` RDF; a model-average figure + notebook. | The divergence view ships an averaged OS/PFS curve with its decomposition and exports it. |
| **4 — Report & audit wiring** | The dataset-health report gains a per-context model-selection-uncertainty summary; high-fraction contexts are flagged as the contexts where adding a better-validated model has the most value (curation triage, as sensitivity is for parameters). | `onkos report` ranks contexts by irreducible model-choice risk. |

Step 1 alone is a self-contained, citable contribution: an open, validated tool that splits
a composed oncology survival forecast into estimable parameter noise versus irreducible
model-selection risk is something the field reasons about qualitatively and nobody ships
openly as a tested artifact.
