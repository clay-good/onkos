# Onkos — research spec: the model-selection budget (structural variance decomposition)

**Status:** implemented in v0.26.0 (`onkos.budget`). This is the design-of-record; the
methodological source anchors of §7 are documented but their Crossref-verified citation
curation is still pending, honest by design. Written in the v0.1 house style; every value
is illustrative and `unverified` by design — the infrastructure is the contribution.

**The capstone.** Onkos has, one axis at a time, made every structural choice in a composed
survival forecast first-class: *which TGI model* (the divergence view + model-selection
fraction, v0.21), *whether the parameters are even estimable* (identifiability, v0.22),
*which interaction model* for combinations (v0.23), *which resistance mechanism* (v0.24),
and *which on-treatment metric / survival model* drives the hazard (v0.25). Each was scored
in isolation. This spec accounts for them **together**: a single **model-selection budget**
that splits the total uncertainty of a forecast into *parameter noise* and the *between-axis
variance contributed by each structural choice* — and names which assumption is the biggest
single driver of the forecast, i.e. where standardization or validation buys the most.

> The project's parameter sensitivity (v0.9) answers "which *parameter*'s uncertainty moves
> the prediction, so verify it first." This is the same question one level up: "which
> *structural assumption* moves the prediction, so standardize or validate it first." It is
> variance-based sensitivity analysis applied to the discrete modeling choices rather than to
> the continuous parameters — the structural analog of the tornado, and the honest synthesis
> of the whole model-selection arc.

---

## 1. The problem this extends

A composed Onkos forecast is a stack of choices: `TGI model → on-treatment metric → survival
link → OS`, each parameterized with its own inter-individual variability. Onkos already
quantifies the pieces, but only one factor at a time:

| Uncertainty source | Question | Scored by |
| --- | --- | --- |
| **Parameter (within)** | how much does the forecast move under the IIV? | `uncertainty`, `sensitivity` (v0.7, v0.9) |
| **TGI-model choice** | how much depends on which growth/kill/resistance model? | `combine` — model-selection fraction (v0.21) |
| **Survival-link choice** | how much depends on which metric / survival structure (week-8 vs k_g; Weibull vs Cox)? | `simulate` survival_link / link_metric (v0.13, v0.25) — surfaced, **not yet decomposed jointly** |

The v0.21 decomposition split variance into *within* (parameter) and *between* (TGI model)
for a **fixed** survival link. But v0.25 showed the survival-link choice can *invert* the
answer — so for a real go/no-go forecast it is a second, co-equal structural axis, and the
two have never been put on one ledger. The deepening is the **joint** decomposition: total
variance attributed across *all* the structural axes at once, so their relative weights — and
their interaction — are visible in a single budget.

**Why this is the right deepening (and the right scope).** It (1) is the synthesis the arc
was building toward — it advances the thesis by *unifying* the axes, not adding a new one;
(2) is pure post-processing over the existing `compare` + `ensemble_samples` machinery, no
new dataset subsystem and near-zero schema change (mirrors `combine`/`identify`); (3) is
textbook **two-way ANOVA variance components / variance-based sensitivity** (§7), so the math
is standard and the landmarks are closed-form; (4) is *safe by construction* — a trial-level
variance accounting, no individual prediction, no therapy ranking; and (5) sharpens the
honest message to its strongest, most decision-relevant form: *of everything uncertain in
this go/no-go forecast, here is the share that is irreducible structural-choice risk, broken
out by which choice — and a bigger trial shrinks only the parameter share.*

---

## 2. The decomposition

Fix a context (tumor type, line), an endpoint, and a scalar target `Q` (median OS in weeks).
Two **structural factors** index the forecast:

* **A — TGI model**, with levels `a = 1..M` = exactly `compare().included` (the in-context,
  transportability-passing models; out-of-context models stay excluded — averaging never
  rehabilitates a transported model);
* **B — survival link**, with levels `b = 1..L` = every curated survival link for the
  context and endpoint (default *and* non-default: e.g. week-8 Weibull, Cox, k_g), since each
  is an alternative survival model a user could pick.

For each cell `(a, b)` run the existing per-model IIV ensemble (`ensemble_samples`) and reduce
its finite `Q` draws to a conditional mean and within-cell variance:

```
μ_ab = E[Q | model a, link b]            (over parameter IIV)
σ²_ab = Var[Q | model a, link b]         (parameter noise within the cell)
```

The **balanced two-way variance-component identity** then splits the total:

```
WITHIN   = mean_ab σ²_ab                                  parameter noise (reducible by more data)
V_model  = (1/M) Σ_a ( mean_b μ_ab − μ̄ )²                 TGI-model main effect
V_link   = (1/L) Σ_b ( mean_a μ_ab − μ̄ )²                 survival-link main effect
V_inter  = Var_ab(μ_ab) − V_model − V_link                model×link interaction (≥ 0)
BETWEEN  = Var_ab(μ_ab) = V_model + V_link + V_inter       total structural-choice variance
TOTAL    = WITHIN + BETWEEN
```

where `μ̄ = mean_ab μ_ab`. Each component is non-negative and they **sum exactly to TOTAL**
(the ANOVA sum-of-squares identity, balanced design). The budget is the vector of fractions:

```
budget = ( WITHIN, V_model, V_link, V_inter ) / TOTAL          ∈ simplex
```

and the single headline is the **dominant axis** — the structural factor with the largest
share, the one where standardizing the assumption (or validating to eliminate the choice)
removes the most forecast uncertainty. A high `V_inter` is itself informative: it means the
right survival link *depends on* which TGI model you picked (or vice versa), so the axes
cannot be standardized independently.

**Composition with v0.21.** Collapsing factor B to a single level (`L = 1`) recovers exactly
the v0.21 within/between split (`V_link = V_inter = 0`, `V_model = BETWEEN`). The budget is a
strict generalization, enforced by a landmark.

---

## 3. Weighting & honesty (inherited from v0.21–v0.23)

The default is the **equal/uniform** layout (the maximally-agnostic balanced design). Where a
declared weighting is wanted, the same forecast-combination weights as `combine` apply to
factor A (tier / evidence), reported alongside — never posterior model probabilities (the
models are fit to different trials; a posterior is not identifiable, §3 of the
model-selection-uncertainty spec). The budget **cannot raise a tier** (it carries the worst
included tier), and a near-degenerate design is flagged, not hidden: `M = 1` ⇒ `V_model = 0`
with a `single_tgi_model` note; `L = 1` ⇒ `V_link = 0` with a `single_survival_link` note (a
zero is an absence of cross-checks, not a clean bill of health).

---

## 4. The expected result — the structural axes can dominate

For NSCLC first line, factor B is rich (week-8 Weibull, Cox, k_g), and v0.25 showed those
links disagree enough to *invert* the model ranking. So the budget is expected to put a
**large share on the survival-link axis** — frequently larger than the parameter share — i.e.
the most consequential uncertainty in the OS forecast is not the parameters a bigger trial
would pin down, nor even the tumor-growth model, but *which survival model you assumed*. That
is a decision-grade, and uncomfortable, finding: it says the standardization with the highest
leverage is the survival endpoint/metric, not the tumor model everyone argues about. The
budget makes that claim quantitative instead of rhetorical.

---

## 5. Reference kernel & validation landmarks

No new kernel — the budget is pure post-processing over existing trajectories. It gets its own
**landmark suite** (in the spirit of `test_combine.py`): closed-form properties of the
variance-component algebra, proven on constructed grids with known moments.

| Landmark | Closed-form condition |
| --- | --- |
| **Sum identity** | `WITHIN + V_model + V_link + V_inter = TOTAL`; the four fractions sum to 1. |
| **Non-negativity** | every component ≥ 0 (incl. the residual interaction) on any grid. |
| **Single-factor collapse** | `L = 1` ⇒ `V_link = V_inter = 0`, `V_model = BETWEEN` (recovers the v0.21 split); symmetric for `M = 1`. |
| **Identical cells** | all `μ_ab` equal ⇒ every between component = 0; budget = (1, 0, 0, 0). |
| **Pure main effect** | `μ_ab` depends only on `b` ⇒ `V_model = V_inter = 0`, `V_link = BETWEEN`. |
| **Additive layout** | `μ_ab = r_a + c_b` ⇒ `V_inter = 0` (no interaction for an additive grid). |
| **ANOVA identity vs SS** | the component variances match a direct balanced two-way sum-of-squares decomposition (÷ N) to ≤1e-9. |
| **Convex-hull bound** | `μ̄ ∈ [min_ab μ_ab, max_ab μ_ab]`. |
| **Zero-within degeneracy** | all `σ²_ab = 0` ⇒ `WITHIN = 0`, budget is purely structural. |
| **Tier floor** | the budget's tier equals the worst included TGI-model tier (cannot be raised). |

This mirrors the project's discipline: the estimator *is* the two-way variance-component
decomposition (a balanced ANOVA / first-order Sobol over the structural factors), not an
unconstrained attribution.

---

## 6. API, CLI, and surface

```python
b = onkos.model_selection_budget(ds, context=dict(tumor_type="NSCLC", line="first"),
                                 endpoint="OS", target="median_os_weeks", n=200)

b.total                      # total variance of median OS (week²)
b.within, b.v_model, b.v_link, b.v_inter      # the four components
b.fractions                  # {"parameter":…, "tgi_model":…, "survival_link":…, "interaction":…}
b.dominant                   # the largest-share axis — standardize/validate this first
b.models, b.links            # the factor levels (the grid axes)
b.tier, b.warnings           # worst included tier; single-level / degeneracy notes
b.to_dict()                  # carries the clinical-use prohibition
```

**CLI.**

```bash
onkos budget --tumor-type NSCLC --line first --endpoint OS        # the structural budget
onkos budget --json
```

prints the four-way budget, the dominant axis, and the grid (TGI models × survival links).

**Report.** `onkos report` gains a per-context **model-selection budget** section naming each
context's dominant structural axis and the parameter-vs-structural split (binned for a
byte-stable diff), so curation triage now ranks not just *which parameter* (sensitivity) and
*which context* (model-selection fraction) but *which structural axis* most needs
standardization.

**No new export model** — the budget is an analysis over models, not a model.

---

## 7. Source anchors (methodological; DOIs added at curation time)

- **Variance-based sensitivity (Sobol).** Sobol, I.M. (2001), *Global sensitivity indices for
  nonlinear mathematical models and their Monte Carlo estimates* — the first-order index over
  factors that the budget computes for *structural* factors (Onkos already uses it for
  parameters in `sensitivity`).
- **Two-way ANOVA variance components.** The classical balanced sum-of-squares decomposition
  (`SS_total = SS_A + SS_B + SS_AB`) the budget is, expressed as variances.
- **Law of total variance.** The within/between identity (already cited for `combine`), here
  generalized to multiple between-factors.
- **Model averaging in drug development.** MCP-Mod (Bretz 2005) and NLME model averaging
  (Buatois 2018) — the regulatory precedent for treating structural model choice as a
  first-class, quantified uncertainty (already cited for `combine`).

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not a third+ factor in v0.x.** Two structural factors (TGI model, survival link) plus
  parameter within. Adding exposure-response choice or the interaction model as further ANOVA
  factors is a clean, separable extension.
- **Not posterior model probabilities.** The layout is a declared (uniform or
  combination-weighted) design, never a fitted model posterior.
- **Not optimal standardization.** The budget *identifies* the dominant axis; it does not
  prescribe a fix. Recommending which model to adopt is out of scope (it edges toward picking
  a winner, which the project refuses).
- **No individual-level variance.** Trial-level only.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only.** Every quantity is a trial-level variance of a composed
  forecast over published models; nothing is a person's prognosis or its uncertainty.
- **No therapy ranking, no model recommendation.** The budget attributes variance across
  *assumptions*; it never ranks treatments and never declares a model the winner.
- **The structural share is irreducible by data.** The headline distinction — parameter
  (shrinkable by a bigger trial) vs structural (not) — is the honest opposite of false
  precision; the budget never implies the structural share will go away with more patients.
- **Cannot raise a tier; degeneracy is flagged.** The worst included tier governs; a
  single-level factor reports a zero *with* a warning.
- **The line, restated.** Any feature that takes a real patient's measurement and returns a
  forecast or a therapy choice **does not get built.** Putting all the structural uncertainty
  on one ledger changes none of this; it makes the assembled risk legible, and stops there.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Variance-component core** | `onkos.budget.variance_components` (balanced two-way) + the landmark suite (§5); the v0.21 collapse. | The four components compute, sum to total, and every landmark passes. |
| **2 — Binding to the forecast** | `model_selection_budget` over `compare().included` × the eligible survival links, reusing `ensemble_samples`; the `dominant` axis; degeneracy/tier guards. | The NSCLC budget computes with all four components and names a dominant axis. |
| **3 — Surfaces** | `onkos budget` CLI; `to_dict` with the clinical-use flag; a stacked-budget figure + a CI-executed notebook. | The budget is reachable and visualized; the notebook runs in CI. |
| **4 — Report & triage** | `onkos report` ranks contexts by their dominant structural axis and the parameter-vs-structural split (binned). | `onkos report` shows, per context, where standardization has the most leverage. |

Step 1 alone is a self-contained, citable contribution: an open, landmark-tested tool that
decomposes a composed oncology survival forecast into parameter noise and the variance
contributed by each structural modeling choice — turning the field's qualitative "it depends
on your assumptions" into a quantitative budget that names which assumption — is something the
field reasons about only narratively and nobody ships as a tested artifact.
