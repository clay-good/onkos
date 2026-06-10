# Onkos — research spec: model discriminability — can a trial even tell the competing models apart?

**Status:** implemented in v0.38.0 (`onkos.discriminability`). This is the design-of-record; written in
the v0.1 house style. All values are illustrative and `unverified` by design; the infrastructure is the
contribution.

**The model-selection arc quantified how much the survival forecast depends on the modeling choice. This
closes the loop: it asks whether a trial could ever resolve that choice.** Across seventeen versions Onkos
showed that the TGI model, the resistance mechanism and origin, the bridge metric, and the survival
structure each move the predicted OS. The natural, and sharpest, follow-up is the *identifiability* of the
model choice itself: given two models' population OS curves, **what trial would it take to distinguish
them?** When the answer is tens of thousands of events, the model choice is **practically unidentifiable
from the trial** — it can only be assumed, not resolved by the data. That is precisely what makes the
silent-transport risk silent, and here it is quantified, in events.

> This is the model-level twin of `onkos.identify` (v0.22, *can a trial estimate this parameter?*) and
> `onkos.design` (v0.31, *what schedule estimates it best?*). Those ask about parameter identifiability
> within a model; this asks about discriminating *between* models. The unifying message: a quantity —
> a parameter or a whole model — that the trial cannot resolve is one whose uncertainty must be carried,
> not wished away.

---

## 1. The problem this closes

Distinguishing two survival curves is a logrank power calculation. The number of events required to detect
a difference at power `1−β` and two-sided level `α`, for 1:1 allocation, is Schoenfeld's formula:

```
d = 4 (z_{1-α/2} + z_{1-β})² / (ln HR)²
```

where `HR` is the hazard ratio between the curves. `required_events` implements it (it depends only on
`|ln HR|`, so it is symmetric in `HR ↔ 1/HR` and diverges as `HR → 1`). The HR between two simulated
population OS curves is the **follow-up-horizon hazard ratio** (`horizon_hazard_ratio`): the ratio of their
cumulative hazards `H = −ln S` at the trial horizon — exact under proportional hazards, a horizon-average
otherwise. `model_discriminability` runs this over every pair of a context's eligible TGI models and flags
the pairs that need an infeasible trial.

**Why this is the right closer (and the right scope).** It (1) is the rigorous resolution of the whole
divergence arc — not another divergence, but its *detectability*; (2) is **pure post-processing** over the
existing OS curves with a landmark-tested closed-form power core, no new dataset record, kernel, schema, or
export, so every default artifact is byte-identical; (3) completes the identifiability family — parameter
(v0.22), design (v0.31), now **model** discriminability; (4) is *safe by construction* — a trial-level
power calculation over published model structures, the same posture as v0.22/v0.31, no individual quantity;
(5) reframes the load-bearing idea: the silent model-selection risks are silent because they are
practically unidentifiable from a surrogate-driven trial.

---

## 2. The result — the resistance choice is practically unidentifiable under week-8 OS

For NSCLC first line, the required events to distinguish each model pair's OS curves under the default
week-8 survival link (power 0.8, two-sided α 0.05), illustrative:

| Model pair | horizon HR | required events | verdict |
| --- | --- | --- | --- |
| two-population vs Wang (biexp) | 0.48 | 58 | feasible |
| acquired vs Wang | 0.50 | 64 | feasible |
| Claret vs Wang | 0.50 | 67 | feasible |
| Norton-Simon vs two-population | 1.89 | 78 | feasible |
| Norton-Simon vs acquired | 1.82 | 87 | feasible |
| Norton-Simon vs Claret | 1.79 | 92 | feasible |
| Norton-Simon vs Wang | 0.90 | 3,127 | infeasible |
| **Claret vs two-population** | **1.05** | **11,807** | **infeasible** |
| **acquired vs two-population** | **1.03** | **27,031** | **infeasible** |
| **Claret vs acquired** | **1.02** | **102,694** | **infeasible** |

Two findings:

1. **The silent model-selection risks are silent because they are practically unidentifiable.** The pairs
   that diverge only in the regrowth *tail* — the resistance *mechanism* (Claret's phenomenological λ vs
   the two-population subclone, v0.24) and the resistance *origin* (acquired vs pre-existing, v0.32) — need
   **10⁴–10⁵ events** to distinguish under the week-8 surrogate. No realistic trial reaches that. The
   v0.24/v0.32 qualitative observation ("the week-8 surrogate is nearly blind to the resistance-model
   choice") is now a number: you would need ~12,000 to ~100,000 events. The model choice cannot be resolved
   by the data; it can only be assumed, with its tier and transportability attached.

2. **The risk is concentrated exactly where the trial cannot look.** The pairs that differ in early
   *shrinkage* (a deep early responder vs a minimal or complete one) are easily distinguished — ~60–90
   events, a small trial. And the survival-*metric* choice (week-8 vs k_g for the **same** model) produces
   a large OS swing detectable in <500 events. So the consequences a surrogate-driven trial *can* see (early
   shrinkage, the metric's effect) are identifiable; the model choice the surrogate is *blind* to (the tail
   mechanism) is not. The unquantified risk lives precisely in the blind spot.

**The honest framing.** This is a power calculation over illustrative models, not a real trial design. The
horizon HR is one summary of a possibly non-proportional comparison (documented as such), and the
feasibility bounds (≈500 / 3000 events) are declared phase-3-scale references, not clinical thresholds. The
headline is robust to those choices: the resistance-model pairs need orders of magnitude more events than
the shrinkage-distinct pairs.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** The curves carry their propagated tier (out-of-context models floor to D
  upstream in `compare`); a discriminability analysis cannot raise a tier.
- **A power calculation, not a trial.** The design parameters are declared; nothing is fit, and no trial is
  designed or recommended.
- **Design / trial level only.** Required events to distinguish *population* curves — the same posture as
  `identify`/`design`. Never an individual quantity, never a recommendation.
- **Default view untouched.** No record, kernel, schema, or export changes; every default artifact is
  byte-identical.

---

## 4. Validation landmarks

The power core is closed-form; the binding adds the finding (`tests/test_discriminability.py`):

| Landmark | Condition |
| --- | --- |
| **Schoenfeld benchmark** | `HR = 0.5` at 80% power, two-sided α 0.05 needs ~65 events (the textbook logrank value). |
| **Symmetric in HR↔1/HR** | `required_events(HR) == required_events(1/HR)` (depends only on `|ln HR|`). |
| **Identical ⇒ infinite** | `HR = 1` (or 0, or nan) needs infinite events. |
| **Smaller divergence ⇒ more events** | monotone: `HR` nearer 1 needs more events. |
| **More power / smaller α ⇒ more events** | the power and significance levers move events the right way. |
| **Horizon HR recovers PH** | for genuinely proportional curves `S_b = S_a^HR`, the horizon ratio recovers `HR` exactly; identical curves give `HR = 1`. |
| **Resistance models indistinguishable** | under week-8, the Claret/two-population/acquired pairs need >3000 events, while Claret-vs-Norton needs <500 — and the resistance pairs need >10× the shrinkage-distinct pairs. |
| **Several pairs flagged** | ≥3 of the NSCLC pairs are practically indistinguishable; not all are. |
| **Metric consequence is detectable** | the week-8-vs-k_g swing (same model) needs <500 events — the contrast. |
| **Tier & clinical-use** | the propagated tier rides through; infinite events serialize as `null`, not a NaN/inf literal. |

---

## 5. API, CLI, and surface

```python
from onkos.discriminability import required_events, horizon_hazard_ratio, model_discriminability

required_events(0.75)                       # events to distinguish a HR=0.75 difference
horizon_hazard_ratio(curve_a, curve_b)      # follow-up-horizon HR between two OS curves

md = model_discriminability(ds, context=ctx)   # pairwise required events for the context's models
md.n_indistinguishable                          # how many pairs need an infeasible trial
md.indistinguishable_pairs                      # the model choices the data cannot resolve
```

```bash
onkos discriminability --tumor-type NSCLC --line first
```

**No new module dependencies beyond scipy (already required), record, kernel, or export** —
`onkos.discriminability` is a pure post-processing module surfaced through a CLI command, a figure, and a
CI-executed notebook.

---

## 6. Source anchors (methodological; values illustrative)

- **Logrank sample size.** Schoenfeld, D.A. (1981/1983), the required-events formula for the proportional-
  hazards logrank test — the power core.
- **Restricted-mean / non-proportional caveats.** The literature on non-proportional-hazards trial design
  — context for the horizon-average HR used when the curves are not proportional (e.g. the joint link).
- **Parameter identifiability.** The Onkos `identify` (v0.22) and `design` (v0.31) modules — the
  within-model identifiability family this extends to between-model discrimination.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not a trial design.** It reports required events to distinguish models; it does not design, size, or
  recommend a real trial.
- **Not a full power simulation.** The shipped core is the closed-form Schoenfeld events from a horizon HR;
  a simulation-based logrank power (with accrual, censoring, non-PH weighting) is a clean follow-on.
- **Not external validation.** It quantifies *internal* discriminability (can a trial tell the models'
  predicted curves apart), not whether either model is correct — that is the `predictive_performance`
  external-validation field's job.
- **Not individual.** No per-patient survival, no event times emitted; the calculation is over population
  curves.

---

## 8. Safety & scope (unchanged hard line)

- **Design / trial level only.** Required events to distinguish population curves; nothing is an individual
  prediction.
- **No trial recommendation.** The analysis characterizes detectability; it does not size or recommend a
  study.
- **Cannot raise a tier, fits nothing.** The power parameters are declared; the underlying models' tier
  governs.
- **The line, restated.** Any feature that takes a real patient's data and returns a prognosis, a trial
  recommendation, or a model verdict for that person **does not get built.** Quantifying how unidentifiable
  a model choice is changes none of this.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Power core** | `required_events` (Schoenfeld) + `horizon_hazard_ratio` (cumulative-hazard ratio). | the Schoenfeld benchmark and PH-recovery landmarks pass. |
| **2 — Pairwise discriminability** | `discriminating_events` / `model_discriminability` over a context's models under a survival link. | the §2 table holds; the resistance pairs are flagged infeasible. |
| **3 — The finding** | the resistance-choice-unidentifiable result and the detectable-metric contrast. | the resistance-indistinguishable and metric-detectable landmarks hold; default view byte-identical. |
| **4 — Surfaces** | CLI `onkos discriminability`, a curves + required-events-heatmap figure, a CI-executed notebook, README + changelog + API contract. | the discriminability question is visualized and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested logrank-power calculation over
the simulated OS curves, with the Schoenfeld benchmark as its correctness anchor. Step 3 is the payload:
showing that the resistance mechanism/origin choice needs 10⁴–10⁵ events to detect under the week-8
surrogate — that the silent model-selection risk is silent *because* it is practically unidentifiable from
the trial — is the quantitative core of this spec, and the rigorous close of the model-selection arc,
shipped as a tested artifact.
