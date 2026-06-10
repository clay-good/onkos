# Onkos — research spec: joint longitudinal–survival modeling (the current-value link)

**Status:** implemented in v0.34.0 (`onkos.joint`). This is the design-of-record; written in the
v0.1 house style. `α` and the baseline hazard are illustrative, declared parameters — never fitted
here; the infrastructure (and the structural contrast it makes computable) is the contribution.

**Every survival link Onkos ships assumes proportional hazards — and that is itself a modeling
choice.** The TGI→OS surrogate has, since v0.12, been **two-stage**: collapse the simulated tumor
trajectory to one scalar covariate — the week-8 change (v0.12), the growth-rate constant `k_g`
(v0.25), the integrated burden (v0.33) — then apply a parametric (Weibull) or Cox baseline with that
*static* covariate. A static covariate means the hazard ratio between two tumors is **constant over
time**: a proportional hazard. The joint longitudinal–survival model — the rigorous, two-stage-free
formulation the pharmacometric literature treats as the gold standard for linking a longitudinal
biomarker to survival — relaxes exactly that assumption. This spec adds its canonical **current-value**
link and shows that "two-stage vs joint" is a model-selection axis at the survival-link layer, with
the disagreement concentrated in precisely the resistance models whose tail the week-8 surrogate is
blind to.

> The two-stage surrogate and the joint model are the two standard ways to connect tumor dynamics to
> survival. The two-stage approach is ubiquitous because it is simple (extract a metric, fit a hazard);
> the joint model is the statistically rigorous one (it propagates the whole trajectory and its
> uncertainty into the hazard, avoiding the regression-dilution and immortal-time biases of a landmark
> covariate). Onkos does not adjudicate — it makes both computable over the same trajectory and shows
> where they disagree.

---

## 1. The problem this extends

A two-stage Weibull link is `S(t) = exp(-(t/scale)^shape · exp(β·x))`, with `x` a **scalar** read off
the trajectory once. Its cumulative hazard is `H(t) = (t/scale)^shape · exp(β·x)` — the covariate
multiplies the baseline by a **time-constant** factor `exp(β·x)`. Two tumors with covariates `x₁`, `x₂`
have hazard ratio `exp(β(x₁−x₂))` at *every* time. That is the proportional-hazards assumption, and
the Cox link (v0.13) shares it (it only swaps the parametric baseline for a tabulated one).

The current-value joint link makes the hazard track the **instantaneous** tumor size:

```
λ(t) = λ₀(t) · exp(α · log(v(t)/y0)),     S(t) = exp(-∫₀ᵗ λ(u) du)
```

with `λ₀` the Weibull baseline hazard and `α` the association between log tumor size and log hazard.
The implementation integrates the time-varying hazard ratio as a Stieltjes sum against the **analytic**
baseline cumulative hazard `H₀(t) = (t/scale)^shape`, which makes it a *strict generalization*, exact
in two limits:

- a tumor held at baseline (`v ≡ y0`) ⇒ `HR ≡ 1` ⇒ the Weibull **baseline** survival exactly;
- a tumor held at constant `c·y0` ⇒ `HR ≡ c^α` constant ⇒ the **two-stage** Weibull-PH curve exactly,
  with `x = log c`, `β = α`. So the v0.33 burden link is the *constant-trajectory* special case, and
  more generally it is the joint model's "average-HR-in-the-exponent" approximation.

**Why this is the right deepening (and the right scope).** It (1) is the reserve list's "joint
longitudinal–survival modeling — the rigorous version" of the two-stage surrogate, and the natural
culmination of the bridge-metric arc (v0.25 metric-choice → v0.33 *static* whole-trajectory burden →
v0.34 *dynamic* whole-trajectory hazard); (2) is **pure post-processing** over the existing kernels —
a module with a landmark-tested pure core, no new dataset record, kernel, schema, or export, so the
default view and every artifact are byte-identical; (3) adds the one survival-link structure that is
**non-proportional**, beside the two PH structures already shipped (Weibull v0.12, Cox v0.13); (4) is
*safe by construction* — a population survival curve over a published trajectory, `α` declared not
fitted; (5) reinforces the load-bearing message: the regrowth tail the week-8 surrogate cannot see is,
in the joint model, a rising hazard that **inverts** the resistance-model ranking.

---

## 2. The result — a non-proportional hazard, and a ranking inversion

For NSCLC first line at unit drug effect, association `α = 1.0`, Weibull baseline from the default
week-8 link (`shape = 1.3`, `scale = 60`):

| Model | two-stage week-8 mOS | **joint** mOS | HR(8 wk) | HR(end) | PH-violation HR(end)/HR(8wk) |
| --- | --- | --- | --- | --- | --- |
| Norton-Simon (complete responder) | 58 | **n/r** (>260) | 0.71 | →0 | **< 1** (hazard keeps falling) |
| Claret (phenom. resistance) | 91 | **199** | 0.18 | 1.7 | ~10× |
| two-population (mechanistic) | 94 | **144** | 0.13 | 12.2 | ~96× |
| acquired resistance | 92 | 103 | 0.16 | 41.1 | **~255×** |
| Wang biexponential (minimal resp.) | 54 | 47 | 0.80 | 106.8 | ~134× |

Two findings:

1. **The hazard is non-proportional, and the violation scales with the regrowth tail.** Each two-stage
   link has a constant hazard ratio (PH-violation ≡ 1 by construction). The joint link's hazard ratio
   is suppressed during the deep early response (HR ≈ 0.13–0.18) then **rises by 10× to 255×** as the
   resistant clone regrows — largest for the acquired-resistance and two-population models. The complete
   responder is the mirror image: its hazard *keeps falling* (PH-violation < 1) as the tumor is
   eradicated. This time-varying hazard ratio is structurally inexpressible by any two-stage (PH) link,
   parametric or Cox.

2. **The ranking inverts, for a structural reason.** The week-8 two-stage surrogate ranks the deep-
   early-shrinking two-population model *above* the phenomenological Claret model (mOS 94 vs 91). The
   joint link, weighting the regrowth tail as a rising hazard, **inverts** this (Claret 199 vs
   two-population 144): the lighter-tail model wins. Same trajectories, opposite order — so the choice
   of survival-link *structure* (two-stage PH vs joint current-value) is a model-selection axis, and it
   is the structurally correct counterpart to the v0.25/v0.33 *metric* axis.

**The honest framing.** The absolute joint medians are not directly comparable to the two-stage medians
(they are different survival models with their own baseline anchoring); what is comparable, and is the
point, is the **shape** (a non-proportional hazard) and the **ranking** (inverted relative to week-8).
`α` is a declared, illustrative association, never fitted — a real joint model estimates it jointly from
longitudinal and survival data; Onkos simulates the structure forward and shows its consequences.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Inherits the trajectory's tier.** A joint analysis is post-processing over `simulate`; it carries
  the propagated tier and transport warnings unchanged. Out-of-context transport floors to D + warns.
- **Cannot raise a tier, cannot fit.** `α` and the baseline are declared inputs; nothing is estimated,
  so no result can promote a record or claim identifiability of the association.
- **Population / trial level only.** A population survival curve over a published trajectory; never an
  individual hazard, prognosis, or therapy ranking. The comparison ranks *models under a context*.
- **Default view untouched.** No record, kernel, schema, or export changes; every default artifact and
  every v0.12–v0.33 number is byte-identical.

---

## 4. Validation landmarks

The pure core (`current_value_survival`) is landmark-tested in isolation; the binding adds the finding
(`tests/test_joint.py`):

| Landmark | Condition |
| --- | --- |
| **Constant HR ⇒ two-stage exactly** | a constant hazard ratio recovers the two-stage Weibull-PH curve to machine precision — the strict-generalization check (PH is the constant-trajectory special case). |
| **Unit HR ⇒ Weibull baseline** | `HR ≡ 1` (tumor at baseline size) recovers the Weibull baseline survival exactly. |
| **α = 0 removes the association** | the hazard ratio is identically 1 and the survival is the baseline — no marker effect. |
| **Monotone, normalized** | `S(0)=1`, `S` non-increasing, and a uniformly larger HR lowers survival pointwise. |
| **Constant-size ⇒ burden link** | a constant-size trajectory's joint curve equals the two-stage Weibull-PH curve with `x=log c`, `β=α` — the analytic bridge to v0.33. |
| **Regrowing tumor ⇒ non-proportional** | a resistant-regrowth trajectory has `HR(8wk) < 1 < HR(end)` and `ph_violation ≫ 1`. |
| **Eradication ⇒ falling tail hazard** | a complete responder's HR keeps falling (`ph_violation < 1`); its joint OS exceeds every regrowing model's. |
| **Re-ranks vs two-stage** | the joint link produces ≥1 rank-discordant model pair relative to the week-8 surrogate, and inverts the Claret/two-population order. |
| **Tier & guardrails** | the propagated tier rides through; out-of-context transport floors to D + warns; the clinical-use prohibition is attached. |

---

## 5. API, CLI, and surface

```python
from onkos.joint import joint_survival, compare_joint_vs_two_stage, current_value_survival

# Current-value joint OS for one model, beside its two-stage counterpart.
j = joint_survival(ds, "resistance.nsclc_first_line.two_population", context=ctx, alpha=1.0)
j.median_os, j.two_stage_median_os     # the joint vs the week-8 two-stage median
j.hazard_ratio                         # the time-varying HR(t) — the non-proportional signature
j.ph_violation                         # HR(end)/HR(8wk): 1 for PH, ≫1 for a regrowing tumor

# Across a context's eligible TGI models: the ranking divergence and the PH-violation panel.
cmp = compare_joint_vs_two_stage(ds, context=ctx, alpha=1.0)
cmp.rank_discordant_pairs, cmp.max_ph_violation
```

```bash
onkos joint --tumor-type NSCLC --line first --alpha 1.0
```

**No new module dependencies, record, kernel, or export.** The contribution is the `onkos.joint`
module (a pure core + a binding), the finding, a figure, a CI-executed notebook, and a CLI command.

---

## 6. Source anchors (methodological; values illustrative)

- **Joint longitudinal–survival models.** The current-value association structure
  `λ(t) = λ₀(t)·exp(α·m(t))` linking a longitudinal marker `m(t)` to the hazard is the canonical joint
  model (Rizopoulos and the broader shared-parameter / joint-model literature) — the rigorous
  alternative to a two-stage landmark covariate, which it improves on by avoiding regression-dilution
  and immortal-time bias.
- **Two-stage TGI→OS surrogates.** Wang (2009) (DOI 10.1038/clpt.2009.64), Claret (2009)
  (DOI 10.1200/JCO.2008.21.0807), Stein (2008) (DOI 10.1634/theoncologist.2008-0075) — the
  proportional-hazards links this contrast is measured against.
- **Tumor-dynamics / resistance review.** Bruno et al. (2020) (DOI 10.1158/1078-0432.CCR-19-0287) —
  context for why the resistant regrowth tail (a rising hazard in the joint model) is the part the
  two-stage surrogate misses.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not a fitted joint model.** `α` and the baseline hazard are declared, illustrative inputs; Onkos
  simulates the structure forward and never estimates the association from data.
- **Not a record/kernel.** The current-value link needs the *whole* trajectory, so it does not fit the
  `analytic(t, x)` kernel contract; it is a post-processing module (like `budget`/`response`/`design`),
  which is also why it changes no export.
- **Not a new endpoint.** This is an OS link *structure*, not a new endpoint; PFS routes (v0.30) and
  ORR/DoR (v0.27–v0.28) are untouched.
- **Not a budget factor (yet).** Adding the joint link as a third level of the v0.26 budget's survival-
  link factor (beside the PH links) is a clean follow-on; this spec introduces the structure and its
  contrast first.
- **Not cross-context-calibrated.** The finding uses the NSCLC baseline; the structure applies to every
  context with a default Weibull OS link, and a per-context `α` study is a follow-on.

---

## 8. Safety & scope (unchanged hard line)

- **Population / trial level only.** A population survival curve over a published trajectory; nothing is
  an individual hazard or prognosis.
- **No therapy ranking.** The joint-vs-two-stage comparison ranks *models under a context*, never
  treatments.
- **Cannot raise a tier, cannot fit.** A joint analysis inherits the trajectory's tier and estimates
  nothing.
- **The line, restated.** Any feature that takes a real patient's tumor measurement and returns a
  hazard, a prognosis, or a therapy choice **does not get built.** Making the proportional-hazards
  assumption explicit and its violation computable changes none of this.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Pure core** | `current_value_survival(t, HR, shape, scale)` — Stieltjes integration against the analytic baseline cumulative hazard. | The constant-HR⇒two-stage and unit-HR⇒baseline landmarks pass to machine precision. |
| **2 — Binding** | `joint_survival(ds, record, context, α, baseline_link)` — simulate the trajectory, read the Weibull baseline, build the time-varying HR; tier/warnings ride through. | `α=0` recovers the baseline; out-of-context floors to D. |
| **3 — The finding** | `compare_joint_vs_two_stage` — the per-model PH-violation and the ranking inversion vs the week-8 surrogate. | The §2 table and the dynamics landmarks hold; default view byte-identical. |
| **4 — Surfaces** | a hazard-ratio + survival figure, a CI-executed notebook, a CLI command, README + changelog + API-contract. | The non-proportional-hazard axis is visualized and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested current-value survival link
that contains the proportional-hazards link as an exact special case. Step 3 is the payload: showing
that the joint model's hazard ratio rises with the resistant regrowth — a non-proportional hazard the
two-stage surrogate cannot represent — and that this inverts the resistance-model ranking, is the
quantitative core of this spec, shipped as a tested artifact.
