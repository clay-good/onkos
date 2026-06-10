# Onkos — research spec: exposure-response model choice — the dose-extrapolation model-selection axis

**Status:** implemented in v0.36.0 (`onkos.dose_response`). This is the design-of-record; written in
the v0.1 house style. All values are illustrative and `unverified` by design; the infrastructure is the
contribution.

**The composed survival forecast begins with one upstream modeling choice the rest of Onkos has taken
as given: the exposure-response model.** Every chain so far started from a *given* drug-effect magnitude
(or a `drug_effect` scalar). But that effect is the output of an **exposure-response (ER) model** —
Emax, sigmoid-Emax, or power — that maps a drug exposure to an effect. Those shapes all fit the studied
dose comparably and **diverge when extrapolated** to a dose the trial did not study, which is precisely
what a dose-selection decision asks of them. This spec makes the ER-model choice the project's core
transportability thesis applied one layer upstream: anchor the shapes to agree at the studied dose, then
quantify how far their predicted effect — and the resulting population OS — diverge off it. A
dose-response model fit at one dose carries an unquantified model-selection risk the moment it is used to
pick another.

> This is the dose-response analog of the dataset's load-bearing idea (`derivation_context` +
> `transportability`): a model derived in one context, applied in another, with the predictive validity
> of the move unknown. Here the "context" is the *dose*. The ER model is fit at the studied exposure; the
> go/no-go question — can we de-escalate? do we need to escalate? — lives at exposures the model has
> never seen. Onkos makes the extrapolation, and its survival consequence, visible.

---

## 1. The problem this extends

`onkos.simulate` consumes either a `drug_effect` scalar or an `exposure` + `exposure_response` record,
where the ER record's kernel maps exposure `C` to effect `E`. The curated ER shapes are:

| kernel | form | extrapolation character |
| --- | --- | --- |
| `er_emax` | `E = Emax·C/(EC50+C)` | **saturating** — bounded above by `Emax` |
| `er_sigmoid_emax` | `E = Emax·Cᵍ/(EC50ᵍ+Cᵍ)` | **switch-like** — flat then steep then saturating |
| `er_power` | `E = slope·Cᶿ` | **unbounded** — no ceiling, grows without limit |

Fit to a single studied dose, all three can pass through the same point; their *shapes* — what happens
away from that point — are a modeling assumption. This module **re-anchors** each curated shape so it
passes through a reference `(c_ref, e_ref)` (the studied dose and its effect): the shape parameters
(`EC50`, `gamma`, `theta`) are kept from the curated record and the single scale parameter (`Emax` or
`slope`) is solved to hit `e_ref` at `c_ref`. The curves are then identical at `c_ref` and differ only
in how they extrapolate. The combined effect feeds the *existing* TGI → survival chain unchanged.

**Why this is the right deepening (and the right scope).** It (1) is the reserve list's
"exposure-response model as a further model-selection / budget factor" — here as a self-contained axis;
(2) is **pure post-processing** over the curated ER kernels (reusing the inverse-free forward maps), with
a landmark-tested calibration, no new dataset record, kernel, schema, or export; (3) is the **upstream**
complement to the survival-side axes (metric v0.25/v0.33, link-structure v0.34, additivity-reference
v0.35) — the first modeling link in the chain, until now unquantified; (4) is *safe by construction* — a
population OS over a re-anchored published shape, never a dose recommendation; (5) is the purest
expression of the project's transportability thesis, with the *dose* as the context that the model is
silently transported across.

---

## 2. The result — invisible at the studied dose, a model-selection axis off it

For Claret NSCLC first line, the three ER shapes anchored at `(c_ref = 150, e_ref = 1.0)`, illustrative:

| Dose (rel. to studied) | Emax `E` | power `E` | sigmoid `E` | effect spread | **OS spread** |
| --- | --- | --- | --- | --- | --- |
| 0.25× (38) — deep de-escalation | 0.40 | 0.31 | 0.12 | 0.28 | **19 wk** |
| 0.5× (75) — de-escalation | 0.67 | 0.55 | 0.40 | 0.27 | **14 wk** |
| **1× (150) — studied dose** | **1.00** | **1.00** | **1.00** | **0.00** | **0 wk** |
| 2× (300) — escalation | 1.33 | 1.80 | 1.60 | 0.47 | 5 wk |
| 4× (600) — escalation | 1.60 | 3.25 | 1.88 | 1.65 | 5 wk |

Two findings:

1. **The ER-model choice has zero consequence at the studied dose and a real one off it.** At `c_ref`
   the three shapes agree exactly (effect spread 0, OS spread 0) — the control that proves the divergence
   elsewhere is the *extrapolation*, not the models being globally different. Move off the studied dose
   and the OS prediction depends on which ER shape you assumed: up to ~19 weeks of median OS riding on a
   choice that is invisible in the data at the dose you studied.

2. **The risk is asymmetric — sharpest on de-escalation.** Downward extrapolation (lower dose) lands the
   effect on the *steep* part of the effect→survival relationship, so the ER-model choice moves OS most
   there (19 wk at quarter-dose). Upward extrapolation produces a larger *effect* spread (the unbounded
   power curve runs away — `E = 3.25` at 4×) but a smaller *OS* spread (the tumor is already controlled;
   survival saturates). De-escalation — the dose-finding question of "can we give less?" — is exactly
   where the unquantified ER-model risk bites hardest.

**The honest framing.** Onkos refits nothing and recommends no dose. The re-anchoring is an explicit
analysis ("if these dose-response *shapes* all explained the studied dose equally, how would they
extrapolate?"), the reference point is declared, and the headline is the *spread as a function of how far
you extrapolate*, not a preferred shape. The divergence-at-the-anchor of zero is the built-in control.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** The effect drives the underlying TGI → survival chain; the comparison
  carries that chain's propagated tier (C), and out-of-context transport floors to D. The ER-model choice
  cannot raise a tier.
- **A re-anchoring, not a fit.** The shapes and the reference point are declared; nothing is estimated
  from data, and no shape is endorsed.
- **Population / trial level only.** Population OS over a published model; never a per-patient prediction,
  never a dose or therapy recommendation. The comparison ranks *ER shapes under a dose*, not doses.
- **Default view untouched.** No record, kernel, schema, or export changes; every default artifact is
  byte-identical.

---

## 4. Validation landmarks

The calibration is closed-form; the binding adds the finding (`tests/test_dose_response.py`):

| Landmark | Condition |
| --- | --- |
| **Anchor is exact** | every re-anchored shape passes through `(c_ref, e_ref)` to machine precision — the construction that makes them indistinguishable at the studied dose. |
| **Calibration preserves monotonicity** | a calibrated curve is increasing in dose (more drug ⇒ more effect). |
| **Shapes diverge off the anchor** | anchored at one point, the saturating / unbounded / switch-like shapes give different effects above and below it. |
| **Rejects a non-positive anchor** | `c_ref ≤ 0` or `e_ref ≤ 0` raises. |
| **No divergence at the studied dose** | the OS spread at `c_ref` is 0 (the control). |
| **Divergence appears on extrapolation** | the max OS spread off the anchor is a real, nonzero number of weeks. |
| **De-escalation diverges most in OS** | the OS spread at `0.25·c_ref` exceeds the spread at `4·c_ref`. |
| **Effect spread grows below the anchor** | the effect spread at `0.25·c_ref` exceeds the (zero) spread at `c_ref`. |
| **Single model ⇒ no divergence** | one ER shape has nothing to disagree about; divergence is identically 0. |
| **Tier & guardrails** | the comparison rides the chain's tier; out-of-context transport floors to D; the result carries the clinical-use prohibition. |

---

## 5. API, CLI, and surface

```python
from onkos.dose_response import calibrated_er, compare_er_extrapolation

# A curated ER shape re-anchored to the studied dose.
f = calibrated_er(ds, "exposure_response.power_generic", c_ref=150.0, e_ref=1.0)
f(150.0)   # == 1.0 exactly; f(75.0), f(300.0) differ by shape

# The dose-extrapolation divergence across shapes, through the TGI -> survival chain.
cmp = compare_er_extrapolation(ds, "resistance.claret_2009.tgi", context=ctx, c_ref=150.0, e_ref=1.0)
cmp.reference_os_divergence   # ~0 — the anchored control
cmp.max_os_divergence         # weeks of OS riding on the ER-shape choice on extrapolation
cmp.os_divergence_at(37.5)    # the de-escalation case
```

```bash
onkos dose-response resistance.claret_2009.tgi --c-ref 150 --e-ref 1.0
```

**No new module dependencies, record, kernel, or export** — `onkos.dose_response` is a pure
post-processing module surfaced through a CLI command, a figure, and a CI-executed notebook.

---

## 6. Source anchors (methodological; values illustrative)

- **Emax / sigmoid-Emax exposure-response.** Holford & Sheiner (DOI 10.2165/00003088-198106060-00002) —
  the saturating and Hill forms re-anchored here.
- **Power / log-linear exposure-response.** The unbounded form used in many oncology ER analyses (no
  saturating ceiling), the divergent extrapolation contrast.
- **Exposure-response for dose selection.** The regulatory exposure-response framework (the basis for
  dose-finding and the "is the studied dose the right dose?" question this axis quantifies).

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not a dose recommender.** The comparison ranks ER *shapes* at a fixed dose grid, never doses or
  treatments; the dose grid is descriptive, not an optimization.
- **Not a refit.** The shapes are re-anchored, not estimated from data; no shape is preferred.
- **Not a budget factor (yet).** Folding the ER-model choice into the v0.26 model-selection budget as a
  third structural factor (beside TGI model and survival link) is a clean follow-on; this spec introduces
  the axis and its dose-extrapolation signature first.
- **Not new ER records.** It re-anchors the curated shapes; calibrating context- or drug-specific ER
  curves is a breadth follow-on.

---

## 8. Safety & scope (unchanged hard line)

- **Population / trial level only.** Population OS over a re-anchored published shape; nothing is an
  individual prediction.
- **No dose or therapy ranking.** The axis ranks ER shapes under a dose, never doses or treatments.
- **Cannot raise a tier, refits nothing.** The shapes and reference are declared; the underlying model's
  tier governs.
- **The line, restated.** Any feature that takes a real patient's exposure and returns a dose, an effect,
  or a prognosis **does not get built.** Making the ER-model choice explicit and its dose-extrapolation
  consequence computable changes none of this.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Calibration** | `calibrated_er` re-anchors each curated ER shape (Emax / sigmoid-Emax / power) to `(c_ref, e_ref)` by solving its scale parameter. | the anchor-is-exact and monotonicity landmarks pass. |
| **2 — The extrapolation axis** | `compare_er_extrapolation` runs every shape over a dose grid through the TGI → survival chain and reports the effect and OS spread per dose. | the §2 table holds; the OS spread is 0 at `c_ref` and nonzero off it. |
| **3 — The asymmetry finding** | the de-escalation-diverges-most result and the effect-vs-OS divergence contrast. | the down > up OS-spread landmark holds; default view byte-identical. |
| **4 — Surfaces** | CLI `onkos dose-response`, an ER-curves + OS-spread figure, a CI-executed notebook, README + changelog + API contract. | the dose-extrapolation axis is visualized and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested re-anchoring of the curated ER
shapes with the anchor-is-exact identity as its correctness control. Step 3 is the payload: showing that
the ER-model choice is invisible at the studied dose but a real model-selection axis off it — with the
risk sharpest on de-escalation, the dose-finding question — is the quantitative core of this spec,
shipped as a tested artifact.
