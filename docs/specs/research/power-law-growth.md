# Onkos — research spec: power-law (sub-exponential) growth — and the exponential-overestimation it exposes

**Status:** implemented in v0.41.0 (`growth_power_law` kernel + `growth_laws.power_law` record). This is
the design-of-record; written in the v0.1 house style. Values are illustrative and `unverified`; the
infrastructure is the contribution.

**The growth-law assumption is a silent model-selection choice at the very top of the chain — and the
field's convenient default is the one that overestimates.** v0.40 completed the *bounded* growth laws
with von Bertalanffy (surface-limited, carrying capacity). This adds the other, empirically best-
supported member: the **power-law** `dV/dt = a·V^p` (p < 1) — *unbounded* but sub-exponential. Benzekry
et al. (2014), comparing classical growth laws across many experimental tumor datasets, found the
power-law the best fit overall. Its consequence is sharp: assuming **exponential** growth (the default
in much of TGI modeling) systematically **overestimates** extrapolated tumor burden, because real growth
decelerates. This is the growth-layer analog of the v0.36 exposure-response dose-extrapolation axis — a
modeling choice invisible at the studied timepoint that dominates the extrapolation.

> Power-law growth is mechanistically a fractal/feeding-surface argument: the fraction of proliferating
> cells scales sub-linearly with volume, so `dV/dt ∝ V^p`, `p < 1`. Unlike von Bertalanffy it has no
> carrying capacity — the tumor grows forever, just ever-more-slowly per unit mass. It is the law that
> says "the tumor is bigger than exponential predicts early, and far smaller than exponential predicts
> late."

---

## 1. The model

```
dV/dt = a·V^p ,   p < 1            specific growth rate (1/V)dV/dt = a·V^{p-1}  (falls with size)
```

Separable (`V^{1-p}` is linear in t), giving the closed form

```
V(t) = ( V0^{1-p} + a·(1-p)·t )^{1/(1-p)}            (p ≠ 1; the p → 1 limit is exponential)
```

Its characteristic, analytically derivable properties:

- **Sub-exponential everywhere**: the specific growth rate `a·V^{p-1}` decreases monotonically with size
  — slower than *any* exponential, from the first cell.
- **Unbounded, no carrying capacity**: unlike logistic/Gompertz/von Bertalanffy, `V → ∞` (polynomially in
  `t^{1/(1-p)}`), so it is the *unbounded* sub-exponential law — distinct from von Bertalanffy.
- **Rate-matched overestimation**: matched to an exponential at the same baseline size *and* instantaneous
  rate (`kg = a·V0^{p-1}`), the power-law stays strictly below that tangent exponential for all `t > 0`,
  and the gap explodes on extrapolation (≈90× by two years for the illustrative `a = 0.16, p = 0.75`).

**Why this is the right work (and the right scope).** It (1) completes the growth-law family's
*unbounded sub-exponential* slot — the empirically most-supported law (Benzekry 2014), the one whose
absence was most conspicuous; (2) is a first-class closed-form **reference kernel** that round-trips to
every export and is landmark-validated, the second to use a fractional-power `rhs_infix` (von Bertalanffy
was the first); (3) carries a genuine, **sharp, non-fragile finding** — the exponential-overestimation —
that is the growth-layer member of the project's extrapolation/transportability thesis (cf. the v0.36 ER
dose-extrapolation axis); (4) is *safe by construction* — an unperturbed growth law, no patient data;
(5) is a different *kind* of contribution from the model-selection-analysis arc — a model, not an axis.

---

## 2. The result — the exponential default overestimates

Matched to the same baseline (`V0 = 10`) and the same instantaneous growth rate, the exponential and the
power-law (`a = 0.16, p = 0.75`) diverge on extrapolation (illustrative):

| Week | power-law `V` | exponential `V` (matched) | overestimate |
| --- | --- | --- | --- |
| 26 | 110 | 102 | ~1× (still tangent) |
| 52 | 330 | 1,040 | ~3× |
| 78 | 680 | 10,600 | ~16× |
| 104 | 1,240 | 115,800 | **~93×** |

The finding is sharp and robust (not a tuned near-tie): an exponential growth assumption fit to early
data — which looks fine while the tumor is small — inflates the two-year burden by ~90× relative to the
sub-exponential law that better fits real tumors. The growth-law choice is invisible at the studied
timepoint and dominates the extrapolation, exactly the silent-transport pattern Onkos exists to make
visible, now at the top of the chain (the growth law) rather than its middle (the survival metric) or its
entry (the exposure-response shape).

**The honest framing.** Onkos does not claim power-law is universally correct (Benzekry found it *best on
average*, not always); it adds it as the empirically-favored option and shows the consequence of the
exponential default. The exponent `p` is a declared, illustrative parameter, not fit here.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Tier B**, illustrative, `unverified`, like the other growth-law records; out-of-context use warns.
- **A growth law, not a fit.** `a`, `p` are declared; nothing is estimated.
- **Population / trajectory level only.** An unperturbed growth trajectory; no survival, no individual
  quantity.
- **Round-trip validated.** Binds to NONMEM/SBML/PharmML/rxode2/Pumas, checked analytic↔ODE and
  MathML-per-state; an export bug cannot ship silently.

---

## 4. Validation landmarks

`tests/test_landmarks.py` (scientific landmarks) and `tests/test_roundtrip.py` (export axis):

| Landmark | Condition |
| --- | --- |
| **Sub-exponential** | the specific growth rate `a·V^{p-1}` falls monotonically with size. |
| **Below the rate-matched exponential** | matched at baseline (`kg = a·V0^{p-1}`), `V_power(t) < V_exp(t)` for all `t > 0`, and `V_exp(104) > 20·V_power(104)` (≈90× here) — the exponential-overestimation. |
| **Closed form ↔ integration** | `V(t) = (V0^{1-p} + a(1-p)t)^{1/(1-p)}` matches SciPy integration of `rhs` to ~1e-6. |
| **MathML round-trip** | the fractional-power `rhs_infix` (`a·V^p`) re-parses from the generated SBML and evaluates against `rhs` to ~1e-6 (`ANALYTIC_RECORDS`). |
| **Cross-format** | NONMEM `$THETA`, rxode2/Pumas/SO parameter vectors all read back the dataset values. |

---

## 5. API, CLI, and surface

```python
import numpy as np, onkos
ds = onkos.load()
from onkos.export.registry import get_kernel, kernel_values
spec = get_kernel(ds["growth_laws.power_law"]); v = kernel_values(ds["growth_laws.power_law"]); v["V0"] = 10.0
spec.analytic(np.linspace(0, 104, 209), v)         # sub-exponential trajectory (unbounded)
onkos export --format sbml --output exports/sbml/  # round-trips like any ODE kernel
```

**No new module, CLI command, or export format** — a kernel + a record + the finding, surfaced through
the existing simulate/export paths, a growth-law extrapolation figure, and a CI-executed notebook.

---

## 6. Source anchors (methodological; values illustrative)

- **Power-law / comparative tumor growth.** Benzekry et al. (2014), *Classical Mathematical Models for
  Description and Prediction of Experimental Tumor Growth*, PLoS Comput Biol — the cross-dataset
  comparison finding power-law the best-fitting unperturbed law (curated against the existing
  `bruno-2020-review` tumor-dynamics review pending its own DOI-verified record).
- **Fractal / feeding-surface mechanism.** The sub-linear proliferating-fraction argument underlying
  `dV/dt ∝ V^p`.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not a fitted exponent.** `a`, `p` are illustrative (`p = 0.75`, Benzekry-typical); no data is fit.
- **Not a claim that power-law is always right.** It is added as the empirically-favored option; the
  finding is about the exponential default's overestimation, not power-law's universality.
- **Not a TGI model.** Unperturbed growth; coupling a kill term to power-law growth is a later addition.
- **Not the generalized Bertalanffy-Pütter** (`a·V^p − b·V^q`, two exponents): the shipped kernel is the
  pure power-law; the two-exponent family is a clean follow-on the same separable structure supports.

---

## 8. Safety & scope (unchanged hard line)

- **Population / trajectory level only.** An unperturbed growth trajectory; nothing is an individual
  prediction.
- **No therapy or prognosis.** A growth law describes tumor dynamics, not a patient outcome.
- **Tier B, cannot self-promote.** Illustrative, `unverified`; promotion needs source-PDF review (spec §9).
- **The line, restated.** Adding a growth law changes none of the project's hard boundaries.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Kernel** | `growth_power_law` (`analytic` closed form, `rhs`, `rhs_infix` `a·V^p`). | the analytic↔ODE and MathML round-trips pass. |
| **2 — Record** | `growth_laws.power_law` (`a = 0.16`, `p = 0.75`), tier B. | `onkos validate` passes; the record round-trips and exports. |
| **3 — Landmarks** | the sub-exponential, below-rate-matched-exponential, and closed-form landmarks. | `tests/test_landmarks.py` passes. |
| **4 — Surfaces** | a growth-law extrapolation figure + a CI-executed notebook; README growth-law update; changelog. | the exponential-overestimation is visualized and documented. |

Step 1–3 is the contribution: a canonical, closed-form, round-trip-and-landmark-validated power-law
growth kernel, the empirically best-supported unperturbed law. Step 4 makes its sharp consequence — the
exponential default overestimates extrapolated burden — visible and documented.
