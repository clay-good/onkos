# Onkos — research spec: von Bertalanffy (surface-area-limited) growth — completing the growth-law family

**Status:** implemented in v0.40.0 (`growth_von_bertalanffy` kernel + `growth_laws.von_bertalanffy`
record). This is the design-of-record; written in the v0.1 house style. Values are illustrative and
`unverified` by design; the infrastructure is the contribution.

**The spec's declared growth-law envelope had a hole.** Spec §2 lists the in-scope unperturbed
growth laws as "exponential, logistic, Gompertz, Simeoni exponential→linear, **von
Bertalanffy/power-law**." Onkos shipped the first four but not the last. This release fills the gap
with the canonical **von Bertalanffy** (surface-area-limited / ontogenetic) growth model — a different
*kind* of work from the eighteen model-selection-analysis versions that preceded it: a new
first-class reference kernel, like Gompertz or Norton-Simon, not another post-processing axis.

> Von Bertalanffy is the growth law with a mechanistic story: a tumor proliferates at its **surface**
> (where nutrients and oxygen reach it, `∝ V^{2/3}`) and loses mass throughout its **volume**
> (`∝ V`). The competition `dV/dt = a·V^{2/3} − b·V` is sub-exponential from the first cell and
> saturates at a carrying capacity `V∞ = (a/b)³`. It is the slot in the growth-law family between the
> unbounded exponential and the capacity-limited logistic/Gompertz, and it is conspicuously the one
> the dataset was missing.

---

## 1. The model

```
dV/dt = a·V^{2/3} − b·V          a = surface-proliferation coefficient, b = volume-loss coefficient
```

The substitution `u = V^{1/3}` linearizes it: `du/dt = (a − b·u)/3`, so

```
V(t) = ( c + (V0^{1/3} − c)·e^{−b·t/3} )³ ,     c = a/b = V∞^{1/3}
```

a clean closed form (no numerical integration needed for the single-state case). Its characteristic
properties, each an analytically derivable landmark:

- **Carrying capacity** `V∞ = (a/b)³`, where `dV/dt = 0` (exact).
- **Inflection** (peak absolute growth rate) at `V = (2a/3b)³ = (2/3)³·V∞ ≈ 0.296·V∞` — strictly
  *below* the logistic's `V∞/2`, the surface-limited signature.
- **Sub-exponential everywhere**: the specific growth rate `(1/V)dV/dt = a·V^{−1/3} − b` falls
  monotonically with size, from the first cell — unlike exponential (constant) and unlike
  logistic/Gompertz (which start near their intrinsic rate and only bend near capacity).

**Why this is the right work (and the right scope).** It (1) closes a **declared-in-scope gap**
(spec §2's named growth-law family) — a housekeeping-grade completion with real modeling content;
(2) is a first-class **reference kernel** (closed form + `rhs` + `rhs_infix`), so it round-trips to
every export and is landmark-validated exactly like the other growth laws — the first kernel whose
`rhs_infix` uses a fractional power (`V^{2/3}`), which exercises the MathML `power` round-trip path
that was implemented but unused; (3) is *safe by construction* — an unperturbed growth law, no patient
data; (4) adds another option to the growth-law model-selection axis (different laws extrapolate the
same early data differently — von Bertalanffy predicts a slower late approach to capacity than
logistic/Gompertz); (5) is genuinely a *different kind* of contribution after the model-selection arc,
which the project needed.

---

## 2. The result — the growth-law family, completed

From a common baseline (`V0 = 10`, `V∞ = 200`), the four single-state laws separate by their
**specific growth rate** signature `(1/V)dV/dt` — the fingerprint a correct kernel must reproduce:

| Law | specific growth rate vs size | behavior |
| --- | --- | --- |
| exponential | constant (`kg`) | unbounded |
| logistic | linear decline (`kg(1 − V/V∞)`) | sigmoid to `V∞`, inflection at `V∞/2` |
| Gompertz | log decline (`kg·ln(V∞/V)`) | sigmoid to `V∞`, inflection at `V∞/e` |
| **von Bertalanffy** | **`a·V^{−1/3} − b`** (falls as `V^{−1/3}`) | **sub-exponential from the start, inflection at `(2/3)³·V∞`** |

The finding is modest and honest: von Bertalanffy is the surface-limited member of the family, with a
distinct sub-exponential signature and a lower inflection than the other capacity-limited laws. As a
growth-law option it extrapolates the same early data to a *slower* late approach to capacity — one
more entry in the growth-law model-selection axis, now with the family complete.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Tier B**, like the other illustrative growth-law records (an established structural form with an
  illustrative value), `unverified`. Out-of-context use warns (the growth-law `out_of_context_action`).
- **A growth law, not a fit.** `a`, `b` are declared, illustrative; nothing is estimated.
- **Population / trajectory level only.** An unperturbed growth trajectory; no survival, no individual
  quantity.
- **Round-trip validated.** The kernel binds to NONMEM/SBML/PharmML/rxode2/Pumas and is checked
  analytic↔ODE and MathML-per-state like every ODE kernel; an export bug cannot ship silently.

---

## 4. Validation landmarks

`tests/test_landmarks.py` (the scientific-landmark axis) and `tests/test_roundtrip.py` (the export
axis):

| Landmark | Condition |
| --- | --- |
| **Carrying capacity** | `dV/dt = 0` exactly at `V = (a/b)³` (and `(a/b)³ ≈ 200`, the record's illustrative `V∞`). |
| **Surface-limited inflection** | the absolute growth rate peaks at `(2a/3b)³`, strictly below `V∞/2` (unlike logistic). |
| **Sub-exponential** | the specific growth rate falls monotonically with size — surface-limited from the first cell. |
| **Analytic ↔ ODE** | the closed form matches SciPy integration of `rhs` to ~1e-4 (`ANALYTIC_RECORDS`). |
| **MathML round-trip** | the `rhs_infix` (with the fractional power `V^{2/3}`) re-parses from the generated SBML and evaluates against `rhs` to ~1e-6 — the first kernel to exercise the `power` MathML path. |
| **Cross-format** | NONMEM `$THETA`, rxode2/Pumas/SO parameter vectors all read back the dataset values. |

---

## 5. API, CLI, and surface

The von Bertalanffy law is an ordinary growth-law record — every existing surface works on it:

```python
import numpy as np, onkos
ds = onkos.load()
from onkos.export.registry import get_kernel, kernel_values
spec = get_kernel(ds["growth_laws.von_bertalanffy"]); v = kernel_values(ds["growth_laws.von_bertalanffy"])
v["V0"] = 10.0
spec.analytic(np.linspace(0, 200, 401), v)   # closed-form trajectory; V∞ = (a/b)³

onkos export --format sbml --output exports/sbml/   # round-trips like any ODE kernel
```

**No new module, CLI command, or export format** — the contribution is a kernel + a record + the
completed family, surfaced through the existing simulate/export paths, a growth-law figure, and a
CI-executed notebook.

---

## 6. Source anchors (methodological; values illustrative)

- **Von Bertalanffy growth.** von Bertalanffy, L. (1957), *Quantitative laws in metabolism and growth*
  — the surface-vs-volume (`V^{2/3}` vs `V`) ontogenetic growth law.
- **Power-law / comparative tumor-growth modeling.** Benzekry et al. (2014), *Classical Mathematical
  Models for Description and Prediction of Experimental Tumor Growth* — the modern comparison of
  exponential / logistic / Gompertz / power-law / von Bertalanffy fits (the methodological anchor for
  treating these as a family); curated against the existing `bruno-2020-review` tumor-dynamics review
  citation pending its own DOI-verified record.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not a fitted law.** `a`, `b` are illustrative, set so `V∞ = (a/b)³ ≈ 200` to match the family's
  conventional carrying capacity; no data is fit.
- **Not the generalized Richards / Bertalanffy-Pütter family.** The shipped kernel is the classic
  `p = 2/3` surface law; a general exponent `p` (the power-law / Richards family) is a clean follow-on
  the same closed-form substitution supports.
- **Not a TGI model.** This is unperturbed growth; coupling a drug-effect/kill term to the von
  Bertalanffy law (a von-Bertalanffy-based TGI model) is a separate, later addition.

---

## 8. Safety & scope (unchanged hard line)

- **Population / trajectory level only.** An unperturbed growth trajectory; nothing is an individual
  prediction.
- **No therapy or prognosis.** A growth law describes tumor dynamics, not a patient outcome.
- **Tier B, cannot self-promote.** Illustrative, `unverified`; promotion needs source-PDF review per
  spec §9.
- **The line, restated.** Adding a growth law changes none of the project's hard boundaries.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Kernel** | `growth_von_bertalanffy` (`analytic` closed form, `rhs`, `rhs_infix` with `V^{2/3}`). | the analytic↔ODE and MathML round-trips pass (incl. the fractional-power path). |
| **2 — Record** | `growth_laws.von_bertalanffy` (`a`, `b`; `V∞ = (a/b)³ ≈ 200`), tier B. | `onkos validate` passes; the record round-trips and exports. |
| **3 — Landmarks** | the carrying-capacity, surface-limited-inflection, and sub-exponential landmarks. | `tests/test_landmarks.py` passes. |
| **4 — Surfaces** | a growth-law-family figure + a CI-executed notebook; README kernel-taxonomy + growth-law update; changelog. | the completed family is visualized and documented. |

Step 1–3 alone is the contribution: a canonical, closed-form, round-trip-and-landmark-validated growth
kernel that fills the spec's declared growth-law family. Step 4 documents it and visualizes the family,
now complete.
