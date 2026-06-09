# Onkos — research spec: drug-combination interaction models

**Status:** implemented in v0.23.0 (`onkos.interaction`; steps 1–4 of §10). This is
the design-of-record; the methodological source anchors of §7 are documented but their
Crossref-verified citation curation is still pending, honest by design. Written in the
v0.1 house style; every effect magnitude and interaction parameter is illustrative and
`unverified` by design — the infrastructure is the contribution.

**The assumption that sinks combination programs.** Oncology is overwhelmingly
*combination* therapy, yet a composed survival forecast for a combination silently
depends on one unmeasured modeling choice: *how do the two drugs' effects combine?*
Highest-single-agent, additive (Bliss/Loewe null), or synergistic — these give very
different predicted benefits from the *same* single-agent activity, and the difference
is routinely assumed rather than measured. This module makes the **interaction model a
first-class, quantified model-selection axis**: it combines two drug effects under each
declared interaction rule, propagates the result through the existing TGI → survival
chain, and reports how much the predicted outcome depends on which interaction you
assumed — never inventing a synergy value the data cannot support.

> The project's load-bearing idea is that a silent modeling choice — which model, which
> context — determines a billion-dollar outcome, so the choice should be visible and its
> consequence measured. The drug-effect subsystem already makes the *kill mechanism* a
> divergence axis (log-kill vs Norton-Simon). The **interaction model is the same move
> one layer up**: two trials can show identical single-agent activity yet predict
> different combination outcomes because one assumed synergy and the other additivity.
> Onkos refuses to assume; it shows the spread and labels synergy an *assumption*.

---

## 1. The problem this extends

Onkos composes a survival forecast from a growth law + a drug effect + resistance + an
exposure-response link + a survival link, and quantifies three uncertainty axes over
that composition (parameter, model-selection, design identifiability). All of it assumes
a **single** drug effect `E`. Real regimens combine agents, and the combination's net
effect is not given by the dataset — it is a *modeling assumption* with first-class
consequences.

| Combination question | Status before this spec |
| --- | --- |
| What is the net effect of drug A + drug B at given single-agent activities? | ⚠️ Undefined — `E` is a scalar; there is no way to combine two. |
| How much does the predicted outcome depend on the **interaction model** (HSA / additive / synergy)? | ⚠️ Not surfaced — the interaction assumption is invisible. |
| Can Onkos distinguish synergy from additivity? | ✅ (the honest answer) **No** — and it should say so, not invent a synergy parameter. |

The gap is exactly the project's named failure mode applied to combinations: a silent,
outcome-determining modeling choice presented as if it were settled.

**Why this is the right deepening (and the right scope).** It (1) advances the project's
own thesis — adds a new, explicitly-declared model-selection axis rather than breadth;
(2) is pure post-processing over the existing kernels — the combined effect feeds the
*same* TGI → survival chain, so no new ODE kernel, no dataset subsystem, near-zero schema
change (mirrors `onkos.combine` and `onkos.identify`); (3) has direct precedent — the
Bliss-independence, Loewe-additivity, HSA, and interaction-index references are the
field's standard combination nulls (§7); (4) is *safe by construction* — it simulates one
combination *regimen* under different interaction *assumptions* at the population level,
never ranks regimens and never emits an individual prediction; and (5) sharpens the honest
message: its headline output is *"the predicted benefit depends this much on an
interaction assumption you have not measured,"* which is the opposite of false precision.

---

## 2. The combination framework

Fix a TGI model, a context, and two single-agent **effect magnitudes** `E_A, E_B ≥ 0`
(in the kernel's drug-effect units — the same scalar that drives the kill term). An
**interaction model** maps the pair to a combined effect `E_AB` that then drives the
existing TGI → survival simulation unchanged:

```
E_AB = combine(E_A, E_B ; model, ψ)        then        simulate(..., drug_effect = E_AB)
```

The declared interaction models:

| Model | `E_AB =` | Meaning |
| --- | --- | --- |
| `hsa` | `max(E_A, E_B)` | **Highest single agent** — the most conservative null; the combination does no better than its better component. |
| `additive` | `E_A + E_B` | **Bliss-independence / effect-additive null.** For a log-linear (exponential) kill process the two surviving fractions multiply, `e^{−E_A·Δt}·e^{−E_B·Δt} = e^{−(E_A+E_B)·Δt}`, so *Bliss independence is exactly additive kill rates* — an identity Onkos states rather than hides (§6, landmark). |
| `greco` | `E_A + E_B + ψ·√(E_A·E_B)` | **Interaction-index** (Greco). `ψ = 0` ≡ additive; `ψ > 0` synergy; `ψ < 0` antagonism. The synergy is a *single declared parameter*, never fitted from the dataset. |

**Why effect-level, and the boundary we name.** Onkos combines at the **effect-magnitude
(kill-rate) level**, the quantity its kernels actually consume. Full **Loewe dose-
additivity** is defined at the *dose* level through the two agents' dose-response curves
(it solves the dose-equivalence equation `C_A/EC50_A + C_B/EC50_B = 1` at iso-effect).
Onkos *has* exposure-response curves (`er_emax`, with `Emax`/`EC50`), so dose-level Loewe
is a clean, noted extension (§8); the v0.x scope is effect-level combination, which is
sufficient to make the interaction-model-selection risk legible and is labeled as such.

**The honesty boundary — synergy is an assumption, not a finding.** The combined effect
`E_AB` is a *modeling assumption*, and every result says so: `ψ` is a declared input
(default 0, the additive null), never derived from the dataset, and a non-zero `ψ` carries
a `synergy_is_an_assumption` warning. Onkos does not estimate synergy — distinguishing
synergy from additivity requires a combination trial designed for it, and asserting it
without one is precisely the over-optimism this module exists to expose.

---

## 3. The divergence view — the interaction model as a model-selection axis

The headline is the **interaction-model divergence**: for one combination (the same
`E_A, E_B`), simulate the tumor-size and population-OS trajectories under every
interaction model and report how much the survival prediction moves.

```
interaction_divergence = max pointwise spread of the OS curves across {hsa, additive, greco±}
median_os_range        = (min, max) median OS across the interaction models
```

This is the combination-therapy analog of the virtual-trial divergence view (which spans
*TGI models*) and the kill-mechanism axis (which spans *kill mechanisms*): a number that
says *of everything you might predict for this combination, how much is just the
interaction assumption you have not pinned down.* A large divergence is the signal that a
program's projected combination benefit rests on an unverified synergy assumption.

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **The underlying TGI model's tier governs.** A combination simulation inherits the
  propagated tier of the TGI → survival chain it runs through (worst-input-wins, including
  the out-of-context transport floor). The interaction model does not *raise* a tier; an
  assumed synergy never makes a forecast more trustworthy.
- **The combined effect is flagged as an assumption.** `E_AB` is an illustrative effect
  magnitude, not a curated parameter; a non-zero `ψ` is surfaced as such. No synergy value
  is presented as a measured quantity.
- **Population / regimen level only.** The simulation is of *one combination regimen* at
  the trial level under stated interaction assumptions. It does **not** rank regimens
  (combination vs monotherapy is shown descriptively, never as a recommendation), does
  **not** compare drugs to choose between them, and emits **no** individual-level output.
- **Single-agent degeneracy is exact.** If either effect is zero the combination reduces to
  monotherapy under *every* interaction model — no spurious interaction is manufactured
  (a landmark, §6).

---

## 5. Reference kernel & validation landmarks

No new ODE kernel — combination is algebra over the effect that feeds the existing kernels.
The interaction math gets its own **landmark suite** (in the spirit of `test_combine.py`
and `test_identifiability.py`): closed-form properties of the combination rules themselves.

| Landmark | Closed-form condition |
| --- | --- |
| **Additive null** | `greco(ψ = 0) = additive = E_A + E_B`. |
| **Bliss ≡ additive identity** | Combining fractional survivals `(1−f_A)(1−f_B)` with `f = 1−e^{−E}` yields effect `E_A + E_B` — Bliss independence equals additive kill rates for log-linear TGI. |
| **Ordering** | For `E_A, E_B ≥ 0`: `hsa ≤ additive`, and `greco(ψ<0) ≤ additive ≤ greco(ψ>0)`. |
| **Synergy monotonicity** | `E_AB` is non-decreasing in `ψ`. |
| **Single-agent degeneracy** | `E_B = 0 ⇒ E_AB = E_A` for *every* model (no manufactured interaction). |
| **Symmetry** | `combine(E_A, E_B) = combine(E_B, E_A)` for every model. |
| **Monotonicity** | `E_AB` non-decreasing in each of `E_A, E_B`. |
| **Antagonism floor** | `greco` is clamped at 0 (a combined effect is never negative). |
| **Divergence sign** | A combination with `ψ ≠ 0` and both effects positive yields a strictly positive interaction divergence; `E_B = 0` yields zero divergence. |
| **Tier passthrough** | A combination trajectory's tier equals the underlying TGI → survival chain's propagated tier (interaction cannot raise it). |

This mirrors the project's validation discipline: round-trip proves exports match the
kernel; landmarks prove a kernel *is* the model it names; here the landmarks prove the
combination layer *is* the standard interaction nulls and a monotone interaction index —
not an unconstrained synergy knob.

---

## 6. API, CLI, and surface

**Python.** A new `onkos.interaction` module:

```python
from onkos.interaction import combine_effects, simulate_combination, compare_interactions

E_AB = combine_effects(0.6, 0.6, model="greco", psi=0.5)     # the pure interaction math

tr = simulate_combination(ds, "resistance.claret_2009.tgi",
                          context=dict(tumor_type="NSCLC", line="first"),
                          effect_a=0.6, effect_b=0.6, interaction="additive")
tr.tier, tr.median_os, tr.warnings        # the TGI→survival chain, driven by E_AB

cmp = compare_interactions(ds, "resistance.claret_2009.tgi", context=ctx,
                           effect_a=0.6, effect_b=0.6, psi=0.5)
cmp.combined_effects        # {hsa, additive, greco+ψ, greco−ψ} -> E_AB
cmp.median_os               # per interaction model
cmp.os_divergence           # how much predicted OS depends on the interaction assumption
cmp.median_os_range, cmp.warnings
cmp.to_dict()               # carries the clinical-use prohibition + synergy-is-assumption note
```

**CLI.**

```bash
onkos interactions resistance.claret_2009.tgi --effect-a 0.6 --effect-b 0.6 --psi 0.5
onkos interactions resistance.claret_2009.tgi --effect-a 0.6 --effect-b 0.6 --json
```

prints the combined effect and median OS under each interaction model, the OS divergence,
and the synergy-is-an-assumption warning.

**No new export model.** Combination is an analysis *over* the drug effect that drives a
model, not a model; it adds no NONMEM/SBML/PharmML surface (the underlying TGI model
already exports). The virtual-trial JSON of a combination run is a standard `Trajectory`
serialization with the combined effect recorded.

---

## 7. Source anchors (methodological; DOIs added at curation time)

Well-established methods, not Onkos parameters; added to `dataset/citations/` through the
normal Crossref/PubMed-verified pipeline, honest until a human confirms each.

- **Bliss independence.** Bliss, C.I. (1939), *The toxicity of poisons applied jointly*,
  Annals of Applied Biology — the probabilistic-independence combination null.
- **Loewe additivity.** Loewe, S. (1953), *The problem of synergism and antagonism of
  combined drugs*, Arzneimittel-Forschung — dose-additivity and the isobologram, the
  reference Onkos's effect-level additive approximates and names the boundary to.
- **Highest single agent.** Berenbaum, M.C. (1989), *What is synergy?*, Pharmacological
  Reviews — the survey that frames HSA, Bliss, and Loewe as distinct nulls.
- **Interaction index / response-surface synergy.** Greco, Bravo & Parsons (1995), *The
  search for synergy: a critical review from a response surface perspective*,
  Pharmacological Reviews — the universal-response-surface interaction parameter Onkos's
  `greco` model mirrors.
- **Combinations in oncology TGI/QSP.** Reviews of tumor-dynamics combination modeling —
  the layer one step above where Onkos's interaction-model-selection axis sits.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not synergy estimation.** Onkos does not fit `ψ` from data; it propagates a *declared*
  interaction assumption and shows the divergence. Estimating synergy needs a combination
  trial designed for it.
- **Not dose-level Loewe (yet).** v0.x combines at the effect level. Dose-level Loewe over
  the `er_emax` curves (solving the dose-equivalence equation) is a clean, separable
  extension, noted not hidden.
- **Not therapy ranking.** The module simulates one regimen under interaction assumptions;
  it never recommends a combination over a monotherapy or ranks regimens.
- **No three-plus-drug surface in v0.x.** Pairwise interaction only; higher-order
  combinations are a later, explicitly enumerated step.

---

## 9. Safety & scope (unchanged hard line)

- **Population / regimen level only.** Everything is a property of a *published TGI model
  under a stated combination assumption*. It is not an estimate of any person's tumor or
  survival, and adds no individual-level output.
- **No therapy ranking, no recommendation.** Combination-vs-monotherapy curves are
  descriptive simulation of what an interaction model predicts, never a choice between
  treatments.
- **Synergy is labeled an assumption.** A non-zero interaction parameter is a declared
  input carrying a warning, never a fitted or curated quantity; the headline is the
  *spread* across interaction assumptions, not a single synergistic answer.
- **The underlying tier governs and cannot be raised** by any interaction assumption.
- **The line, restated.** Any feature that takes a real patient's measurement and returns a
  combination outcome or a therapy choice **does not get built.** Making the interaction
  assumption legible and its consequence measurable changes none of this; it shows the
  trial-level disagreement, and stops there.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Interaction core** | `onkos.interaction.combine_effects` with `hsa` / `additive` / `greco(ψ)`; the closed-form landmark suite (§5). | The three nulls compute, the Bliss≡additive identity holds, and every landmark passes. |
| **2 — Simulation bridge** | `simulate_combination` feeding `E_AB` through the existing TGI → survival chain; tier passthrough; the synergy-is-an-assumption warning. | A combination trajectory simulates with the underlying model's tier and warnings, plus the assumption flag. |
| **3 — Divergence view & surfaces** | `compare_interactions` with the interaction-model OS/tumor divergence; `onkos interactions` CLI; `to_dict` with the clinical-use flag; a combination figure + a CI-executed notebook. | The interaction-model divergence computes and is visualized; the notebook runs in CI. |
| **4 — Documentation & honesty wiring** | README section framing the interaction model as a model-selection axis; the Bliss≡additive identity and the synergy-as-assumption discipline documented; roadmap + cheat sheets updated. | The combination axis is documented with the same rigor and the same honesty guardrails as the rest of the project. |

Step 1 alone is a self-contained, citable contribution: an open, validated tool that takes
a curated oncology TGI model and two single-agent activities and shows, with the standard
interaction nulls, how much a combination's predicted benefit depends on an unmeasured
synergy assumption — turning the field's qualitative *"assume additivity (or synergy)"*
into a measured divergence, is something nobody ships openly as a tested artifact.
