# Onkos — research spec: mechanistic (two-population) resistance

**Status:** implemented in v0.24.0 (`two_population_resistance` kernel + record; the
resistance mechanism as a model-selection axis). This is the design-of-record; the
methodological source anchors of §7 are documented but their Crossref-verified citation
curation is still pending, honest by design. Written in the v0.1 house style; every
parameter value is illustrative and `unverified` by design — the infrastructure is the
contribution.

**The Hydra term, made mechanistic.** Resistance is the central modeling challenge of
oncology (spec §1: the λ "Hydra" term — cut one head, two grow back). Onkos has modeled
it *phenomenologically*: the Claret model encodes resistance as an **exponential decay of
the drug effect** (`kD·E·e^{−λt}`), a curve-fitting device whose λ has no cellular
referent and is ~90%-CV unidentifiable (v0.22). This spec adds the *mechanistic*
alternative — a tumor of **sensitive + pre-existing resistant subclones**, where the drug
kills only the sensitive cells and the resistant clone outgrows — and makes the **choice
between the two resistance models a quantified model-selection axis**, the same move the
project already makes for the kill mechanism (log-kill vs Norton-Simon) and the
drug-interaction model.

> Two trials can show the same early tumor shrinkage and the same nadir, yet predict
> different long-term outcomes because one fit resistance as a fading drug effect and the
> other as an expanding resistant subpopulation. The phenomenological and mechanistic
> models are observationally close on a short trial and divergent on the long-horizon tumor
> burden that should gate go/no-go — exactly the silent, outcome-determining modeling choice Onkos
> exists to surface. And the mechanistic form gives resistance a *biologically
> interpretable* parameter (the initial resistant burden `R0`) in place of an
> unidentifiable rate.

---

## 1. The problem this extends

| Resistance representation | Form | Parameter meaning | Status |
| --- | --- | --- | --- |
| **Phenomenological** (Claret) | kill term decays: `dV/dt = (kL − kD·E·e^{−λt})·V` | `λ` = rate the drug effect fades — *no cellular referent*; ~90% CV, unidentifiable (v0.22). | ✅ since v0.1 |
| **Mechanistic** (two-population) | sensitive clone killed, resistant clone grows: `V = S + R` | `R0` = initial resistant burden (a *biological* quantity); the regrowth is the resistant clone outgrowing. | ⚠️ **this spec** |

The phenomenological model is a fit; the mechanistic model is a hypothesis about *why*
the tumor regrows. They are not interchangeable: the same nadir can be produced by a
fading effect or by a small resistant clone, and **which you assume determines the
predicted tail** — the part of the curve that drives overall survival. Onkos already
proves resistance is poorly identifiable from short trials (v0.22); this spec shows that
the *structural* choice of resistance model is a second, compounding uncertainty on top of
the parameter one.

**Why this is the right deepening (and the right scope).** It (1) advances the project's
own thesis on its most load-bearing term — resistance — by making the *mechanism* an
explicit divergence axis rather than breadth; (2) reuses the existing multi-state ODE
machinery (the Simeoni transit chain and the IO QSP are already two-plus-state kernels
with round-trip export), so it is a new *kernel*, not new infrastructure; (3) has direct
precedent — the pre-existing-resistant-subclone model is the canonical Goldie-Coldman
formulation (§7); (4) is *safe by construction* — population/trial-level forward
simulation, no individual prognosis, no therapy ranking; and (5) sharpens the honest
message: it replaces an unidentifiable phenomenological rate with an interpretable
biological one *and* shows that even so, the resistance-model choice moves the survival
forecast.

---

## 2. The mechanistic model

A solid tumor of size `V` is two clones — a drug-**sensitive** population `S` and a
drug-**resistant** population `R`, observed together as `V = S + R`:

```
dS/dt = (kg − kd·E)·S          sensitive: net growth kg, killed at potency kd by effect E
dR/dt =  kgr·R                 resistant: grows at kgr, NOT killed by the drug
S(0) = V0 ,  R(0) = R0         a small pre-existing resistant burden R0 (Goldie-Coldman)
```

The kill potency `kd` matches the Claret parameterization (`kd·E` is the per-time kill on
sensitive cells), so when the two models are compared the *effective single-agent kill is
identical* and the divergence is purely the **resistance mechanism**, not an effect-scale
artifact. The system is linear and decoupled, so it is closed-form
(`S(t) = V0·e^{(kg−kd·E)t}`, `R(t) = R0·e^{kgr·t}`); Onkos integrates it numerically through
the same multi-state path as Simeoni, and exports both compartments to SBML/NONMEM with
the standard round-trip validation.

**Behavioral signature (the landmarks, §5).** With `R0 = 0` the model is pure sensitive
exponential kill — under an effective kill `kd·E > kg` the tumor is *eradicated*
(`V → 0`). With `R0 > 0` and `kd·E > kg`, the sensitive clone collapses, the tumor reaches
a **nadir**, and then **regrows driven entirely by the resistant clone** — so the
asymptotic log-growth rate is `kgr` (the resistant clone dominates), and the resistant
*fraction* `R/V` rises monotonically under treatment from `R0/V0` toward 1. This is the
mechanistic origin of the nadir-then-regrowth that the phenomenological λ approximates by
hand.

---

## 3. The divergence view — the resistance model as a model-selection axis

Because the mechanistic record has `purpose = "tgi"` in the NSCLC first-line context where
the Claret model lives, it joins the **virtual-trial divergence view automatically**: the
overlay now spans *two resistance mechanisms* (phenomenological decay-of-effect vs
mechanistic resistant-subclone) alongside the biexponential and Norton-Simon models. The
headline is the same as the kill-mechanism and interaction axes: *of everything you might
predict for this tumor, how much is just the resistance model you assumed* — now made
legible for the project's most load-bearing term.

**The sharp, honest finding — and where it hides.** The two models are tuned to share the
early kill (matched `kd`), so they agree closely at week 8 (≈−87% vs −82% change) and hence
on the **week-8-driven OS** (median ≈94 vs ≈91 wk). Yet they diverge **≈5× in the tumor
tail** (≈74 vs ≈15 mm at 3 years), because one regrowth is a fading effect and the other a
compounding resistant clone. This is precisely the short-trial-indistinguishable,
long-horizon-divergent failure mode the project exists to expose — and it carries a second,
uncomfortable lesson: a **week-8-based OS surrogate is nearly blind to the resistance-model
choice**, which is exactly how a short-trial-fit resistance model transports silently into a
late-phase prediction it cannot support. The divergence axis is real; the standard surrogate
does not see it. (An endpoint sensitive to the regrowth tail — e.g. a late landmark or a
TTP-style metric — would surface the full split; building one is a clean follow-on.)

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins still governs.** The mechanistic record is tier C (single illustrative
  context, no external validation, an unidentifiable resistant fraction); a composed
  forecast through it inherits that, and an out-of-context transport floors to D exactly as
  for every other model.
- **The interpretable parameter is still honestly uncertain.** `R0` (initial resistant
  burden) is biologically meaningful but, like λ, poorly identified from a short trial; it
  carries a high IIV CV and is surfaced by the identifiability analyzer (v0.22) — mechanistic
  does not mean measured.
- **Population/trial level only.** `S`, `R`, and `V` are trial-level trajectories of a
  published model, never a patient's clones; no individual-level output is added, no therapy
  is ranked.

---

## 5. Reference kernel & validation landmarks

A new multi-state ODE kernel `two_population_resistance`, validated on two independent axes
exactly like every kernel: the **round-trip** (SBML MathML and NONMEM `$THETA` re-parsed
and checked against the reference rhs, per compartment) and a **landmark suite** of the
characteristic, analytically-derivable properties of the *published model* it implements.

| Landmark | Closed-form condition |
| --- | --- |
| **Closed form** | `V(t) = V0·e^{(kg−kd·E)t} + R0·e^{kgr·t}` matches the numerically integrated observable. |
| **No-resistance reduction** | `R0 = 0` ⇒ `V(t) = V0·e^{(kg−kd·E)t}` (pure sensitive exponential kill). |
| **Eradication** | `R0 = 0` and `kd·E > kg` ⇒ `V → 0` monotonically (no regrowth). |
| **Resistant dominance** | `R0 > 0`, `kd·E > kg` ⇒ the late-time log-slope of `V` → `kgr` (the resistant clone sets the tail). |
| **Nadir then regrowth** | `R0 > 0`, `kd·E > kg`, `kgr > 0` ⇒ `V` has an interior minimum (a nadir), then rises. |
| **Resistant-fraction monotonicity** | `R/V` increases monotonically under treatment from `R0/V0`. |
| **Untreated reduction** | `E = 0` ⇒ both clones grow; `V(t) = V0·e^{kg·t} + R0·e^{kgr·t}`. |
| **Round-trip per compartment** | exported SBML/NONMEM rate laws re-evaluate to the reference rhs for both `S` and `R` (≤1e-6). |

This mirrors the project's two-axis validation discipline: round-trip proves exports match
the kernel; landmarks prove the kernel *is* the Goldie-Coldman two-population model it
names, not an unconstrained two-exponential fit.

---

## 6. API, CLI, and surface

No new module — the mechanistic model is a dataset record bound to a new kernel, so it flows
through the *existing* surfaces:

```python
# Forward simulation of the mechanistic model (sensitive + resistant clones).
tr = onkos.simulate(ds, "resistance.nsclc_first_line.two_population", context=ctx, drug_effect=1.0)
tr.tumor_size, tr.os_curve, tr.tier        # V = S + R; population OS; propagated tier

# It joins the virtual-trial divergence view automatically — now spanning two
# resistance mechanisms (phenomenological decay-of-effect vs mechanistic subclone).
cmp = onkos.compare(ds, purpose="tgi", context=ctx, drug_effect=1.0)
cmp.os_divergence                          # absorbs the resistance-model choice

# Exports: both compartments to SBML/NONMEM/rxode2/Pumas, round-trip validated.
onkos.identifiability(ds, "resistance.nsclc_first_line.two_population", context=ctx)  # R0 is poorly identified
```

```bash
onkos simulate resistance.nsclc_first_line.two_population
onkos simulate --compare        # the divergence view now includes the mechanistic model
onkos export --format sbml --output exports/sbml/   # two species: sensitive, resistant
```

---

## 7. Source anchors (methodological; DOIs added at curation time)

Well-established methods, not Onkos parameters; added to `dataset/citations/` through the
normal Crossref/PubMed-verified pipeline, honest until a human confirms each.

- **Pre-existing resistant subclones.** Goldie, J.H. & Coldman, A.J. (1979), *A mathematical
  model for relating the drug sensitivity of tumors to their spontaneous mutation rate*,
  Cancer Treatment Reports — the canonical sensitive/resistant two-population formulation
  this kernel implements (a pre-DOI reference; PMID anchor).
- **Resistance dynamics under therapy.** Foo, J. & Michor, F. (2014), *Evolution of acquired
  resistance to anticancer therapy*, Journal of Theoretical Biology — the modern treatment
  of resistant-subpopulation outgrowth under selective drug pressure.
- **Phenomenological contrast.** Claret et al. (2009) — the decay-of-drug-effect model
  (already in the dataset) that the mechanistic form is compared against.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not acquired-resistance switching.** v0.x models *pre-existing* resistance (a fixed
  initial resistant burden, Goldie-Coldman). An explicit sensitive→resistant transition rate
  (acquired resistance / phenotype switching) is a clean, separable extension, noted not
  hidden.
- **Not partial resistance.** The resistant clone is taken fully resistant (`kill = 0`); a
  resistant-but-attenuated kill is a later refinement.
- **Not clonal-evolution multi-state.** Two clones only; `n`-clone / mutation-network models
  are a different, larger subsystem.
- **No estimation of the resistant fraction.** `R0` is a declared illustrative parameter; the
  identifiability analyzer (v0.22) shows a short trial cannot pin it down.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only.** `S`, `R`, `V`, and the OS curve are trajectories of a
  published model for a trial-level context, never an estimate of any person's tumor, clones,
  or survival. No individual-level output is added.
- **No therapy ranking, no recommendation.** The mechanistic model is one more model in the
  divergence view; it ranks nothing.
- **Mechanistic does not mean measured.** The interpretable `R0` is still tier-C and
  unidentifiable from a short trial; the honesty fields say so.
- **The line, restated.** Any feature that takes a real patient's measurement and returns a
  clonal composition, a prognosis, or a therapy choice **does not get built.** Making the
  resistance mechanism an explicit, comparable modeling choice changes none of this.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Kernel** | `two_population_resistance` (sensitive + resistant clones, kill on sensitive only); rhs + rhs_infix + observable; multi-state seeding via `init_inputs`. | The kernel integrates, the closed form matches, and the round-trip passes for both compartments. |
| **2 — Record + landmarks** | An NSCLC first-line mechanistic record (kd matched to Claret) + the Goldie-Coldman citation; the landmark suite (§5). | The record validates, joins the divergence view, and every landmark passes. |
| **3 — Divergence axis & presentation** | The resistance-model divergence (phenomenological vs mechanistic) made visible; a figure + a CI-executed notebook; refreshed virtual-trial numbers. | The two resistance mechanisms overlay in the divergence view with a measured OS divergence. |
| **4 — Documentation** | README section framing the resistance mechanism as a model-selection axis; roadmap, layout, and the perturbed NSCLC illustration numbers updated. | The mechanistic-resistance axis is documented with the same rigor and honesty guardrails as the rest of the project. |

Step 1 alone is a self-contained, citable contribution: an open, round-trip-validated,
landmark-tested implementation of the Goldie-Coldman two-population resistance model that
replaces an unidentifiable phenomenological rate with an interpretable biological one — and,
placed beside the phenomenological model, shows that the resistance-mechanism choice is its
own measurable model-selection risk on the survival forecast.
