# Onkos — research spec: dose-level Loewe additivity — the additivity *reference* as a model-selection axis

**Status:** implemented in v0.35.0 (`onkos.interaction` extension: `ERCurve`, `loewe_effect`,
`combine_doses`, `compare_additivity_references`). This is the design-of-record; written in the v0.1
house style. All values are illustrative and `unverified` by design; the infrastructure is the
contribution.

**There is no single "no-interaction" reference for a drug combination — and the choice changes the
predicted benefit.** v0.23 made the *interaction model* a model-selection axis, but it operated at the
**effect level**: combine two single-agent effect magnitudes `E_A, E_B` under HSA (max), Bliss/effect-
additivity (`E_A + E_B`), or Greco. The classical gold-standard null it deferred — flagged in the v0.23
notes as "dose-level Loewe over the ER curves" — operates at the **dose level**: Loewe additivity
combines two *doses* through each drug's dose-response curve, via the isobole
`d_A/D_A(E) + d_B/D_B(E) = 1`. This spec adds it and shows that *which additivity reference you call
"no interaction"* is itself a model-selection axis, with a survival consequence, and that Loewe is the
only one of the three that is self-consistent.

> Loewe vs Bliss is the central, decades-old controversy in synergy quantification. They answer
> different questions — Loewe: "is the combination better than giving more of one drug?" (dose
> additivity, the right null for same-mechanism drugs); Bliss: "are the drugs acting independently?"
> (effect/probabilistic additivity). They coincide only for special curve shapes. Onkos does not
> adjudicate the controversy — it makes both references computable over the same doses and the same
> TGI→survival chain, and shows the spread.

---

## 1. The problem this extends

`onkos.interaction` (v0.23) combines effect magnitudes:

```
E_AB = combine_effects(E_A, E_B, model=...)   # hsa: max ; additive(Bliss): E_A+E_B ; greco: +psi·√(E_A E_B)
```

That is fine when you already have the two single-agent *effects*. But the effect a dose produces is a
point on a dose-response (exposure-response) curve, and the **dose-additive** null asks a different
question: holding the combination at doses `(d_A, d_B)`, what effect would dose-additivity predict? The
Loewe answer is the `E` on the isobole

```
d_A / D_A(E)  +  d_B / D_B(E)  =  1
```

where `D_x(E)` is the dose of drug `x` *alone* producing effect `E` — the **inverse** of its ER curve.
Onkos already curates the ER curves (`er_emax`, `er_sigmoid_emax`, `er_power`), all analytically
invertible, so the isobole is solved exactly. The combined effect then drives the *existing* TGI→
survival chain, identically to the v0.23 effect-level path.

**The defining property — the sham-combination identity.** Combine a drug with *itself*, splitting a
dose `d = d_A + d_B`. Loewe must return the single-agent effect at the summed dose: the isobole
collapses to `(d_A + d_B)/D(E) = 1 ⟹ D(E) = d_A + d_B ⟹ E = f(d_A + d_B)`. **Loewe satisfies this
exactly; Bliss does not** — for any saturating curve `f(d_A) + f(d_B) > f(d_A + d_B)`. That single
identity is why Loewe is the principled dose-additive reference and why the choice between it and Bliss
is not cosmetic.

**Why this is the right deepening (and the right scope).** It (1) is the pre-committed v0.23 follow-on
("dose-level Loewe over the er_emax curves"); (2) is **pure post-processing** — an extension of the
existing interaction module reusing the curated ER kernels, with a landmark-tested pure core
(`loewe_effect` over an `ERCurve`), no new dataset record, kernel, schema, or export; (3) ties the
exposure-response subsystem (until now only a transform inside `simulate`) into the model-selection
story for the first time; (4) is *safe by construction* — a population regimen simulation, the
reference a declared choice, never an estimated synergy; (5) reinforces the load-bearing message one
layer further out: even the *null* against which interaction is judged is a silent modeling choice.

---

## 2. The result — three references, one dose pair, three survival curves

For Claret NSCLC first line, doses `d_A = 150` (via `emax_generic`, `E = 1.4·C/(150+C)`) and
`d_B = 90` (via the dacomitinib Emax, `E = 1.8·C/(90+C)`), illustrative:

| Reference | combined effect `E_AB` | median OS | satisfies sham test? |
| --- | --- | --- | --- |
| HSA (`max`) | 0.90 | 88 | n/a (conservative bound) |
| **Loewe** (dose-additive isobole) | **1.07** | **92** | **yes** (exact) |
| Bliss (effect-additive `E_A+E_B`) | 1.60 | 101 | **no** (overstates) |

Two findings:

1. **The reference is a model-selection axis with an OS consequence.** From the *same* doses, the three
   "no-interaction" nulls give combined effects spanning 0.90–1.60 and median OS spanning 88–101 weeks.
   The ordering is structural for saturating curves: `HSA ≤ Loewe ≤ Bliss`. Effect-additivity (Bliss)
   *overstates* the combined effect — it can even exceed either drug's maximum (`E_AB = 1.60` above the
   shared ceiling `1.4`), which is physically impossible for a single saturating target and is the
   classic critique of Bliss for same-mechanism drugs. HSA *understates* it. Loewe is the
   self-consistent middle.

2. **The disagreement grows with dose.** At low doses (the steep, near-linear part of the curves) the
   three references nearly agree; as doses climb into saturation they diverge, because effect-additivity
   keeps adding effects that dose-additivity knows are redundant (you are climbing the same saturating
   curve). So the reference choice is negligible for a low-dose combination and large for a high-dose
   one — a dose-dependent model-selection risk, exactly where combination dose-finding lives.

**The honest framing.** Onkos estimates no synergy and crowns no reference. All three are illustrative,
tier-C; the reference is a *declared* choice (like `psi` in v0.23), and the headline is the spread.
Loewe's distinction is not that it is "right" but that it is the only one of the three that is
internally consistent (a drug with itself is additive) — a property the dataset makes checkable, not
asserted.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** The combined effect drives the underlying TGI→survival chain; the
  comparison carries that chain's propagated tier (C), and out-of-context transport floors to D + warns.
  The reference choice cannot raise a tier.
- **The reference is an assumption, not a measurement.** Like the v0.23 interaction model, the
  additivity reference is declared; nothing is fitted, no synergy is estimated.
- **Population / regimen level only.** Combination OS over published model structures; never an
  individual prediction, never a therapy or dose recommendation. The comparison ranks *references under
  a regimen*, not treatments.
- **Default view untouched.** No record, kernel, schema, or export changes; every default artifact is
  byte-identical.

---

## 4. Validation landmarks

The isobole core (`loewe_effect` over an `ERCurve`) is landmark-tested in isolation; the binding adds
the finding (`tests/test_loewe.py`):

| Landmark | Condition |
| --- | --- |
| **Sham-combination identity** | Loewe of a drug with itself equals `f(d_A + d_B)` to `1e-6` — the defining dose-additivity property. |
| **Bliss fails the sham test** | for a saturating curve `f(d_A) + f(d_B) > f(d_A + d_B)`, so Bliss ≠ Loewe; the gap is the reason the choice matters. |
| **Single-agent limits** | one dose zero returns the other drug's single-agent effect exactly. |
| **ER inverse round-trips** | `D_x(f_x(d)) = d` for the Emax, sigmoid-Emax, and power kernels (the analytic inverse is exact). |
| **Effect-ceiling clamp** | two saturating curves cannot jointly express more than `min(emax_A, emax_B)`; a huge dose pair clamps there, finite. |
| **Monotone in dose** | scaling both doses up raises the Loewe combined effect. |
| **Record-free core** | the isobole solves on plain callables (an `ERCurve` from lambdas), so it is testable without the dataset. |
| **Reference ordering** | `HSA ≤ Loewe ≤ Bliss` for the saturating curves, with a genuine spread. |
| **OS divergence** | the same dose pair gives different OS across references (`bliss > loewe > hsa`), the survival consequence. |
| **Tier & guardrails** | the comparison rides the chain's tier; out-of-context transport floors to D + warns; unknown reference/kernel raise. |

---

## 5. API, CLI, and surface

```python
from onkos.interaction import er_curve, loewe_effect, combine_doses, compare_additivity_references

ca = er_curve(ds, "exposure_response.emax_generic")        # forward f(d), inverse D(E), emax
loewe_effect(150, 90, curve_a=ca, curve_b=er_curve(ds, "exposure_response.dacomitinib_egfr.emax"))

# Combined effect under any reference, and the OS divergence across all three.
combine_doses(ds, 150, 90, er_a=..., er_b=..., reference="loewe")   # hsa | bliss | loewe
cmp = compare_additivity_references(ds, "resistance.claret_2009.tgi", context=ctx,
                                    dose_a=150, dose_b=90, er_a=..., er_b=...)
cmp.os_divergence, cmp.median_os    # the reference axis as a survival spread
```

```bash
onkos loewe resistance.claret_2009.tgi --dose-a 150 --dose-b 90 \
            --er-a exposure_response.emax_generic --er-b exposure_response.dacomitinib_egfr.emax
```

**No new module, record, kernel, or export** — an extension of `onkos.interaction`, surfaced through a
CLI command, a figure, and a CI-executed notebook.

---

## 6. Source anchors (methodological; values illustrative)

- **Loewe additivity / isobole.** Loewe, S. (1953), *The problem of synergism and antagonism of combined
  drugs* — the dose-additive null and the isobologram (the foundational reference for the
  sham-combination identity).
- **Bliss independence.** Bliss, C.I. (1939), *The toxicity of poisons applied jointly* — the
  effect/probabilistic-additive null contrasted here.
- **Synergy-reference controversy.** The modern reviews of Loewe-vs-Bliss (and the Greco interaction
  index already in v0.23) — context for treating the reference as a model-selection axis rather than a
  settled convention.
- **Exposure-response curves.** Holford & Sheiner Emax framework (DOI 10.2165/00003088-198106060-00002),
  the curated `er_emax` / `er_sigmoid_emax` / `er_power` kernels the isobole inverts.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not an estimated synergy.** Loewe here is the additive *null*, not a fitted interaction; a non-null
  interaction index (Greco `psi`) remains the v0.23 declared-assumption surface.
- **Not a dose recommender.** The comparison ranks references under a fixed dose pair, never doses or
  treatments; the dose-scaling sweep is descriptive, not an optimization.
- **Not a full response-surface.** The shipped surface is the combined effect at given doses under each
  reference; a full 2-D interaction surface / synergy-score map is a clean follow-on.
- **Not new ER records.** It uses the curated ER kernels as-is; calibrating context-specific combination
  ER curves is a breadth follow-on.

---

## 8. Safety & scope (unchanged hard line)

- **Population / regimen level only.** Combination OS over published structures; nothing is an individual
  prediction.
- **No therapy or dose ranking.** The reference comparison ranks *nulls*, never treatments or doses.
- **Cannot raise a tier, estimates nothing.** The reference is declared; the underlying model's tier
  governs.
- **The line, restated.** Any feature that takes a real patient's data and returns a combination
  regimen, a dose, or a prognosis **does not get built.** Making the additivity reference explicit and
  its OS consequence computable changes none of this.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — ER curves + inverses** | `ERCurve` + `er_curve` with analytic inverses for the Emax / sigmoid-Emax / power kernels. | the inverse round-trips to `1e-9`. |
| **2 — Isobole core** | `loewe_effect` solving the dose-additive isobole (bisection, no new dependency). | the sham-combination identity holds exactly; single-agent and ceiling-clamp limits pass. |
| **3 — The reference axis** | `combine_doses` / `compare_additivity_references` — HSA / Bliss / Loewe combined effect and OS over the existing chain. | the §2 table and the `bliss > loewe > hsa` OS ordering hold; default view byte-identical. |
| **4 — Surfaces** | CLI `onkos loewe`, a dose-scaling + OS figure, a CI-executed notebook, README + changelog + API contract. | the reference axis is visualized and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested Loewe isobole solver over the
curated ER curves, with the sham-combination identity as its correctness anchor. Step 3 is the payload:
showing that the "no-interaction" reference is a model-selection axis — that Bliss overstates, HSA
understates, and only Loewe is self-consistent — and that the choice moves the survival prediction, is
the quantitative core of this spec, shipped as a tested artifact.
