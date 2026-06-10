# Onkos — research spec: early-surrogate readout timing — *when* you read the surrogate is a model-selection axis

**Status:** implemented in v0.37.0 (`onkos.early_surrogate`). This is the design-of-record; written in
the v0.1 house style. All values are illustrative and `unverified` by design; the infrastructure is the
contribution.

**The metric-choice work asked *which* on-treatment quantity predicts survival; the ctDNA era forces the
orthogonal question — *when* do you read it?** The field is pushing the surrogate readout ever earlier:
circulating-tumor-DNA (ctDNA) "molecular response" at week 2-4, well before a RECIST tumor-size change is
reliable at week 8, on the premise that an earlier signal is an earlier go/no-go. This spec makes the
**readout landmark time** an explicit model-selection axis, distinct from the metric (v0.25/v0.33) and
survival-structure (v0.34) axes, and shows that **earliness trades against fidelity**: the earlier you
read, the more the surrogate over-rewards deep-but-doomed early responders.

> Onkos models ctDNA molecular response as proportional to tumor burden — the standard first-order
> shedding assumption (ctDNA ∝ tumor volume). The modeled distinction between a ctDNA readout and a
> RECIST-size readout is therefore purely the **landmark time** (and assay noise). Genomic / mutational
> ctDNA content is out of scope (spec §2: Onkos stays at the tumor-dynamics scale). That deliberate
> reduction is what isolates the timing question cleanly — it is the honest scope, not a limitation
> hidden.

---

## 1. The problem this extends

Every surrogate Onkos has shipped reads the trajectory at a **fixed** landmark — week 8 (v0.12, v0.25,
v0.33). But the landmark is a choice, and the ctDNA program's entire value proposition is *moving it
earlier*. `landmark_response(t, v, week)` generalizes the fixed week-8 covariate to an arbitrary landmark
— the relative tumor-burden change `(v(week) − v0)/v0` — which, under the proportional-shedding
assumption, is also the modeled ctDNA molecular response at that week. It recovers
`week8_relative_change` exactly at `week = 8`.

The timing axis is then measured against a **durable-benefit reference**: the ranking of a context's
models by median OS under the tail-aware k_g survival link (which rewards slow regrowth, v0.25/v0.29).
`surrogate_timing_fidelity` ranks the models by their early-surrogate response at each landmark week and
counts how many model pairs that ranking orders *oppositely* to the durable-benefit ranking
(a Kendall-style discordance).

**Why this is the right deepening (and the right scope).** It (1) is the reserve list's "ctDNA
early-surrogate" direction, scoped honestly as the timing question it really is; (2) is **pure
post-processing** over the existing trajectories and the existing k_g link — a module with a
landmark-tested pure core, no new dataset record, kernel, schema, or export, so every default artifact is
byte-identical; (3) adds the **time** dimension to the bridge-metric story (v0.25 = which metric, v0.33 =
which integrated quantity, v0.37 = at which time), orthogonal to both; (4) is *safe by construction* — a
population ranking over published trajectories, no individual molecular response, no go/no-go
recommendation; (5) reinforces the depth-vs-durability thread (v0.27/v0.28/v0.34): the earlier the
landmark, the more a transient deep response masquerades as durable benefit.

---

## 2. The result — earliness trades against fidelity

For NSCLC first line, the durable-benefit (k_g link OS) reference ranks the models
`Norton-Simon > Wang > Claret > two-population > acquired` (the complete responder first; the
fast-deep-then-fast-regrow mechanistic-resistance models last). The discordance of the early-surrogate
ranking against that reference, by landmark week (illustrative):

| Landmark week | discordant pairs / 10 | top-ranked by the surrogate |
| --- | --- | --- |
| **2** (ctDNA era) | **9** | acquired (deep, doomed) |
| 4 | 8 | two-population (deep, doomed) |
| 8 (RECIST) | 8 | two-population |
| 12 | 7 | two-population |
| 16–24 | 7 | two-population |
| 36 | 5 | Claret |
| 52 | **3** | Norton-Simon (durable) |

Two findings:

1. **Discordance falls monotonically as the landmark moves later.** At the ctDNA-era week-2 readout the
   surrogate ranking is almost completely inverted relative to durable benefit (9/10 pairs wrong); by week
   52 it is mostly aligned (3/10). The earlier you read, the worse the agreement with the outcome that
   matters — because an early landmark sits at or before the nadir, before any resistant regrowth, so it
   cannot see the tail that determines durable benefit.

2. **The bias has a direction: it over-rewards deep-but-doomed early responders.** At the earliest
   landmarks the surrogate's top model is a mechanistic-resistance model (acquired / two-population) —
   precisely the models the durable-benefit reference ranks *last*. Their fast, deep early response (the
   ctDNA molecular response a week-2-4 readout would celebrate) is exactly what the regrowth tail will
   undo. This is the depth-vs-durability surrogate failure (v0.27/v0.28), localized in *time*: shifting
   the readout earlier maximizes it.

**Cross-context.** The earliness-fidelity trade-off reproduces in the other solid-tumor contexts (breast,
CRC, HCC, melanoma), each of which has had a two-population model and a k_g link since v0.29 — it is not
an NSCLC artifact.

**The honest framing.** Onkos does not claim an early readout is wrong, only that it carries a
*timing-specific model-selection bias* that a single early go/no-go hides. The reference is one tail-aware
choice (the k_g link); the headline is the monotone *trade-off*, not an optimal landmark. ctDNA's modeled
advantage here is timing alone — the deliberate scope that lets the timing axis be measured cleanly.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** The rankings ride the trajectories' propagated tier; out-of-context
  models floor to D upstream in `compare`. A timing analysis cannot raise a tier.
- **A timing choice, not a fit.** The landmark grid and the reference link are declared; nothing is
  estimated, and no landmark is endorsed as correct.
- **Population / trial level only.** A ranking of *models under a context* across readout times; never an
  individual molecular response, prognosis, or go/no-go recommendation.
- **Default view untouched.** No record, kernel, schema, or export changes; every default artifact is
  byte-identical.

---

## 4. Validation landmarks

The pure core is closed-form; the binding adds the finding (`tests/test_early_surrogate.py`):

| Landmark | Condition |
| --- | --- |
| **Recovers the week-8 metric** | `landmark_response` at week 8 equals `week8_relative_change` to `1e-9`. |
| **Record-free, signed** | works on a synthetic trajectory; 0 at baseline, negative under shrinkage. |
| **Deepens then recovers** | for a shrink-then-regrow trajectory the response is deeper at the nadir landmark than at a late one — why an early readout flatters a doomed responder. |
| **Discordance counts inversions** | `discordant_pairs` is 0 for identical, `n(n−1)/2` for fully reversed rankings. |
| **Fidelity improves with a later landmark** | the earliest-landmark discordance exceeds the latest — the core trade-off. |
| **Early readout is badly discordant** | the earliest landmark's discordance is ≥ 70% of all pairs. |
| **Weakly monotone** | discordance never increases as the landmark moves later. |
| **Reference demotes the fast-doomed models** | the durable-benefit ranking puts the mechanistic-resistance models near the bottom, while the early-landmark ranking puts one on top. |
| **Reproduces across contexts** | the trade-off holds for breast, CRC, HCC, melanoma. |
| **Tier & clinical-use** | the propagated tier rides through; the result carries the clinical-use prohibition. |

---

## 5. API, CLI, and surface

```python
from onkos.early_surrogate import landmark_response, surrogate_timing_fidelity

landmark_response(tr.t, tr.tumor_size, 4.0)   # the ctDNA-era (week-4) molecular-response readout

st = surrogate_timing_fidelity(ds, context=ctx)   # discordance vs durable benefit, by landmark week
st.discordance_at(2.0), st.discordance_at(52.0)    # earliest vs latest — earliness trades against fidelity
st.reference_ranking                                # the tail-aware durable-benefit ranking
```

```bash
onkos early-surrogate --tumor-type NSCLC --line first
```

**No new module dependencies, record, kernel, or export** — `onkos.early_surrogate` is a pure
post-processing module surfaced through a CLI command, a figure, and a CI-executed notebook.

---

## 6. Source anchors (methodological; values illustrative)

- **ctDNA molecular response as an early endpoint.** The circulating-tumor-DNA literature on early
  molecular response (week 2-4) as a candidate surrogate that precedes radiographic response — the
  motivation for the timing push this axis quantifies.
- **ctDNA ∝ tumor burden.** The first-order shedding assumption (ctDNA concentration proportional to
  tumor volume) under which the modeled ctDNA-vs-RECIST distinction reduces to readout timing.
- **Tumor growth-rate constant as durable-benefit reference.** Stein (2008)
  (DOI 10.1634/theoncologist.2008-0075) — the tail-aware k_g link used as the durable-benefit ranking.

---

## 7. Deliberate non-goals (so the scope stays honest)

- **Not genomic ctDNA.** Mutational / clonal ctDNA content is out of scope (§2); the model is burden-
  proportional, so the axis is timing, stated plainly.
- **Not an optimal landmark.** The module quantifies the trade-off; it does not recommend a readout week.
- **Not a refit, not assay noise.** The landmark grid and reference are declared; measurement noise is not
  modeled (a clean follow-on would add a noise model to study earliness-vs-precision as well as
  earliness-vs-fidelity).
- **Not a new endpoint record.** It reuses the existing trajectories and the k_g link; per-context ctDNA
  shedding records are a breadth follow-on.

---

## 8. Safety & scope (unchanged hard line)

- **Population / trial level only.** A model ranking over readout times; nothing is an individual
  molecular response or prognosis.
- **No go/no-go recommendation.** The axis ranks *models under a context*, never a trial decision or a
  landmark.
- **Cannot raise a tier, fits nothing.** The landmark and reference are declared; the underlying models'
  tier governs.
- **The line, restated.** Any feature that takes a real patient's ctDNA and returns a molecular-response
  call, a prognosis, or a go/no-go **does not get built.** Making the readout-timing bias explicit and
  computable changes none of this.

---

## 9. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Readout core** | `landmark_response(t, v, week)` generalizing the week-8 covariate; `discordant_pairs` ranking distance. | the week-8 recovery and the deepens-then-recovers landmarks pass. |
| **2 — The timing axis** | `surrogate_timing_fidelity` — discordance of the early-landmark ranking vs the durable-benefit (k_g) ranking, over a landmark grid. | the §2 table holds; discordance is weakly monotone in the landmark. |
| **3 — The finding** | earliness-trades-against-fidelity and the over-rewards-deep-but-doomed direction, reproduced cross-context. | the fidelity-improves-with-later-landmark and cross-context landmarks hold; default view byte-identical. |
| **4 — Surfaces** | CLI `onkos early-surrogate`, a trajectories + discordance-curve figure, a CI-executed notebook, README + changelog + API contract. | the readout-timing axis is visualized and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested generalization of the surrogate
readout to an arbitrary landmark, with the week-8 recovery as its correctness anchor. Step 3 is the
payload: showing that the earlier you read the surrogate the more its model ranking contradicts durable
benefit — that the ctDNA-driven push to week-2-4 endpoints maximizes the depth-vs-durability bias — is the
quantitative core of this spec, shipped as a tested artifact.
