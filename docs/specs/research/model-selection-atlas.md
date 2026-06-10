# Onkos — research spec: the model-selection atlas — one index and one survey of every axis

**Status:** implemented in v0.39.0 (`onkos.atlas`). This is the design-of-record; written in the v0.1
house style. The atlas is a synthesis/navigation layer over the existing axes; no values are introduced
beyond what the underlying modules already compute.

**Eighteen versions turned each silent modeling choice into its own quantified axis, in its own module —
and the project lost the forest for the trees.** Onkos now makes explicit: which TGI model, which
resistance mechanism and origin, which bridge metric, which survival structure, which exposure-response
shape, which readout timing — and whether a trial could even tell the models apart. Each lives in a
separate module with its own headline and its own CLI. There is no single place that says, for a context,
*which* of these choices matters and by how much. This spec adds that synthesis layer: a declarative
**registry** of the model-selection axes (the source of truth for "what model-selection risk does Onkos
quantify") and a one-call **per-context survey** that runs each applicable axis and returns its native
headline.

> The atlas is deliberately a *navigational survey*, not a variance decomposition. The axes are not
> orthogonal; their headlines are in different units (weeks of median-OS spread, discordant model pairs,
> required trial events); and the magnitudes depend on the operating point. So the atlas reports each axis
> in its **own** unit, flags `comparable = False`, and points to `onkos.model_selection_budget` (v0.26)
> for the rigorous, orthogonal two-factor partition. The atlas answers "where is the risk, and which
> module do I open?"; the budget answers "what fraction of the variance is structural?"

---

## 1. The problem this synthesizes

The model-selection arc shipped these axes, each as a standalone module:

| Axis | varies | module | since |
| --- | --- | --- | --- |
| TGI model | which growth model | `compare` | v0.13 |
| survival link / bridge metric | which on-treatment metric drives the hazard | `simulate` (`survival_link=`) | v0.25 |
| survival structure | two-stage PH vs the joint current-value link | `joint` | v0.34 |
| exposure-response shape | which ER curve maps dose → effect | `dose_response` | v0.36 |
| readout timing | when the early surrogate is read | `early_surrogate` | v0.37 |
| model discriminability | whether a trial can tell the models apart | `discriminability` | v0.38 |
| additivity reference | which combination "no-interaction" null | `interaction` / `loewe` | v0.35 |

`onkos.atlas` adds (1) `AXES`, a frozen registry of `Axis(key, label, varies, finding, module, cli, scope,
version)` records — the canonical list the README cheat-sheet mirrors and the landmark suite checks for
completeness; and (2) `model_selection_atlas(ds, context)`, which runs each applicable **single-agent**
axis and returns an `Atlas` of `AtlasEntry(key, label, headline, unit, detail)`. Combination axes (the
additivity reference) need a regimen, not just a context, so they are catalogued in the registry but not
run by the per-context survey.

**Why this is the right synthesis (and the right scope).** It (1) is the navigation/index layer the
project needs after eighteen axis-versions — a discoverability and "where do I look" win; (2) is **pure
orchestration** over the existing modules — no new dataset record, kernel, schema, or export, so every
default artifact is byte-identical; (3) is the registry as **docs-from-code** — the README's axis
cheat-sheet now has a single source of truth, so docs cannot silently drift from the shipped modules;
(4) is *safe by construction* — it surfaces only what the underlying modules already compute, inherits
their tier, and adds no individual quantity; (5) is scrupulously honest about what it is **not** — it
refuses to claim cross-axis comparability and routes the rigorous question to the budget.

---

## 2. The result — one map of where the model-selection risk lies

For NSCLC first line, the per-context survey (illustrative):

| Axis | headline | unit |
| --- | --- | --- |
| **survival structure (PH vs joint)** | **~108** | weeks (median-OS spread) |
| **survival link / bridge metric** | **~97** | weeks (median-OS spread) |
| TGI model | ~41 | weeks (median-OS spread) |
| exposure-response shape | ~22 | weeks (median-OS spread; max on extrapolation) |
| early-surrogate readout timing | 8 / 10 | discordant model pairs at the earliest readout |
| model discriminability | 4 / 10 | model pairs needing an infeasible trial |

Read it as a map, not a ledger:

1. **The survival-side *structure* and *metric* choices dominate the OS swing** (~100+ weeks each), above
   the TGI-model choice (~41) and the exposure-response shape (~22, and only on dose extrapolation — it is
   zero at the studied dose). The four weeks-unit axes form a loosely-comparable leaderboard
   (`os_swing_axes`): "weeks of median OS riding on this one choice." They are *not* orthogonal and do not
   sum to a total — that is the budget's job — but as magnitudes they are directly interpretable and
   directly comparable.

2. **The detectability axes reframe the leaderboard.** The same survey reports that, at the ctDNA-era early
   readout, 8 of 10 model pairs are misranked relative to durable benefit, and 4 of 10 model pairs are
   practically indistinguishable by a realistic trial. So the largest OS swings sit beside the sobering
   fact that the model choice driving them often cannot be resolved by data — the project's whole thesis,
   in one view.

**The honest framing.** The atlas computes nothing new; it is a faithful, tested orchestration of the
existing modules (the landmark suite checks each headline equals its module's output for the same inputs).
Its only claims are (a) here are the axes, (b) here is each one's native headline, (c) the weeks-unit
group is loosely comparable, and (d) for the rigorous orthogonal partition, see the budget. It crowns no
axis and emits no recommendation.

---

## 3. Tier & guardrail propagation (unchanged invariants)

- **Inherits, never moves, the tier.** The survey carries the context's propagated tier (out-of-context
  models floor to D upstream); the atlas adds nothing that could raise it.
- **A survey, not a verdict.** `comparable = False` and the budget pointer are structural fields, not
  prose — the atlas cannot be mistaken for a decomposition.
- **Population / trial level only.** Every headline is a population/trial quantity from an underlying
  module; the atlas introduces no individual output.
- **Default view untouched.** No record, kernel, schema, or export changes; every default artifact is
  byte-identical.

---

## 4. Validation landmarks

`tests/test_atlas.py`:

| Landmark | Condition |
| --- | --- |
| **Registry well-formed** | every `Axis` has all fields, unique keys, a valid scope, a `v0.` version, and an importable module. |
| **Registry covers the shipped axes** | the seven model-selection axes are all registered (catches an un-registered new axis or a renamed module). |
| **Survey runs the single-agent axes** | the six single-agent axes appear; the combination axis (needs a regimen) does not. |
| **Headlines agree with the modules** | the survey's discriminability / readout-timing / TGI-model headlines equal the underlying modules' outputs for the same horizon — a faithful orchestration. |
| **OS-swing group is correct** | `os_swing_axes` is exactly the four weeks-unit axes, each a non-negative spread. |
| **Survey, not decomposition** | `comparable = False` and the note points to the budget. |
| **Tier & clinical-use** | the propagated tier rides through; the JSON carries the clinical-use prohibition and `comparable: false`. |
| **Reproduces across contexts** | the survey is well-formed for the other solid-tumor contexts. |

---

## 5. API, CLI, and surface

```python
from onkos.atlas import AXES, model_selection_atlas

AXES                       # the registry: every model-selection axis, its module, CLI, finding, version
a = model_selection_atlas(ds, context=ctx)
a.entries                  # per-axis native headline for the context
a.os_swing_axes            # the loosely-comparable weeks-unit leaderboard
a.comparable, a.note       # False; "see onkos.model_selection_budget for the rigorous partition"
```

```bash
onkos atlas --tumor-type NSCLC --line first
```

**No new module dependencies, record, kernel, or export** — `onkos.atlas` is a pure orchestration layer
surfaced through a CLI command, a figure, and a CI-executed notebook. This release also includes a
**housekeeping pass**: the architecture diagram's analyses box and notebook count, and the repository-
layout module list, were refreshed to include every module shipped since (a doc-drift fix).

---

## 6. Deliberate non-goals (so the scope stays honest)

- **Not a variance decomposition.** The atlas surveys; the budget partitions. The atlas refuses to sum or
  rank across units and says so in a structural field.
- **Not new science.** It introduces no new metric, model, or finding; it indexes and surfaces the
  existing ones. (That is the point — the synthesis was the gap.)
- **Not the combination axes (per context).** The additivity reference needs a regimen (two drugs, two
  doses), so it is registered but not run by the context survey.
- **Not auto-generated docs (yet).** The registry is the source of truth and the README mirrors it by
  hand, checked by the completeness landmark; emitting the README table directly from `AXES` is a clean
  follow-on.

---

## 7. Safety & scope (unchanged hard line)

- **Population / trial level only.** Every headline is a population/trial quantity; the atlas adds no
  individual output.
- **No recommendation.** It maps where the risk lies; it does not rank choices as "right" or recommend a
  model, dose, or trial.
- **Cannot raise a tier, computes nothing new.** It inherits the underlying modules' tiers and outputs.
- **The line, restated.** Any feature that takes a real patient's data and returns a prognosis or a
  therapy choice **does not get built.** Indexing the model-selection axes changes none of this.

---

## 8. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Registry** | `Axis` + `AXES`, the declarative source of truth for the model-selection axes. | the well-formed and coverage landmarks pass. |
| **2 — Survey** | `model_selection_atlas` running each applicable single-agent axis to its native headline (worst-case-over-models for the OS-swing axes, so a complete responder cannot mask a swing). | the headlines-agree-with-modules landmark passes. |
| **3 — Honesty guardrails** | `comparable = False`, the `os_swing_axes` group, and the budget pointer. | the survey-not-decomposition landmark holds; default view byte-identical. |
| **4 — Surfaces + housekeeping** | CLI `onkos atlas`, an OS-swing + detectability figure, a CI-executed notebook, README atlas section + cheat-sheet, and the architecture/layout doc-drift fix. | the atlas is visualized and documented, and the architecture docs match the shipped modules. |

Step 1–2 alone is a self-contained contribution: a tested registry and a faithful per-context survey of
the model-selection axes. Step 3 is what keeps it honest — the atlas is the navigation layer, the budget is
the decomposition, and the atlas says so in code, not just prose. Step 4 closes the doc drift that
eighteen rapid versions accumulated, so the architecture diagram and module list once again match the
repository.
