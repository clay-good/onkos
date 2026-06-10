# Onkos — research spec: duration of response (depth is not durability)

**Status:** implemented in v0.28.0 (`onkos.response` — DoR over the ensemble). This is the
design-of-record; the methodological source anchors of §7 are documented but their
Crossref-verified citation curation is still pending, honest by design. Written in the v0.1
house style; every value is illustrative and `unverified` by design — the infrastructure is
the contribution.

**The durability ORR cannot see.** v0.27 added the objective response rate (ORR) and showed
it is a *conditional* OS surrogate — faithful when survival tracks early shrinkage, inverted
when it tracks the regrowth tail. This spec supplies the missing dimension and the
*mechanism* of that failure: **duration of response (DoR)** — how *long* responses last.
ORR measures response **breadth** (how many patients respond); DoR measures response
**durability** (how long each response holds). They are different, and the difference is
exactly where the ORR surrogate breaks: a drug can produce a response in *everyone* and have
those responses *all* collapse in months. Onkos computes DoR over the same ensemble as ORR
and shows the model with the **highest ORR has the shortest DoR** — and that durability
deficit is precisely why it has the worst tail-driven survival.

> The clinical lesson this encodes is the immunotherapy lesson: a checkpoint inhibitor with a
> modest response rate but durable responses can beat a targeted agent with a high response
> rate but short ones, on overall survival. ORR alone ranks them the wrong way; DoR is the
> dimension that explains it. Onkos makes "depth is not durability" a measured, tested
> quantity rather than a slogan.

---

## 1. The problem this extends

A RECIST response has two independent properties, and the dataset only surfaced one:

| Property | What it asks | Endpoint | Status before this spec |
| --- | --- | --- | --- |
| **Breadth / depth** | how many patients respond (and how deeply)? | ORR, depth of response | ✅ (v0.27, v0.8) |
| **Durability** | how *long* does a response last before progression? | **DoR** | ⚠️ a per-trajectory metric only (`duration_of_response_weeks`), never a population endpoint |
| **Breadth × durability → OS** | does a response rate that ignores durability predict survival? | the surrogate mechanism | ⚠️ shown to fail (v0.27) but the *cause* was not isolated |

ORR's conditional surrogacy (v0.27) is not mysterious once durability is on the table: ORR
counts responders at their best timepoint and forgets them. A model whose responses are deep
but brief (a resistant subclone regrowing) scores a perfect ORR and a short survival under
any tail-sensitive endpoint. DoR is the quantity that was missing — the population durability
of response — and it converts the v0.27 *observation* (ORR mis-ranks tail-driven OS) into a
*mechanism* (the high-ORR model's responses are not durable).

**Why this is the right deepening (and the right scope).** It (1) completes the response
endpoint the project just opened and isolates the cause of the v0.27 surrogate failure; (2)
is pure post-processing — DoR is read from the *same* RECIST episode as the response category,
over the *same* ensemble, with no new kernel or dataset record (mirrors `response`/`budget`);
(3) has direct precedent — DoR is a standard RECIST-derived secondary endpoint (§7); (4) is
*safe by construction* — a trial-level median with explicit censoring, no individual
prediction, no therapy ranking; and (5) sharpens the honest message: *the highest response
rate can be the least durable*, the computational core of "high ORR, no OS benefit."

---

## 2. Duration of response, consistently defined

DoR is computed from the **same observed-baseline RECIST episode** as the best-response
category, so the two are mutually consistent (`response_episode` returns both). For a
responder (CR/PR):

```
PR onset  = first time the SLD falls to ≤ 70% of the observed baseline v[0]
progression = first post-nadir time the SLD rises to ≥ 120% of the nadir
DoR       = progression − onset            (weeks)
```

DoR is `nan` for a non-responder (never reached PR) and **right-censored** for a responder
whose response never progresses within the simulation horizon (a durable responder — exactly
the patients DoR most wants to count). Over the IIV ensemble:

```
n_responders          = # samples achieving CR or PR
median_dor_weeks      = median of the *observed* (uncensored) responder DoRs
dor_censored_fraction = responders without progression / n_responders
```

The observed median is a **lower bound** when censoring is high (the censored responders are
the most durable); a `dor_heavily_censored` warning fires above 50% so the bound is never
mistaken for the truth. This is the honest, simple summary; a Kaplan-Meier median is a clean
later refinement.

---

## 3. The result — depth is not durability

For NSCLC first line, ORR (breadth) and median DoR (durability) point in *different*
directions (illustrative, weeks):

| Model | ORR (breadth) | median DoR (durability) | censored | k_g-OS (tail-driven) |
| --- | --- | --- | --- | --- |
| Two-population (mechanistic resistance) | **1.00** (highest) | **32** (shortest) | 0% | **34** (shortest) |
| Claret (phenomenological resistance) | 0.96 | 62 | 4% | 46 |
| Norton-Simon (complete responder) | 0.68 | 56 (+ 24% censored, i.e. durable) | 24% | **102** (longest) |
| Wang biexponential | 0.45 | 22 | 0% | 43 |

The headline is the first row: the model with the **highest ORR has the shortest DoR**. Its
responses are universal but brief — a deep early shrink followed by a fast resistant
regrowth — so under a survival endpoint that sees the tail (the v0.25 k_g link) it has the
*worst* OS despite the *best* response rate. DoR is the dimension that reconciles them: it
converts the v0.27 surrogate inversion from a paradox into a mechanism. Conversely the
complete responder (Norton-Simon) has a lower ORR but its responses are durable (high
censoring = no progression), matching its best tail-driven OS — the immunotherapy pattern.

**The honest framing.** Onkos does not claim DoR is a *better* surrogate (it sees only
responders and carries censoring). It claims something stronger and safer: ORR and DoR are
*orthogonal* properties of a response, and the ORR-surrogate failures are *durability*
failures — visible the moment durability is measured.

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** DoR carries the propagated tier of the simulated chain;
  out-of-context transport floors to D, exactly as for ORR and OS.
- **Population / trial level only.** DoR is a trial-level median over responders, never an
  individual's response duration.
- **No therapy ranking.** The depth-vs-durability comparison ranks *models under a context*,
  never treatments; "durable" is a property of a model's responses, not a recommendation.
- **Censoring is surfaced, not hidden.** The observed median is reported with its censored
  fraction and a lower-bound warning; the durable responders are never silently dropped to
  make a model look worse (or a finite median look complete).

---

## 5. Validation landmarks

No new kernel — DoR is RECIST arithmetic over the existing trajectories and ensemble. The
landmark suite (`tests/test_duration.py`) pins the behavior:

| Landmark | Condition |
| --- | --- |
| **Episode consistency** | `best_response(t, v) == response_episode(t, v)[0]` for every trajectory. |
| **Non-responder DoR is nan** | an SD/PD trajectory (no PR) has `dor = nan`. |
| **Closed-form DoR** | a constructed PR-at-onset → PD trajectory yields `dor = t_PD − t_PR` to the grid. |
| **Censoring** | a responder that never regrows ⇒ `dor = nan` (censored), and a high censored fraction raises the warning. |
| **Durability orders by regrowth** | a slower-regrowing model has a longer median DoR than a fast-regrowing one. |
| **Depth ≠ durability** | the highest-ORR NSCLC model has a *shorter* median DoR than a lower-ORR model — breadth and durability dissociate. |
| **DoR explains the surrogate** | the model that is ORR→OS discordant under k_g (v0.27) is the short-DoR model. |
| **Rates carry DoR honestly** | `0 ≤ dor_censored_fraction ≤ 1`; `n_responders` matches the CR+PR count; tier passthrough. |

---

## 6. API, CLI, and surface

DoR rides on the existing response surfaces (no new module):

```python
rr = onkos.objective_response_rate(ds, "resistance.nsclc_first_line.two_population", context=ctx)
rr.orr                        # 1.00 — breadth (how many respond)
rr.median_dor_weeks           # 32   — durability (how long responses last)
rr.dor_censored_fraction      # responders without progression (right-censored)
rr.n_responders

rs = onkos.response_vs_survival(ds, context=ctx)   # rows now carry median_dor_weeks
[(r["record_id"].split(".")[1], r["orr"], r["median_dor_weeks"]) for r in rs.rows]
```

**CLI.** `onkos response <id>` now reports ORR, DCR, **and median DoR (with censoring)**, and
`onkos response --durability` prints the per-model breadth-vs-durability table.

**No new export model** — DoR is an analysis over a trajectory, not a model.

---

## 7. Source anchors (methodological; DOIs added at curation time)

- **RECIST 1.1 & DoR.** Eisenhauer, E.A. et al. (2009), *New response evaluation criteria in
  solid tumours: revised RECIST guideline (version 1.1)*, EJC — the response and progression
  definitions DoR is built from (already cited for `response`).
- **Depth vs durability and OS.** The literature on duration of response as an OS-relevant
  endpoint distinct from response rate (e.g. durable-response endpoints in immuno-oncology) —
  the basis for treating DoR as the durability dimension ORR omits.
- **Censored time-to-event summaries.** Standard survival-analysis treatment of right
  censoring — why the observed median DoR is reported as a lower bound under heavy censoring.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not a Kaplan-Meier DoR (yet).** The observed median + censored fraction is the v0.x
  summary; a KM median that uses the censored durable responders is a clean refinement.
- **Not a claim that DoR is a valid OS surrogate.** DoR sees only responders and carries
  censoring; Onkos uses it to *explain* the ORR-surrogate failure, not to replace it.
- **Not individual response duration.** Trial-level median only.
- **No therapy ranking.** The breadth-vs-durability comparison ranks models, never treatments.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only.** ORR, DoR, and the comparison are trial-level quantities
  over published models; nothing is an individual's response or its duration.
- **No therapy ranking, no recommendation.** "Durable" describes a model's responses, not a
  treatment choice.
- **Censoring and lower bounds are explicit.** The observed median is never presented as the
  complete truth when responders remain unprogressed.
- **Cannot raise a tier.** DoR carries the worst tier of the simulated chain.
- **The line, restated.** Any feature that takes a real patient's tumor measurement and
  returns a response duration or a therapy choice **does not get built.** Making population
  durability and its link to the ORR-surrogate failure computable changes none of this.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Consistent episode** | `response_episode(t, v)` returning best response *and* DoR from one observed-baseline trajectory; `best_response` delegates to it. | Category and DoR are mutually consistent; the closed-form DoR landmark passes. |
| **2 — Population DoR** | `median_dor_weeks`, `dor_censored_fraction`, `n_responders` on `ResponseRates`, computed over the ensemble; the heavy-censoring warning. | ORR and DoR are reported together; censoring is surfaced. |
| **3 — Depth ≠ durability** | DoR on the `response_vs_survival` rows; the depth-vs-durability finding (highest ORR = shortest DoR = the k_g-discordant model). | The breadth/durability dissociation and its tie to the surrogate are shown. |
| **4 — Surfaces** | `onkos response` reports DoR; `--durability` table; a breadth-vs-durability figure + a CI-executed notebook; README section. | DoR is reachable, visualized, and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested population duration
of response with honest censoring. Step 3 is the payload: showing, on the same simulated
trial, that the highest response *rate* can be the least durable — and that this durability
deficit is the mechanism behind the contested ORR → OS surrogate — is the quantitative core
of "depth is not durability," shipped as a tested artifact.
