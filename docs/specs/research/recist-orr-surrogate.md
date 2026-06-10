# Onkos — research spec: RECIST response & the ORR → OS surrogate

**Status:** implemented in v0.27.0 (`onkos.response`). This is the design-of-record; the
methodological source anchors of §7 are documented but their Crossref-verified citation
curation is still pending, honest by design. Written in the v0.1 house style; every value
is illustrative and `unverified` by design — the infrastructure is the contribution.

**The dominant phase-2 endpoint, and its contested surrogacy.** The objective response
rate (ORR) — the fraction of a trial achieving a RECIST objective response — is the
workhorse go/no-go readout of early oncology. It is also a *contested* surrogate for
overall survival: drugs with a high response rate routinely fail to extend survival, and
accelerated approvals on ORR are sometimes not confirmed. Onkos has OS and PFS endpoints
but no response endpoint. This adds it — and, because a model's ORR *and* its OS come from
the **same simulated trial**, it makes the ORR → OS relationship a *measured* quantity. The
sharp result: whether ORR faithfully ranks survival is **conditional on the survival
mechanism** — under an early-shrinkage survival model ORR is a perfect surrogate; under a
regrowth-tail survival model it is actively misleading. That is the response-endpoint
analog of the v0.25 survival-metric finding, and the computational explanation for why
high-ORR drugs fail confirmatory OS trials.

> ORR and the week-8 survival surrogate are *both* shrinkage-based, so they agree by
> construction — which is exactly why ORR "works" as a surrogate when survival is assumed
> to track early shrinkage. The moment survival is driven by the regrowth tail (the v0.25
> k_g link), a deep early responder that regrows fast has a *high* ORR and a *short* OS.
> Onkos makes that conditional surrogacy legible instead of leaving it to post-hoc
> disappointment in phase 3.

---

## 1. The problem this extends

Onkos composes `TGI model → metric → survival link → OS/PFS`. The response endpoint sits
*before* the survival link — it reads the tumor-size trajectory directly — and has been
missing:

| Endpoint | What it is | Status before this spec |
| --- | --- | --- |
| **OS / PFS** | population survival from an on-treatment metric | ✅ (v0.12–v0.13, v0.25) |
| **ORR / DCR / best response** | population RECIST response rate from the trajectory | ⚠️ **absent** — yet it is the dominant phase-2 endpoint |
| **ORR → OS surrogacy** | does the response rate faithfully rank survival? | ⚠️ a famous open question, never made computable here |

The gap is not just a missing number; it is a missing *honesty surface*. ORR is the metric
on which the most early go/no-go decisions are made, and its surrogate validity is the
single most consequential assumption behind those decisions. Onkos already computes the
RECIST thresholds (depth of response, the 30%/20% bounds in `onkos.metrics`); the deepening
is to lift them from a per-trajectory metric to a **population rate** and to put that rate
beside the OS read off the same trial.

**Why this is the right deepening (and the right scope).** It (1) fills the dominant
phase-2 endpoint, a genuine gap, and advances the project's surrogate-honesty thesis; (2)
is pure post-processing over the existing IIV ensemble (`ensemble_samples`) and RECIST
constants — no new kernel, no dataset subsystem, near-zero schema change (mirrors
`combine`/`identify`/`budget`); (3) has direct precedent — RECIST 1.1 and the ORR-OS
surrogate literature (§7); (4) is *safe by construction* — population/trial-level rates, no
individual response probability, no therapy ranking; and (5) sharpens the honest message:
its headline is *the response rate can rank treatments the opposite way survival does*,
which is the computational core of every ORR-surrogate failure.

---

## 2. RECIST best response and the population rates

Best overall response is classified per RECIST 1.1 on the tumor-size (sum-of-longest-
diameters) trajectory, from the patient's **observed baseline** `v[0]` (as RECIST measures
it), with `CR > PR > PD > SD` precedence:

```
depth = (v[0] − nadir) / v[0]                          fractional shrinkage at nadir
CR  if depth ≥ 0.95   (near-complete disappearance — an SLD-continuous proxy)
PR  if depth ≥ 0.30   (objective partial response)
PD  if no PR and max-post-nadir ≥ 1.20 × nadir         (≥20% regrowth = progression)
SD  otherwise                                          (stable disease)
```

Over the stored inter-individual variability (the same lognormal ensemble
`onkos.uncertainty` runs), each sampled trajectory yields a best response, and the
population rates are the category fractions:

```
ORR = P(CR or PR)                 the objective response rate
DCR = P(CR, PR, or SD) = 1 − P(PD)   the disease-control rate
distribution = (CR, PR, SD, PD)   sums to 1
```

These are **trial-level rates**, not an individual's response probability. A model with no
reported IIV yields a degenerate distribution (every sample identical → ORR ∈ {0, 1}); the
spread across categories *is* the inter-patient response heterogeneity the dataset encodes.

---

## 3. The ORR → OS surrogate, made computable

Because the OS curve is read off the *same* ensemble, every in-context TGI model has both an
ORR and a median OS. `response_vs_survival` counts the **discordant model pairs** — pairs
where one model has a *higher* ORR yet *shorter* OS than another:

```
discordant_fraction = (# pairs with ORR order ≠ OS order) / (# comparable pairs)   ∈ [0, 1]
```

A nonzero discordant fraction is direct evidence that ORR does not faithfully rank survival
in this context — the exact failure mode (high response, no survival benefit) that sinks
confirmatory trials. The result for NSCLC first line is the spec's sharpest:

| Survival model (the OS link) | ORR ranks OS? | discordant pairs | reading |
| --- | --- | --- | --- |
| **week-8 change** (early-shrinkage surrogate) | ✅ perfectly | 0 / 6 | ORR and OS are both shrinkage-based — ORR "works". |
| **k_g** (regrowth-tail metric, v0.25) | ❌ badly | 4 / 6 (67%) | the highest responder (deep early shrinkage, ORR ≈ 1.0) has the *shortest* OS; the eradicating drug (lower ORR, no regrowth) the *longest*. |

So **ORR's surrogate validity is conditional on the survival mechanism**, and the condition
is unobservable from the early trial. This is not a contradiction in Onkos; it is Onkos
making explicit the assumption that every ORR-based go/no-go silently makes — that survival
tracks early shrinkage — and showing the decision reverses when it does not.

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** The rates carry the propagated tier of the simulated chain;
  an out-of-context transport floors to D, exactly as for OS.
- **Population / trial level only.** ORR, DCR, and the distribution are trial-level rates of
  a published model; no individual response probability is emitted.
- **No therapy ranking.** The discordance is a statement about *models under one context*,
  never a recommendation between treatments, and never "this drug responds better."
- **The surrogate is shown, not asserted.** Onkos does not declare ORR valid or invalid; it
  reports the discordance under each survival assumption and lets the conditional nature
  speak. The headline output is a *disagreement*, the opposite of false precision.

---

## 5. Validation landmarks

No new kernel — the classifier is RECIST arithmetic over existing trajectories, and the
rates are population fractions over the existing ensemble. The landmark suite
(`tests/test_response.py`) pins the behavior:

| Landmark | Condition |
| --- | --- |
| **PR boundary** | a trajectory with exactly 30% shrinkage and no regrowth ⇒ `PR`; 29% ⇒ `SD`. |
| **CR threshold** | ≥95% shrinkage ⇒ `CR`. |
| **PD requires no PR + regrowth** | <30% shrink then ≥20% regrowth from nadir ⇒ `PD`; a PR that later regrows is still `PR` (best overall response). |
| **Monotone stable** | a flat trajectory ⇒ `SD`. |
| **Rates are a simplex** | the CR/PR/SD/PD distribution sums to 1; `0 ≤ ORR ≤ DCR ≤ 1`. |
| **ORR monotone in effect** | a larger drug effect never lowers ORR. |
| **Degenerate ensemble** | zero IIV ⇒ every sample identical ⇒ ORR ∈ {0, 1}. |
| **Surrogate concordance under week-8** | ORR perfectly ranks OS under the shrinkage-based link (0 discordant pairs). |
| **Surrogate discordance under k_g** | ORR mis-ranks OS under the tail-sensitive link (>0 discordant pairs) — the conditional-surrogacy result. |
| **Tier passthrough** | the rates carry the chain's propagated tier (out-of-context ⇒ D). |

---

## 6. API, CLI, and surface

```python
rr = onkos.objective_response_rate(ds, "resistance.claret_2009.tgi",
                                   context=dict(tumor_type="NSCLC", line="first"), n=300)
rr.orr, rr.dcr               # objective-response and disease-control rates
rr.distribution              # {"CR":…, "PR":…, "SD":…, "PD":…}  (sums to 1)
rr.median_os_weeks           # OS from the SAME trial — for the surrogate question
rr.tier, rr.warnings

rs = onkos.response_vs_survival(ds, context=ctx)                       # week-8 OS link
rs.discordant_fraction, rs.orr_predicts_os                            # the surrogate verdict
rs_kg = onkos.response_vs_survival(ds, context=ctx,
                                   survival_link="survival_link.nsclc_os_growth_rate")
rs_kg.discordant_fraction    # >0 — ORR mis-ranks tail-driven survival
```

**CLI.**

```bash
onkos response resistance.claret_2009.tgi              # ORR / DCR / RECIST distribution + OS
onkos response --surrogate --tumor-type NSCLC          # the ORR -> OS discordance table
onkos response --survival-link survival_link.nsclc_os_growth_rate --surrogate
```

**No new export model** — response is an analysis over a trajectory, not a model. (The
RECIST distribution could ride along in the virtual-trial JSON; deferred until a consumer
needs it.)

---

## 7. Source anchors (methodological; DOIs added at curation time)

- **RECIST 1.1.** Eisenhauer, E.A. et al. (2009), *New response evaluation criteria in solid
  tumours: revised RECIST guideline (version 1.1)*, European Journal of Cancer — the
  response-category definitions (the 30% / 20% thresholds) the classifier implements.
- **ORR as a contested OS surrogate.** The trial-level surrogacy literature (e.g.
  meta-analyses correlating ORR with OS across randomized trials) — the basis for treating
  ORR → OS as a measured, context-dependent relationship rather than an assumed one.
- **Tumor-size endpoints and survival.** Claret/Bruno model-based OS prediction (already
  cited) — the framework into which the response endpoint slots.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not individual response prediction.** ORR is a trial-level rate; no per-patient response
  probability is produced.
- **Not full RECIST lesion bookkeeping.** A single SLD trajectory, not target/non-target
  lesion rules, new-lesion appearance, or confirmation scans; CR is an SLD-continuous proxy
  for disappearance. The categories are faithful to the size dynamics, not to the full
  radiology workflow.
- **Not a surrogate *validation*.** Onkos reports discordance under stated survival
  assumptions; it does not certify ORR as a valid or invalid surrogate (that needs
  randomized-trial meta-data, not a simulator).
- **No therapy ranking.** The discordance ranks *models*, never treatments.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only.** Every quantity is a trial-level rate or a model-level
  comparison; nothing is an individual's response or prognosis.
- **No therapy ranking, no recommendation.** ORR and the discordance are statements about
  published models under a context, never a choice between treatments.
- **The surrogate is shown, not asserted.** The headline is a *disagreement* that depends on
  an unobservable survival assumption — the opposite of false precision.
- **Cannot raise a tier.** The rates carry the worst tier of the simulated chain.
- **The line, restated.** Any feature that takes a real patient's tumor measurement and
  returns a response or a therapy choice **does not get built.** Making the population
  response rate and its conditional surrogacy computable changes none of this.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — RECIST classifier + rates** | `best_response` (RECIST 1.1 best overall response) and `objective_response_rate` (ORR / DCR / distribution over the IIV ensemble); the classification + simplex landmarks. | The rates compute and every classification landmark passes. |
| **2 — The surrogate** | `response_vs_survival` with the discordant-pair count; the OS read off the same trial; the conditional-surrogacy result under the week-8 vs k_g links. | The discordance is 0 under week-8 and >0 under k_g. |
| **3 — Surfaces** | `onkos response [--surrogate --survival-link]`; `to_dict` with the clinical-use flag; a RECIST-distribution + ORR-vs-OS figure and a CI-executed notebook. | The response endpoint and the surrogate are reachable and visualized. |
| **4 — Documentation** | README section framing ORR and its conditional surrogacy; roadmap, cheat sheet, and layout updated. | The response endpoint is documented with the project's rigor and guardrails. |

Step 1 alone is a self-contained contribution: an open, landmark-tested population RECIST
response rate from a curated TGI model. Step 2 is the payload: showing, on the same
simulated trial, that the field's dominant phase-2 endpoint can rank treatments the
opposite way survival does — and that whether it does is conditional on the survival
mechanism nobody observes early — is the computational core of the ORR-surrogate problem,
shipped as a tested artifact.
