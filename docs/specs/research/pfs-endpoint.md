# Onkos — research spec: the PFS endpoint (two routes to progression-free survival)

**Status:** implemented in v0.30.0 (`onkos.response` — mechanistic PFS over the ensemble,
beside the statistical PFS link). This is the design-of-record; the methodological source
anchors of §7 are documented but their Crossref-verified citation curation is still pending,
honest by design. Written in the v0.1 house style; every value is illustrative and
`unverified` by design — the infrastructure is the contribution.

**PFS is the endpoint that gates accelerated approvals, and it can be computed two ways that
disagree.** v0.27 made ORR a *measured* OS surrogate; v0.28 supplied DoR as the durability
dimension ORR omits. Both live inside the **response** endpoint. This spec opens the
*progression-free survival* endpoint as the response module's sibling — and in doing so
exposes a model-selection axis that sits *inside a single endpoint*. Onkos already produces a
PFS number from the **parametric PFS survival link** (a week-8-keyed hazard model). It can
*also* produce a PFS number **mechanistically**, by reading the RECIST progression time
directly off the simulated tumor trajectory. These two routes are calibrated from different
data and need not agree — and for shrink-then-regrow resistance dynamics they invert, in
*every* solid-tumor context the dataset covers.

> The clinical stakes are specific to PFS. PFS gates accelerated and conditional approvals;
> it is the endpoint a sponsor reaches for when OS is immature. A PFS that looks long because
> a week-8-keyed hazard model never saw the resistant regrowth — while the tumor model itself
> says progression comes early — is the model-selection risk that, unquantified, turns a
> "positive" PFS trial into a negative confirmatory OS trial. Onkos makes the two routes
> to PFS, and their disagreement, a measured, tested quantity.

---

## 1. The problem this extends

Onkos already has a PFS endpoint — but only one route to it. The progression event has two
legitimate computational definitions, and the dataset surfaced only the statistical one:

| Route to PFS | What computes it | Sees the regrowth tail? | Status before this spec |
| --- | --- | --- | --- |
| **Statistical** | the parametric PFS survival link — a hazard model keyed on the week-8 tumor change | **No** — it samples week 8 and extrapolates a hazard | ✅ (v0.12, every context) |
| **Mechanistic** | RECIST time-to-progression read off the simulated tumor-size trajectory | **Yes** — it watches the SLD cross +20% of its nadir | ⚠️ not computed (the progression logic existed only inside DoR) |

The statistical route is the standard one (PFS hazard models are fit to week-8 or early
landmark change because that is what an early read affords). But it inherits the same blind
spot the v0.25/v0.27/v0.28 arc kept finding: **a week-8-keyed link cannot see a tail it never
sampled.** A model whose tumor shrinks deeply by week 8 and then regrows (a resistant subclone
— the two-population mechanism of v0.24) is scored *long* PFS by the statistical link and
*short* PFS by the mechanism that actually watches the regrowth. The progression arithmetic to
do this already lived inside `response_episode` (DoR needs a progression time); this spec
promotes it to a first-class endpoint and reads off the disagreement.

**Why this is the right deepening (and the right scope).** It (1) completes the endpoint
family — OS (v0.12), ORR (v0.27), DoR (v0.28), and now **PFS computed mechanistically** beside
the statistical link — so every endpoint a phase-2/3 trial reports is in the dataset; (2) is
pure post-processing — the mechanistic route is RECIST arithmetic over the *same* trajectory
and ensemble, the statistical route is the *existing* PFS link, no new kernel or record
(mirrors `response`/`budget`/`identify`); (3) lands the model-selection thesis *inside a single
endpoint* — v0.25 showed the OS metric is a choice, v0.26 budgeted it; this shows the **PFS
route** is a choice, and which route you take inverts the model ranking; (4) is *dataset-wide
from the first commit* — every solid-tumor context already has both a PFS link (v0.12) and a
two-population model (v0.29), unlike the NSCLC-only ORR surrogate; (5) is *safe by
construction* — a trial-level median and a landmark rate with explicit censoring, no individual
prediction, no therapy ranking.

---

## 2. Mechanistic PFS, consistently defined

Time-to-progression is read from the **same tumor-size trajectory** the rest of the response
endpoint uses, with the standard RECIST 1.1 progression rule applied against the **running
nadir** (the smallest SLD recorded up to each time, baseline included):

```
running_nadir(t) = min over s ≤ t of SLD(s)
progression      = first t > 0 with SLD(t) ≥ 1.20 × running_nadir(t)
TTP (mech. PFS)  = that progression time, in weeks      (nan if none within the horizon)
```

A trajectory that never crosses +20% of its running nadir within the simulation horizon is
**right-censored** (`nan`) — a durable non-progressor, exactly the patient PFS most wants to
count. Over the IIV ensemble two honest summaries are reported:

```
median_ttp_weeks        = median of the *observed* (uncensored) TTPs        (a lower bound under censoring)
ttp_censored_fraction   = samples without progression / n
mechanistic_pfs_rate    = P(progression-free at the landmark)  = mean(TTP is nan OR TTP > landmark_weeks)
```

The **landmark rate** is the censoring-robust summary: it is a fixed-horizon
progression-free probability (default 24 weeks ≈ 6 months, the canonical PFS landmark) that
counts the censored durable non-progressors correctly, whereas the observed median is a lower
bound when censoring is high (a `ttp_heavily_censored` warning fires above 50%, mirroring DoR).

The statistical route is unchanged: it is the context's default PFS survival link, summarized
as the ensemble median of `s.median["PFS"]`.

---

## 3. The result — the route is the model-selection axis

For NSCLC first line, the two PFS routes disagree, and they **invert the model ranking**
(illustrative, weeks; n = 300, horizon 156 wk — the module defaults):

| Model | mechanistic median TTP | statistical median PFS (week-8 link) | route ratio |
| --- | --- | --- | --- |
| Claret (phenomenological resistance) | **60** | 34 | 1.74 |
| Norton-Simon (complete responder) | 36 | 22 | 1.59 |
| Two-population (mechanistic resistance) | **34** | **35** | **0.95** |
| Wang biexponential | 28 | 21 | 1.33 |

Read the Claret vs two-population pair. **Mechanistically**, Claret progresses far *later* than
the two-population model (60 vs 34 wk) — Claret's phenomenological resistance regrows gently,
the two-population's resistant clone regrows fast. **Statistically**, the order flips: the
two-population model gets the *longer* PFS (35 vs 34 wk), because at week 8 it is deeply shrunk
and the week-8-keyed hazard link reads that deep early shrinkage as durable benefit — blind to
the regrowth it never sampled. Its route ratio (0.95) is the only one near 1.0 not because the
routes *agree* but because a short mechanistic PFS and a long statistical PFS happen to cross;
every other model has a mechanistic PFS well *above* its statistical PFS (ratio > 1.3). **The
PFS endpoint's verdict on which model looks better depends entirely on which route computes
it** — the two-population model is dead last mechanically (tied shortest) but top statistically.
This is the v0.25 OS-metric inversion, now living *inside* the PFS endpoint.

The finding is **dataset-wide**. The two-population model is the consistent culprit — shortest
or near-shortest mechanistic PFS, yet among the longest statistical PFS — and the route
inverts its rank against the Claret model in *every* context:

| Context | mech. TTP: Claret / two-pop | stat. PFS: Claret / two-pop | route-discordant pairs |
| --- | --- | --- | --- |
| NSCLC, 1L | 60 / 34 | 34 / 35 | 2/6 |
| Breast, 1L | 78 / 39 | 65 / 67 | 1/3 |
| CRC, 1L | 67 / 32 | 53 / 55 | 1/3 |
| HCC, 1L | 48 / 27 | 30 / 31 | 1/3 |
| Melanoma, 1L | 62 / 28 | 45 / 46 | 1/3 |

In every row the mechanistic route ranks Claret well above the two-population model, and the
statistical route ranks them level or reversed. The cause is the same one the whole arc keeps
surfacing: **the week-8 link is blind to the tail; the mechanism is not.** Because the
two-population model (v0.29) and the PFS link (v0.12) are present in all five contexts, the
finding reproduces dataset-wide on the day it ships — it is not an NSCLC artifact.

**The honest framing.** Onkos does not claim the mechanistic route is the *correct* PFS (the
statistical link encodes real trial-fit hazard information the mechanism lacks; the mechanism
assumes the tumor model is true). It claims something stronger and safer: **PFS is not a single
number** — it carries a route choice, that choice is a model-selection axis exactly like the
OS-metric choice, and for resistance dynamics the two routes invert. The disagreement is the
quantity; neither route is privileged.

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins governs.** Both PFS routes carry the propagated tier of the simulated
  chain; out-of-context transport floors to D, exactly as for ORR, DoR, and OS.
- **Population / trial level only.** Mechanistic TTP is a trial-level median and a landmark
  rate over the ensemble, never an individual's time to progression.
- **No therapy ranking.** The route comparison ranks *models under a context*, never
  treatments; a "longer PFS" is a property of a model under a route, not a recommendation.
- **Censoring is surfaced, not hidden.** The observed median is reported with its censored
  fraction and a lower-bound warning; the landmark rate counts the durable non-progressors.
- **Cannot raise a tier.** Neither route can lift the worst tier of the chain.

---

## 5. Validation landmarks

No new kernel — mechanistic PFS is RECIST arithmetic over the existing trajectories and
ensemble. The landmark suite (`tests/test_pfs_routes.py`) pins the behavior:

| Landmark | Condition |
| --- | --- |
| **Closed-form TTP** | a constructed shrink-then-regrow trajectory crosses +20% of its nadir at a known time ⇒ `time_to_progression` returns it to the grid. |
| **Running-nadir rule** | a monotone-growth trajectory progresses against baseline; a deep-shrink-then-regrow trajectory progresses against the *nadir*, not baseline (the running-nadir distinction). |
| **Censoring** | a trajectory that never regrows past +20% of its nadir ⇒ `nan` (censored), and a high censored fraction raises the warning. |
| **Landmark rate bounds** | `0 ≤ mechanistic_pfs_rate ≤ 1`; an all-progressing-early ensemble ⇒ low rate, an all-durable ensemble ⇒ rate 1. |
| **Both routes present** | for a context with a PFS link, `progression_free_survival` returns both a finite `median_ttp_weeks` and a finite `median_pfs_link_weeks` (`has_pfs_link`). |
| **Route inversion** | for NSCLC 1L the two-population model has a *shorter* mechanistic TTP than Claret but a *longer* statistical PFS — the routes disagree. |
| **Dataset-wide discordance** | `pfs_route_divergence` reports ≥1 route-discordant model pair in every solid-tumor context. |
| **Tier passthrough & guardrails** | the propagated tier rides through; the clinical-use flag is on every payload. |

---

## 6. API, CLI, and surface

Mechanistic PFS rides on the existing response module (no new module):

```python
pf = onkos.progression_free_survival(ds, "resistance.nsclc_first_line.two_population", context=ctx)
pf.median_ttp_weeks          # 34 — mechanistic: RECIST progression off the trajectory
pf.median_pfs_link_weeks     # 35 — statistical: the week-8-keyed PFS hazard link, same trial
pf.mechanistic_pfs_rate      # P(progression-free at 24 wk), censoring-robust
pf.ttp_censored_fraction
pf.route_ratio               # mechanistic / statistical median (1.0 = the routes agree)

div = onkos.pfs_route_divergence(ds, context=ctx)   # across the in-context TGI models
div.discordant_pairs, div.total_pairs               # route disagreements in the model ranking
```

**CLI.** `onkos pfs <id>` reports the mechanistic median TTP, the statistical PFS, the landmark
rate, and their ratio; `onkos pfs --routes` prints the per-model two-route table and the
route-discordance count across the in-context models.

**No new export model** — mechanistic PFS is an analysis over a trajectory, not a model; the
PFS *link* is already exported.

---

## 7. Source anchors (methodological; DOIs added at curation time)

- **RECIST 1.1 progression.** Eisenhauer, E.A. et al. (2009), *New response evaluation criteria
  in solid tumours: revised RECIST guideline (version 1.1)*, EJC — the +20%-from-nadir
  progression rule mechanistic TTP applies (already cited for `response`).
- **PFS as a regulatory endpoint and its OS-surrogate caveats.** The literature on PFS as an
  accelerated-approval endpoint and the conditions under which it does and does not predict OS
  — the basis for treating the PFS *route* as a model-selection risk, not a settled number.
- **Landmark / fixed-horizon event rates.** Standard survival-analysis treatment of landmark
  analysis and right censoring — why the landmark progression-free rate is reported beside the
  observed median.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not a Kaplan-Meier PFS (yet).** The observed median + landmark rate + censored fraction is
  the v0.x summary; a KM median that uses the censored durable non-progressors is a clean
  refinement (the same refinement noted for DoR).
- **Not a claim that either route is the true PFS.** The statistical link carries trial-fit
  hazard information; the mechanism assumes the tumor model. Onkos reports the *disagreement*,
  privileging neither.
- **Not a new survival-link kernel.** This adds no record and no kernel; it reads the existing
  PFS link and the existing trajectory.
- **Not individual time-to-progression.** Trial-level median and landmark rate only.
- **No therapy ranking.** The route comparison ranks models, never treatments.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only.** Both PFS routes and the comparison are trial-level
  quantities over published models; nothing is an individual's progression-free time.
- **No therapy ranking, no recommendation.** "Longer PFS" describes a model under a route, not
  a treatment choice.
- **Censoring and lower bounds are explicit.** The observed median is never presented as the
  complete truth when non-progressors remain.
- **Cannot raise a tier.** Mechanistic PFS carries the worst tier of the simulated chain.
- **The line, restated.** Any feature that takes a real patient's tumor measurement and returns
  a progression-free time or a therapy choice **does not get built.** Making the two routes to
  population PFS, and their disagreement, computable changes none of this.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Mechanistic TTP** | `time_to_progression(t, v)` applying the RECIST +20%-from-running-nadir rule to one trajectory; censoring as `nan`. | The closed-form and running-nadir landmarks pass. |
| **2 — Population PFS** | `progression_free_survival(...)` returning the mechanistic median TTP, landmark rate, and censored fraction over the ensemble, beside the statistical PFS-link median; the heavy-censoring warning. | Both routes are reported together; censoring is surfaced. |
| **3 — The route is an axis** | `pfs_route_divergence(...)` across the in-context TGI models; the route-inversion finding (two-population shortest mechanically, longest statistically). | The route disagreement and its dataset-wide reproduction are shown. |
| **4 — Surfaces** | `onkos pfs` reports both routes; `--routes` table; a two-route figure + a CI-executed notebook; README section. | PFS routes are reachable, visualized, and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested mechanistic
population PFS with honest censoring, sitting beside the statistical PFS link. Step 3 is the
payload: showing, on the same simulated trial, that the PFS endpoint carries a *route* choice
that inverts the model ranking — that PFS is not one number but two that disagree for exactly
the resistance dynamics that matter — is the quantitative core of this spec, shipped as a
tested artifact and reproduced across all five solid-tumor contexts.
