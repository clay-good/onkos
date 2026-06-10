# Changelog

All notable changes to Onkos are documented here. Versions follow the phased
roadmap (spec §11). All parameter values are illustrative and `unverified` by
design; the infrastructure is real and tested.

## [0.36.0] — Exposure-response model choice: the dose-extrapolation model-selection axis

Implements the research-track spec `docs/specs/research/exposure-response-extrapolation.md`. Every chain
so far started from a *given* drug effect; that effect is the output of an exposure-response (ER) model
(Emax / sigmoid-Emax / power) that all fit the studied dose comparably but **diverge** when extrapolated
to an untested dose — which is exactly what dose selection asks of them. This makes the ER-model choice
the project's transportability thesis one layer upstream, with the *dose* as the context.

- **New module `onkos.dose_response`** (pure post-processing — no new record, kernel, schema, or export;
  every default artifact byte-identical). `calibrated_er` re-anchors each curated ER shape so it passes
  through a reference `(c_ref, e_ref)` (the studied dose): the shape parameters (EC50, gamma, theta) are
  kept from the curated record, the single scale parameter (Emax or slope) is solved to hit the anchor.
  The curves are then identical at the studied dose and differ only in how they extrapolate.
- **`compare_er_extrapolation`** runs the shapes over a dose grid through the existing TGI→survival chain
  and reports the effect and OS spread per dose, with `reference_os_divergence` (~0, the control) and
  `max_os_divergence`.
- **The finding — invisible at the studied dose, a model-selection axis off it.** Anchored at
  `(150, 1.0)` for Claret NSCLC, the OS spread across ER shapes is **0 at the studied dose** and grows on
  extrapolation: **19 wk at quarter-dose**, 14 wk at half-dose, 5 wk at 2–4×. The risk is **asymmetric** —
  sharpest on de-escalation (the effect sits on the steep part of the survival relationship), where
  dose-finding actually lives; upward extrapolation gives a larger *effect* spread (the unbounded power
  curve runs to E≈3.25 at 4×) but a smaller *OS* spread (survival saturates). A dose-response model fit at
  one dose carries an unquantified model-selection risk the moment it picks another.
- **11 landmarks** (`tests/test_dose_response.py`): the exact anchor, calibration monotonicity, off-anchor
  divergence, the zero-divergence control at the studied dose, the de-escalation-diverges-most asymmetry,
  the single-model-zero-divergence degenerate case, and the inherited tier/transport guardrails.
- **CLI `onkos dose-response`**; an ER-curves + OS-spread figure (`docs/images/dose_response_extrapolation.png`);
  a CI-executed notebook (`notebooks/29_dose_response_extrapolation.ipynb`); README section; public-API
  surface + contract test extended. No new dataset records. 384 tests, 56 records.

## [0.35.0] — Dose-level Loewe additivity: the additivity reference as a model-selection axis

Implements the research-track spec `docs/specs/research/loewe-additivity.md`. v0.23 made the
*interaction model* a model-selection axis at the **effect level** (combine two effect magnitudes under
HSA / Bliss-additive / Greco). This adds the classical gold-standard null it deferred — **dose-level
Loewe additivity** — which combines two *doses* through the dose-response curves and is the only
"no-interaction" reference that is self-consistent. So *which reference you call "no interaction"* is
itself a model-selection axis.

- **`onkos.interaction` extension** (pure post-processing — no new record, kernel, schema, or export;
  every default artifact byte-identical). `ERCurve` + `er_curve` expose each curated ER kernel as a
  forward map and its **analytic inverse** (Emax, sigmoid-Emax, power). `loewe_effect` solves the
  dose-additive isobole `d_A/D_A(E) + d_B/D_B(E) = 1` by bisection (no new dependency). `combine_doses`
  combines two doses under `hsa` / `bliss` / `loewe`; `compare_additivity_references` runs all three
  through the existing TGI→survival chain and reports the OS divergence.
- **The sham-combination identity** is the correctness anchor: a drug combined with itself is *exactly*
  Loewe-additive (`Loewe(d_A,d_B) = f(d_A+d_B)`), which **Bliss fails** for any saturating curve
  (`f(d_A)+f(d_B) > f(d_A+d_B)`). This is why Loewe is the principled reference and the choice is not
  cosmetic.
- **The finding — three references, one dose pair, three survival curves.** For Claret NSCLC (`d_A=150`,
  `d_B=90`): combined effect HSA 0.90 / Loewe 1.07 / Bliss 1.60, median OS 88 / 92 / 101 wk. The ordering
  is structural for saturating curves (`HSA ≤ Loewe ≤ Bliss`): Bliss *overstates* (its combined effect
  1.60 even exceeds the shared effect ceiling 1.4 — the classic effect-additivity artifact), HSA
  *understates*, Loewe is the self-consistent middle. The disagreement **grows with dose** (negligible at
  low dose, large in saturation) — a dose-dependent model-selection risk where combination dose-finding
  lives.
- **11 landmarks** (`tests/test_loewe.py`): the sham identity (exact), Bliss failing it, single-agent
  limits, ER-inverse round-trips, effect-ceiling clamp, dose monotonicity, a record-free core, the
  reference ordering, the OS divergence, and the inherited tier/transport guardrails.
- **CLI `onkos loewe`**; a dose-scaling + OS figure (`docs/images/loewe_additivity.png`); a CI-executed
  notebook (`notebooks/28_loewe_additivity.ipynb`); README section; public-API surface + contract test
  extended. No new dataset records. 373 tests, 56 records.

## [0.34.0] — Joint longitudinal–survival modeling: the current-value link

Implements the research-track spec `docs/specs/research/joint-survival.md`. Every survival link Onkos
has shipped is **two-stage**: collapse the tumor trajectory to one static covariate (week-8 v0.12, `k_g`
v0.25, integrated burden v0.33), then apply a parametric/Cox baseline. A static covariate means a
**proportional hazard** — a constant hazard ratio over time. The joint longitudinal–survival model (the
rigorous, gold-standard alternative to a two-stage landmark) relaxes exactly that. This adds its
canonical **current-value** link and shows "two-stage vs joint" is a model-selection axis at the
survival-link layer.

- **New module `onkos.joint`** (pure post-processing — no record, kernel, schema, or export change, so
  every default artifact is byte-identical). The current-value link makes the instantaneous hazard track
  the *current* tumor size: `λ(t) = λ₀(t)·exp(α·log(v(t)/y0))`, `S(t)=exp(-∫λ)`. The pure core
  `current_value_survival` integrates the time-varying hazard ratio as a Stieltjes sum against the
  **analytic** baseline cumulative hazard `H₀(t)=(t/scale)^shape`, so it is a strict generalization,
  exact in two limits: a constant HR recovers the two-stage Weibull-PH curve to machine precision, and
  `HR≡1` recovers the Weibull baseline. The v0.33 burden link is the constant-trajectory special case.
- **`joint_survival` / `compare_joint_vs_two_stage`** bind it to a record + context: the trajectory and
  its two-stage (week-8) OS come from `simulate`; the Weibull baseline (`shape`/`scale`) from the
  context's default Weibull OS link; the association `α` is a declared argument (never fitted). Tier and
  transport warnings ride through unchanged.
- **The finding — a non-proportional hazard, and a ranking inversion.** For NSCLC first line (`α=1`) the
  joint hazard ratio is suppressed during the deep early response (HR ≈ 0.13–0.18) then **rises 10× to
  255×** as the resistant clone regrows (largest for acquired resistance and two-population), while a
  complete responder's hazard keeps *falling* — a time-varying hazard ratio no two-stage (PH) link,
  parametric or Cox, can represent. And it **inverts** the week-8 ranking: the surrogate ranks the
  deep-early-shrinking two-population model above Claret (94 vs 91); the joint link, weighting the
  regrowth tail, reverses it (Claret 199 vs two-population 144). The survival-link *structure* choice is
  the structural counterpart to the v0.25/v0.33 *metric* axis.
- **11 landmarks** (`tests/test_joint.py`): the two exact-recovery limits, `α=0`, monotonicity, the
  burden-link bridge, the non-proportional-hazard signature, the eradication mirror, the ranking
  inversion, and the unchanged tier/transport guardrails.
- **CLI `onkos joint`** (the two-stage-vs-joint median + PH-violation table); a hazard-ratio + survival
  figure (`docs/images/joint_survival.png`); a CI-executed notebook (`notebooks/27_joint_survival.ipynb`);
  README section; public-API surface + contract test extended. No new dataset records. 362 tests, 56
  records.

## [0.33.0] — The integrated tumor burden: a third TGI→OS bridge metric

Implements the research-track spec `docs/specs/research/burden-auc-bridge-metric.md`. v0.25 made the
on-treatment metric that drives a survival link a declared, swappable field (`structure.link_metric`)
with two values: the default **week-8 relative change** (depth-only, blind to the tail) and the
**growth-rate constant `k_g`** (tail-only, blind to depth) — and showed the metric choice *inverts*
the model ranking. This adds the natural third option, the **integrated tumor burden**, the one
summary that sees *both* depth and tail, and shows it produces a third, distinct ranking — so even a
"comprehensive" metric is still a consequential choice.

- **New metric `log_burden_auc`** in the model-agnostic Stein/Bruno panel (`extract_tgi_metrics`) —
  the time-averaged log relative tumor size over the horizon (the AUC of the log-size curve, i.e. the
  log geometric-mean relative burden). Depth lowers it; a regrowth tail raises it. Eradication is
  floored at the detection limit (`v/y0 = 1e-3`, a complete response) so the integral is finite and
  stable, never `−∞`-dominated. Version-agnostic trapezoid (no numpy 2.x `trapezoid`/`trapz`
  dependency). It is horizon-dependent by construction — a cumulative-burden summary.
- **New record `survival_link.nsclc_os_burden_auc`** — a non-default (`default=false`) Weibull-PH OS
  link with `link_metric=log_burden_auc`, calibrated `beta`/`scale`, tier C, primary citation
  Wang (2009) (the registrational longitudinal-tumor-size→OS family the integrated burden belongs to).
  Reached only via `survival_link=`, so **every default view and export is byte-identical**.
- **The finding — a third distinct ranking, and a pathology repaired.** For NSCLC first line the three
  bridge metrics rank the five eligible models three different ways:
  - **week-8**: two-pop > acquired > Claret > Norton-Simon > Wang (the complete responder buried 4th;
    deep *early* shrinkers on top).
  - **k_g**: Norton-Simon > **Wang** > Claret > two-pop > acquired (the complete responder 1st — but
    the *minimal* responder Wang ranks 2nd on a slow regrowth slope, blind to the fact it never shrank).
  - **burden**: Norton-Simon > Claret > two-pop > acquired > **Wang** (the complete responder 1st, the
    deep-but-doomed resistance models demoted to the middle, and Wang dropped to last).
  The integrated burden agrees with neither extreme, and it **repairs k_g's depth-blindness**:
  tail-sensitivity without depth-sensitivity ranks a never-responding tumor second-best; an integrated
  metric ranks it last, where it belongs.
- **Budget V_link axis enriched.** NSCLC/first now has **4 eligible OS links** (week-8, Cox, `k_g`,
  burden-AUC), so the model-selection-budget's survival-link factor samples a fourth alternative.
- **11 closed-form + dynamics landmarks** (`tests/test_burden_auc.py`): baseline⇒0, constant⇒`log c`,
  floored eradication, horizon monotonicity, tail-sensitivity (where week-8 is blind), depth-sensitivity
  (where `k_g` is blind), the third distinct ranking, and the unchanged default view + transport→D.
- **Surfaces.** A three-metric slopegraph figure (`docs/images/burden_auc.png`), a CI-executed
  notebook (`notebooks/26_burden_auc_bridge_metric.ipynb`), a README section, and the
  metric-panel-keys API contract extended. No new module, kernel, or CLI command; pure post-processing
  over the existing kernels plus one record. 351 tests, 56 records.

## [0.32.0] — Acquired resistance: the resistance *origin* as a model-selection axis

Implements the research-track spec `docs/specs/research/acquired-resistance.md`. v0.24 made
the resistance *mechanism* a model-selection axis (phenomenological Claret λ vs a mechanistic
resistant subclone); that mechanistic model encodes one *origin* — a **pre-existing** resistant
clone. This adds the other canonical origin, **acquired** (drug-induced) resistance, and shows
the origin is a model-selection axis one layer below the mechanism — and one the week-8 OS
surrogate cannot see.

- **New `acquired_resistance` kernel** — a two-clone ODE (like two-population) where resistance
  is *generated by treatment*: `dS/dt = (kg − kd·E)·S − α·E·S`, `dR/dt = kgr·R + α·E·S`, with
  `R(0) = R0 = 0` (none pre-exists). The new parameter `α` is the drug-induced switching rate
  (Goldie-Coldman acquired-mutation kinetics). Setting `α = 0, R0 > 0` recovers the pre-existing
  two-population model exactly (strict generalization). Round-trip validated like every ODE
  kernel (SBML MathML, NONMEM `$THETA`, rxode2/Pumas — added to `ODE_RECORDS`).
- **New record `resistance.nsclc_first_line.acquired`** — `R0 = 0`, `α = 0.02`, with `kg`/`kd`/
  `kgr` **matched to the pre-existing two-population model**, so the only difference is the
  resistance origin. Tier C, illustrative, `unverified`; primary citation Foo-Michor (2014),
  *Evolution of acquired resistance to anti-cancer therapy* (the apt, DOI-verified anchor).
- **The finding — same week-8, different tail.** Matched on every shared parameter, the acquired
  and pre-existing origins **agree at week 8 and on the week-8-driven OS surrogate** (median OS
  92 vs 94 wk) but **diverge in the regrowth tail**: the acquired model has a markedly shallower
  nadir (8.0 vs 2.8 mm — the drug that kills the sensitive clone simultaneously *generates* the
  resistant one) and reaches RECIST progression earlier (mechanistic TTP 26 vs 32 wk). The
  resistance origin is a silent assumption the week-8 surrogate misses and the mechanistic PFS
  (v0.30) catches; adding the model raises NSCLC's PFS route-discordant pairs to 4/10.
- **No new module** — the acquired model is an ordinary TGI record; it joins the NSCLC compare
  set and is read by the existing `simulate` / `compare` / `pfs` / `response` / export paths
  unchanged. Guardrails ride through: worst-input-wins tier (C), out-of-context transport floors
  to D, population/trial level only, no therapy ranking; `α` is honestly flagged poorly
  identifiable (~110% CV, the acquired analog of the pre-existing `R0`).
- **Surfaces** — an `acquired_resistance` figure (clone decomposition + the origin divergence);
  `notebooks/25_acquired_resistance.ipynb` (CI-executed); 8 landmarks in `tests/test_acquired.py`
  (α=0 recovery, switching-flux conservation, no-drug-no-resistance, shallower nadir, earlier
  progression, week-8 OS agreement, tier/guardrails). Housekeeping: the duration-of-response
  figure's model-label lookup is now `.get`-robust so future model additions don't break it.
  README section + cheat sheets + research-track row + design-decision. Version 0.32.0
  (55 records).

## [0.31.0] — D-optimal trial design: can a better trial rescue the flat parameter?

Implements the research-track spec `docs/specs/research/optimal-design.md`. v0.22
(`onkos.identify`) evaluates whether a *given* sampling schedule can estimate a model's
parameters; this adds the schedule *choice* — the **D-optimal design** under a fixed
measurement budget — and answers what v0.22 left open: it separates the parameters a
better-designed trial could rescue (circumstantial flatness) from the parameters no trial
of this budget can pin down (structural flatness).

- **`onkos.design`** — new module, pure post-processing over the v0.22 Fisher-information
  core (no new kernel, record, or export). The design information `M = Σᵢ sᵢsᵢᵀ` is
  additive over timepoints, so the parameter sensitivities are computed once on a dense
  candidate grid and schedule optimization is pure linear algebra (no re-simulation per
  candidate).
- **`d_optimal_rows(scaled_sens, n, seed_rows)`** — the pure, landmark-tested selection
  core: greedy forward D-optimal row selection (maximize `log det M`), baseline-anchored.
- **`optimal_schedule(...)`** — binds the core to a record + a measurement budget, scoring
  the D-optimal schedule against a uniform schedule of the same budget. Reports the
  per-parameter Cramér-Rao RSE, the collinearity index, the **D-efficiency** over uniform
  (`(det M_opt / det M_uni)^(1/p)`), and whether the optimal design rescues a parameter
  uniform could not. The reported optimal is the better of greedy/uniform, so
  **D-efficiency ≥ 1 by construction**.
- **The finding** — for the Claret NSCLC model (N=7 over 48 wk) the D-optimal design
  concentrates samples at the kill phase (≈8–13 wk) and regrowth onset (≈30 wk + tail),
  improving every parameter (D-efficiency ≈ 1.14, collinearity γ 25→23). The two flat
  parameters **separate**: the borderline resistance term `λ` (54%) is **rescued** across
  the 50% identifiability line (48%), while the deeply flat growth rate `kL` (228%→199%) is
  only tightened — it stays **structurally unidentifiable** under the best possible
  schedule, so its huge CV is a flat-likelihood artifact, not biological spread. The
  2-parameter Wang biexponential is the control: both parameters identifiable under the
  optimal design (D-efficiency ≈ 1.31), proving the `kL` failure is the model's structure,
  not the optimizer.
- **Guardrails unchanged** — design analysis carries the record's propagated tier (cannot
  move it), is design/population-level only (never a per-patient schedule or a dosing/therapy
  choice), restricted to ODE records (a survival/transform record is not a trajectory to
  design for), and honest about the greedy heuristic (guaranteed only ≥ uniform).
- **Surfaces** — `onkos design <id>` (uniform-vs-optimal RSE table, optimal schedule,
  D-efficiency, structural-flat verdict); an `optimal_design` figure;
  `notebooks/24_optimal_design.ipynb` (executed in CI); 11 landmarks in
  `tests/test_design.py`; `optimal_schedule` added to the public-API contract. No new
  dataset record, kernel, or export. Version 0.31.0 (54 records).

## [0.30.0] — The PFS endpoint: two routes to progression-free survival

Implements the research-track spec `docs/specs/research/pfs-endpoint.md`. PFS is the
endpoint that gates accelerated and conditional approvals, and Onkos can compute it two
legitimate ways that need not agree. This opens the *mechanistic* route beside the existing
*statistical* one and shows the **route choice is a model-selection axis** — the model-
selection thesis now lives *inside a single endpoint*.

- **Mechanistic PFS** — `time_to_progression(t, v)` reads the RECIST 1.1 progression time
  directly off a tumor-size trajectory: the first time the SLD rises ≥20% above its
  *running nadir* (baseline included), `nan` if it never progresses within the horizon
  (right-censored). The progression arithmetic already lived inside DoR (v0.28); this
  promotes it to a first-class endpoint. No new kernel or record — pure post-processing.
- **Both routes, same trial** — `progression_free_survival(...)` returns, over the IIV
  ensemble, the mechanistic median TTP, a censoring-robust **landmark progression-free rate**
  (default 24 wk ≈ 6 months), the censored fraction, *and* the statistical median PFS from
  the context's parametric PFS link, plus their `route_ratio`. A `ttp_heavily_censored`
  warning fires above 50%, mirroring DoR.
- **The route is a model-selection axis** — `pfs_route_divergence(...)` counts the in-context
  model pairs whose mechanistic-PFS ranking contradicts their statistical-PFS ranking. The
  two-population (mechanistic resistance) model is the consistent culprit: **shortest or
  near-shortest mechanistic PFS, yet among the longest statistical PFS**, because at week 8
  it is deeply shrunk and the week-8-keyed hazard link reads that as durable benefit — blind
  to the resistant-clone regrowth the mechanism sees. The PFS endpoint's verdict on which
  model looks better depends on which route computes it (NSCLC 1L: 2/6 route-discordant
  pairs).
- **Dataset-wide on day one** — every solid-tumor context already has both a PFS link (v0.12)
  and a two-population model (v0.29), so the route inversion reproduces in all five contexts
  (NSCLC, breast, CRC, HCC, melanoma: ≥1 route-discordant pair each) — not an NSCLC artifact,
  unlike the v0.27 ORR surrogate which needed the `k_g` link.
- **Guardrails unchanged** — both routes carry the propagated tier (worst-input-wins; out-of-
  context transport floors to D), are trial-level only, never rank therapies, and surface
  censoring explicitly. Neither route is privileged; Onkos reports the *disagreement*.
- **Surfaces** — `onkos pfs <id>` (mechanistic TTP, landmark rate, statistical PFS, ratio)
  and `onkos pfs --routes` (the two-route table + route-discordance count); a `pfs_routes`
  figure; `notebooks/23_pfs_endpoint.ipynb` (executed in CI); 11 closed-form/relational
  landmarks in `tests/test_pfs_routes.py`; `progression_free_survival` and
  `pfs_route_divergence` added to the public-API contract. No new dataset record, kernel, or
  export. Version 0.30.0 (54 records).

## [0.29.0] — Cross-context generalization: the findings are not NSCLC artifacts

Implements the research-track spec `docs/specs/research/cross-context-generalization.md`,
a *breadth* release. The project's headline findings — the resistance-mechanism divergence
axis (v0.24), the model-selection budget's survival-link axis (v0.26), the conditional
ORR → OS surrogate (v0.27), and depth ≠ durability (v0.28) — were all demonstrated only on
NSCLC first line. This gives the four other curated solid-tumor contexts the same two pieces
NSCLC had and shows every finding reproduces.

- **8 new dataset records** — for breast, CRC, HCC, and melanoma first line: a mechanistic
  **two-population** resistance model (`resistance.<tt>_first_line.two_population`, kd
  matched to that context's Claret model, resistant growth `kgr` set faster than the
  sensitive growth → universal response but a fast resistant regrowth) and a non-default
  tail-sensitive **`k_g` OS link** (`survival_link.<tt>_os_growth_rate`, calibrated per
  context). Tier C, illustrative, `unverified`; existing citations. Each first-line context
  goes from a 2-model / 1-link layout to 3 models / 2 links.
- **Every NSCLC-only finding reproduces across all five contexts:** the two-population model
  has the highest ORR but the shortest DoR (depth ≠ durability); under the tail-sensitive
  `k_g` link it has the worst OS, so ORR — faithful under the week-8 surrogate (0 discordant
  pairs) — **mis-ranks** OS (2/3 or 3/3 discordant); and the budget's survival-link axis is
  now real (24–72% of forecast variance, the dominant axis in three of the four new
  contexts). The model-selection-budget report now flags **5 / 6 contexts as
  structure-dominated** (was 4/6).
- **CI-enforced generalization:** `tests/test_response.py` asserts the ORR-surrogate is
  concordant under week-8 and discordant under `k_g` in every first-line solid-tumor
  context; `tests/test_budget.py` asserts every such context has ≥2 OS links and a populated
  `v_link` (and the single-link control moved to NSCLC second line). A cross-context figure
  and `notebooks/22_cross_context_generalization.ipynb` (executed in CI).
- No new methodology, kernel, or module — pure curated content over the v0.24–v0.28 machinery.
  README states the findings as dataset-wide with a per-context table; the tumor-context
  library, budget, and divergence figures/numbers are refreshed. Version 0.29.0 (54 records).

## [0.28.0] — Duration of response: depth is not durability

Implements the research-track spec `docs/specs/research/duration-of-response.md`,
completing the response endpoint and isolating the *mechanism* of the v0.27 ORR → OS
surrogate failure. ORR measures response **breadth** (how many respond); this adds
**duration of response (DoR)** — response **durability** (how long).

- `response_episode(t, v)` returns the RECIST best response *and* the DoR from one
  observed-baseline trajectory, so the category and the duration are mutually consistent
  (`best_response` now delegates to it). DoR = time from PR onset (SLD ≤70% of baseline)
  to progression (SLD ≥120% of nadir); `nan` for a non-responder or a response that never
  progresses (right-censored).
- `objective_response_rate(...)` gains `median_dor_weeks`, `dor_censored_fraction`, and
  `n_responders` over the IIV ensemble — the population durability beside the population
  breadth — with a `dor_heavily_censored` warning when >50% of responders never progress
  (the observed median is then a lower bound, never mistaken for the truth). The durable
  responders are right-censored, not silently zeroed.
- **Depth is not durability:** the NSCLC model with the highest ORR (1.00, two-population)
  has the *shortest* median DoR (~32 wk) — its responses are universal but brief (a deep
  early shrink, then a fast resistant regrowth). That durability deficit is exactly why it
  has the worst tail-driven OS: sorted by survival under the k_g link, the *broadest*
  responder is the *worst* survivor while the longest-lived model's responses are the most
  durable. DoR is the mechanism behind the v0.27 surrogate failure and the immunotherapy
  lesson (a modest-ORR/durable drug beating a high-ORR/brief one on survival) made
  computable.
- Landmark-tested (`tests/test_duration.py`, 10 checks): episode consistency
  (`best_response == response_episode[0]`), nan DoR for non-responders, the closed-form
  DoR, censoring of durable responses, durability ordered by regrowth rate, the depth ≠
  durability dissociation, and that the k_g-discordant (highest-ORR) model is the short-DoR
  one (durability tracks survival where breadth inverts it).
- Surfaces: `onkos response` now reports DoR (with censoring), `onkos response --durability`
  prints the breadth-vs-durability table; a breadth-vs-durability figure and
  `notebooks/21_duration_of_response.ipynb` (executed in CI); README section + roadmap +
  cheat sheet + design-decisions row. Pure post-processing — no new kernel, dataset record,
  or export surface. Version 0.28.0.

## [0.27.0] — RECIST response & the ORR → OS surrogate

Implements the research-track spec `docs/specs/research/recist-orr-surrogate.md`: the
objective response rate (ORR) is the dominant phase-2 go/no-go endpoint, and a famously
contested OS surrogate. Onkos had OS and PFS but no response endpoint; this adds it and
makes the ORR → OS relationship a measured, context-dependent quantity.

- `onkos.response`: `best_response(t, v)` classifies RECIST 1.1 best overall response
  (`CR > PR > PD > SD`) from a tumor-size trajectory, measured from the observed baseline
  (CR ≥95% shrinkage, PR ≥30%, PD ≥20% regrowth from nadir with no PR, else SD).
  `objective_response_rate(...)` lifts it to **population rates** over the stored IIV
  ensemble — ORR = P(CR or PR), DCR = P(CR, PR, or SD), and the CR/PR/SD/PD distribution —
  with the median OS read off the *same* trial for the surrogate question. Trial-level
  only; no individual response probability.
- `response_vs_survival(...)` counts the **discordant model pairs** (one model with a
  higher ORR yet shorter OS than another), making the contested surrogate computable. The
  headline result is conditional: ORR and the week-8 survival surrogate are both
  shrinkage-based, so ORR ranks OS **perfectly under the week-8 link (0/6 discordant)** but
  **inverts under the tail-sensitive k_g link (4/6, 67%)** — the highest responder
  (ORR ≈ 1.0) has the shortest OS, the eradicating drug the longest. Whether ORR is a valid
  OS surrogate is conditional on the (unobservable) survival mechanism — the computational
  core of every "high response, no survival benefit" phase-3 failure.
- Landmark-tested (`tests/test_response.py`, 13 checks): the RECIST classification
  boundaries (exactly −30% → PR, ≥95% → CR, a PR that later regrows stays PR, no-PR +
  regrowth → PD), the rate simplex (`0 ≤ ORR ≤ DCR ≤ 1`, distribution sums to 1), ORR
  monotone in drug effect, the degenerate (zero-IIV) ensemble, tier passthrough / D-floor,
  and the conditional-surrogacy result (concordant under week-8, discordant under k_g).
- Surfaces: `onkos response <id> [--survival-link --surrogate --json]`; a RECIST-distribution
  + ORR-vs-OS figure and `notebooks/20_recist_response_orr.ipynb` (executed in CI); README
  section + roadmap + cheat sheet + design-decisions row; the public-API contract test
  extended. Pure post-processing — no new kernel, dataset record, or export surface.
  Version 0.27.0.

## [0.26.0] — The model-selection budget: structural variance decomposition (capstone)

Implements the research-track spec `docs/specs/research/model-selection-budget.md`, the
synthesis of the model-selection arc (v0.21–v0.25). Each structural choice had been scored
in isolation; this puts them on one ledger and names which assumption drives the forecast.

- `onkos.budget`: a balanced **two-way variance-component decomposition** (ANOVA /
  first-order Sobol over the structural factors — the structural analog of the parameter
  tornado in `onkos.sensitivity`). For a context it splits the total variance of median OS
  into `WITHIN(parameter) + V_model(TGI choice) + V_link(survival choice) + V_inter`, over
  the grid of in-context TGI models (`compare().included`) × every eligible survival link
  (week-8 Weibull, Cox, k_g). `variance_components(cell_means, cell_within)` is the pure,
  landmark-tested core; `model_selection_budget(...)` binds it to the forecast by reusing
  `ensemble_samples` per cell. Collapsing the survival-link factor to one level recovers
  exactly the v0.21 within/between split — the budget is a strict generalization.
- **The capstone finding:** for NSCLC first-line OS, ~68% of the forecast variance is
  irreducible structural-choice risk (only ~32% is parameter noise a bigger trial would
  shrink), and the single largest component is the **model×link interaction** — the v0.25
  result that the survival metric can *invert* which TGI model wins is exactly an
  interaction term, and it dominates. The survival-link axis (~24%) outweighs the
  tumor-growth-model axis (~12%).
- Landmark suite (`tests/test_budget.py`, 13 checks): the sum identity (components sum to
  total, fractions to 1), non-negativity of every component (incl. the residual
  interaction), single-factor collapse to the v0.21 split, identical-cells → zero between,
  pure-main-effect and additive-layout (zero interaction) cases, a match against a direct
  two-way sum-of-squares decomposition, the convex-hull bound, zero-within degeneracy, and
  the worst-included-tier floor; plus integration checks (the NSCLC four-component budget,
  the clinical-use flag, and the single-survival-link flag).
- Surfaces: `onkos budget [--tumor-type --line --endpoint --json]`; `onkos report` gains a
  per-context **model-selection budget** section ranking contexts by structure- vs
  parameter-dominance (binary, edge-safe at 0.5) and flagging contexts with only one
  survival link (the survival-model axis uncross-checked); a stacked-budget figure and
  `notebooks/19_model_selection_budget.ipynb` (executed in CI); README capstone section +
  roadmap + design-decisions row; the public-API contract test extended. Version 0.26.0.

## [0.25.0] — Survival-metric choice: which TGI metric predicts OS is a model-selection axis

Implements the research-track spec `docs/specs/research/survival-metric-choice.md`,
completing the v0.24 finding. A composed forecast is `TGI model → on-treatment metric →
survival link → OS`; the first and third arrows were already model-selection axes, but
the bridge metric was a silent constant (the week-8 change, a shrinkage surrogate blind
to the regrowth tail). This makes the bridge metric an explicit, swappable choice.

- `simulate` reads `structure.link_metric` from a survival link (default
  `"week8_relative_change"`, so **every existing curve is byte-identical**) and feeds that
  metric as the hazard covariate. A metric that did not occur (e.g. the growth-rate
  constant `k_g` for a tumor that never regrows) is `nan` and maps to the **no-effect
  covariate `x=0`** (the baseline hazard) — a complete responder gets the best, not an
  undefined, survival.
- New non-default record `survival_link.nsclc_os_growth_rate` (Weibull-PH,
  `link_metric: tumor_growth_rate_kg`, `default: false`): the tail-sensitive,
  more-prognostic Stein/Bruno growth-rate link, opt-in via `survival_link=` exactly like
  the Cox alternative (so *which metric* and *which baseline structure* are orthogonal
  survival-model choices). Citation `stein-2008-grc`; recorded external C-index 0.69 vs
  the week-8 link's 0.64 (k_g out-discriminates early shrinkage), but a better metric
  cannot raise a tier.
- **The sharp result — the metric choice inverts the answer** (illustrative NSCLC medians):
  under week-8 the mechanistic two-population resistance model looks *better* than the
  phenomenological Claret model (deeper early shrinkage, OS 94 vs 91); under `k_g` it looks
  *worse* (faster regrowth, OS 32 vs 39) — the v0.24 tail divergence made decisive, and
  pointing the other way. And the complete responder (Norton-Simon, eradication) is
  undervalued by week-8 (OS 58, last) but correctly ranked first by `k_g` (OS 102) — an
  early-shrinkage gate systematically penalizes a slow-but-complete responder.
- Landmark-tested (`tests/test_survival_metric.py`, 7 checks): backward compatibility (a
  link with no `link_metric` reads week-8; default curves unchanged), the metric materially
  moves OS, the resistance-model ranking inversion, the complete-responder re-ranking, the
  undefined-metric→baseline floor, non-default/opt-in behavior, and the D-floor on
  out-of-context transport.
- Surfaces: reachable through the existing `survival_link=` argument and the exports (no new
  module or kernel); a week-8-vs-`k_g` figure with the inversion and
  `notebooks/18_survival_metric_choice.ipynb` (executed in CI); README section + roadmap +
  design-decisions row; report counts regenerated. Version bump to 0.25.0.

## [0.24.0] — Mechanistic (two-population) resistance: the resistance model as a model-selection axis

Implements the research-track spec `docs/specs/research/mechanistic-resistance.md`:
resistance is the project's most load-bearing term (the λ "Hydra"), modeled until now
only *phenomenologically* (the Claret decay-of-drug-effect, whose λ has no cellular
referent and is ~90%-CV unidentifiable). This adds the *mechanistic* alternative and
makes the choice between them a quantified model-selection axis.

- New `two_population_resistance` reference kernel (Goldie-Coldman): a drug-**sensitive**
  clone and a pre-existing drug-**resistant** clone, observed together as `V = S + R`
  (`dS/dt = (kg − kd·E)·S`, `dR/dt = kgr·R`). The drug crushes the sensitive clone to a
  nadir; the untouched resistant clone then outgrows — the mechanistic origin of
  nadir-then-regrowth — and resistance is now a **biologically interpretable parameter**
  (`R0`, the initial resistant burden) instead of an unidentifiable rate. Two compartments,
  seeded via the existing multi-state `init_inputs` pattern; closed-form
  `V(t) = V0·e^{(kg−kd·E)t} + R0·e^{kgr·t}`.
- New dataset record `resistance.nsclc_first_line.two_population` (tier C), with the kill
  potency `kd` matched to the Claret model so the contrast isolates the **resistance
  mechanism** rather than an effect-scale difference. It joins the NSCLC first-line
  divergence view automatically — two resistance mechanisms in one context — raising the
  measured model-selection fraction **0.39 → 0.47** (the resistance-model axis is real
  between-model risk). Citations: Foo & Michor 2014 (modern quantitative dynamics, the
  record's primary) and Goldie-Coldman 1979 (the conceptual two-population origin, on `R0`).
- **The honest finding:** tuned to share the early kill, the two models agree at week 8
  (≈−87% vs −82%) and on the week-8-driven OS (median ≈94 vs ≈91 wk) yet diverge **≈5× in
  the tumor tail** (≈74 vs ≈15 mm at 3 yr) — the short-trial-indistinguishable,
  long-horizon-divergent failure mode — surfacing that a **week-8 OS surrogate is nearly
  blind to the resistance-model choice** (the silent-transport risk).
- Validated on both project axes: round-trip (SBML/NONMEM/rxode2/Pumas, both compartments,
  added to `ODE_RECORDS`) and a 10-check landmark suite (`tests/test_two_population.py`):
  the closed form, the `R0=0` reduction + eradication, the late-time log-slope → `kgr`
  (resistant-clone dominance), the interior nadir, resistant-fraction monotonicity, the
  untreated growth reduction, the divergence-axis integration, the D-floor on transport,
  and — composing with v0.22 — that `R0` stays practically unidentifiable (mechanistic
  does not mean measured).
- Surfaces: reachable through the existing `simulate` / `compare` / `identify` / exports
  (no new module); a clone-decomposition + resistance-divergence figure and
  `notebooks/17_mechanistic_resistance.ipynb` (executed in CI); README section, roadmap,
  kernel taxonomy, design-decisions row, and the refreshed NSCLC virtual-trial numbers.

## [0.23.0] — Drug-combination interaction: the interaction model as a model-selection axis

Implements the research-track spec `docs/specs/research/combination-interaction.md`
(steps 1–4): oncology is overwhelmingly combination therapy, and a composed forecast
for a combination silently depends on one unmeasured choice — *how do the two drugs'
effects combine?* This makes that choice a first-class, quantified model-selection axis.

- `onkos.interaction`: `combine_effects(E_A, E_B, model, psi)` combines two single-agent
  effect magnitudes under three declared interaction nulls — `hsa` (highest single
  agent, `max`), `additive` (Bliss-independence / effect-additive, `E_A + E_B`), and
  `greco` (interaction index `E_A + E_B + ψ·√(E_A·E_B)`, ψ>0 synergy, ψ<0 antagonism).
  `simulate_combination(...)` feeds the combined effect through the *existing* TGI →
  survival chain unchanged; `compare_interactions(...)` returns an
  `InteractionComparison` with the per-model OS/tumor trajectories and the
  **interaction-model divergence** — how much the predicted survival depends on the
  interaction assumption alone (≈77–100 wk median OS for the Claret NSCLC model at
  `E_A=E_B=0.6`, driven purely by the assumption).
- **Synergy is an assumption, not a finding.** `ψ` is a *declared* input (default 0, the
  additive null), never estimated from the dataset; a non-zero value carries a
  `synergy_is_an_assumption` warning — distinguishing synergy from additivity needs a
  combination trial designed for it. The underlying TGI model's propagated tier governs
  and cannot be raised; an inactive partner reduces to monotherapy under every model
  (no manufactured interaction).
- Guardrails, enforced by a landmark suite (`tests/test_interaction.py`, 13 closed-form
  checks of the combination rules themselves): the additive null (`greco(ψ=0)=additive`),
  the **Bliss≡additive identity** for log-linear kill (`1−(1−f_A)(1−f_B)=1−e^{−(E_A+E_B)}`),
  the `hsa ≤ additive ≤ greco(ψ>0)` ordering, monotonicity in ψ and in each effect,
  single-agent degeneracy, symmetry, the antagonism floor at 0, plus integration checks
  (divergence positive with synergy / zero with an inactive partner, tier passthrough,
  the synergy warning).
- Surfaces: `onkos interactions <id> [--effect-a --effect-b --psi --json]`; a combination
  figure (combined effect vs ψ + OS divergence by interaction model) and
  `notebooks/16_combination_interaction.ipynb` (executed in CI); README section framing
  the interaction model as a model-selection axis (the kill-mechanism move one layer up),
  cheat sheets, roadmap, and the public-API contract test updated.
- Methodological provenance (Bliss 1939; Loewe 1953 additivity and the dose-level
  boundary Onkos names but does not yet cross; HSA / Berenbaum 1989; the Greco 1995
  interaction index) is documented in the README and spec; Crossref-verified citation
  curation is deferred, consistent with the honest-by-default stance.

## [0.22.0] — Practical identifiability: could a trial even estimate this parameter?

Implements the research-track spec `docs/specs/research/practical-identifiability.md`
(steps 1–4): the project's load-bearing *claim* — that the ~90% IIV CV on
kill/resistance terms is there because "resistance is poorly identifiable from short
trials" — becomes a measured, tested quantity instead of an assertion.

- `onkos.identify`: from the **Fisher information of a clinical observation schedule**
  (`M = SᵀWS`, finite-difference sensitivities over the existing reference kernels,
  combined proportional+additive residual error), it returns the **Cramér–Rao** lower
  bound on each structural parameter's relative standard error (`rse_percent`) and the
  Brun–Reichert–Künsch **collinearity index** `γ_K`. `onkos.identifiability(...)`
  returns an `Identifiability` with the per-parameter predicted RSE paired with the
  stored IIV CV, the `practically_identifiable` verdict (`max RSE < 50%` AND
  `γ_K < 15`), and the least-identifiable parameter (curation triage).
- The headline pairing — **predicted RSE next to stored IIV CV** — operationalizes the
  thesis: for the Claret NSCLC model under a realistic RECIST cadence the kill rate
  `kD` is well identified (RSE ≈ 9%) while the growth rate `kL` and resistance decay
  `lambda` are flat (RSE ≈ 229%, 53%) and confounded (`γ_K ≈ 22`), so `lambda`'s 96%
  CV is flagged `cv_is_identifiability_artifact` — partly a flat-likelihood artifact of
  the originating design, not a clean estimate of biological spread. Identifiability is
  design-relative: lengthen the follow-up past resistance-driven regrowth and `lambda`
  finally crosses below the ceiling.
- Guardrails, enforced by a landmark suite (`tests/test_identifiability.py`, 14
  closed-form checks of the information algebra itself): the exponential single-
  parameter CRLB closed form `RSE_k = σ_prop/(|k|·√Σtᵢ²)`, Fisher-information
  additivity over observations, monotonic precision (more scans never raise an RSE),
  residual-error scaling, singular-design honesty (`inf`, never a fabricated bound),
  the orthogonal-design `γ_K = 1` floor and its scale-invariance, and CRLB
  consistency. Identifiability **cannot move a tier** (it passes the record's
  propagated tier through; an out-of-context transport still floors to D), emits no
  individual-level quantity, and is the individual (fixed-effects) design FIM — not the
  population/NLME FIM, stated explicitly to avoid overclaim.
- Surfaces: `onkos identify <id> [--schedule --sigma-prop --json]`; `onkos report`
  gains a per-model **practical-identifiability** section ranking the clinical TGI
  models a realistic design cannot support (the 2-parameter biexponential models are
  identifiable; every 3-parameter resistance-augmented model is not) — binned to a
  binary verdict so the CI report-in-sync diff stays byte-stable; an identifiability
  figure (RSE-vs-CV bars + RSE-vs-follow-up curves) and `notebooks/15_practical_
  identifiability.ipynb` (executed in CI).
- Methodological provenance (the pharmacometric optimal-design FIM, PFIM/PopED; the
  structural-vs-practical identifiability distinction of Raue et al. 2009; the Brun et
  al. 2001 collinearity index; the Cramér–Rao bound) is documented in the README and
  spec; Crossref-verified citation curation is deferred, consistent with the honest-by-
  default stance.

## [0.21.0] — Model-selection uncertainty: the third uncertainty axis

Implements the research-track spec `docs/specs/research/model-selection-uncertainty.md`
(steps 1–4): the virtual-trial divergence view gains its inferential completion —
from *"the models disagree by this much"* to a variance decomposition plus an
honestly-weighted central forecast that carries its disagreement.

- `onkos.combine`: a new module that splits a composed survival forecast's total
  predictive uncertainty into **within-model** (parameter/IIV noise, Axis 1) and
  **between-model** (model-selection risk, Axis 3) via the law of total variance,
  and reports the headline `model_selection_fraction = BETWEEN / (WITHIN + BETWEEN)`
  — the fraction of forecast uncertainty that more data on any *one* model cannot
  resolve. `Comparison.model_average(...)` returns a `ModelAverage` with the
  averaged OS/PFS curve `S̄(t)`, its pointwise between-model band, the worst
  included tier, and the weights.
- Three **declared** weighting schemes (`equal`, `tier` A:B:C = 4:2:1, `evidence`
  ∝ external C-index − 0.5), the cross-scheme `weight_sensitivity` swing with a
  fragility warning, and — wherever weights appear — the explicit label that these
  are *forecast-combination weights (Bates–Granger), NOT posterior model
  probabilities* (the models are fit to different trials, so a posterior model
  probability is not identifiable and would be invented).
- Guardrails, enforced by a landmark suite (`tests/test_combine.py`, 16 closed-form
  checks of the combination math itself — law of total variance to ≤1e-9, equal-
  weight identity, identical-component zero-between, convex-hull bound, survival-
  function validity, monotone re-weighting, zero-weight inertness, tier floor):
  averaging **cannot raise a tier**; only in-context models are averaged; a
  single-eligible-model context yields fraction 0 *and* a `single_eligible_model`
  warning; the point estimate is structurally inseparable from its fraction.
- Surfaces: `onkos compare --average [--weights …] [--decompose] [--json]`;
  `Comparison.to_dict(model_average=…)` embeds an optional `model_average` block in
  the virtual-trial JSON with the `onkos:modelSelectionUncertainty` predicate (added
  to the JSON-LD context); `onkos report` ranks clinical contexts by irreducible
  model-choice risk (curation triage — where adding a better-validated model has the
  most value); a model-averaging figure and `notebooks/14_model_averaging.ipynb`
  (executed in CI).
- `onkos.uncertainty.ensemble_samples` factored out of `simulate_ensemble` as the
  shared per-model sampling core, so the combiner reuses the Axis-1 machinery for
  `WITHIN` rather than re-implementing it.
- Methodological provenance (Bates–Granger forecast combination; the MCP-Mod and
  NLME model-averaging regulatory precedent; the BMA common-data boundary Onkos
  deliberately does not cross) is documented in the README; Crossref-verified
  citation curation is deferred, consistent with the honest-by-default stance.

## [0.20.0] — Scientific landmark validation (a second validation axis)

- `tests/test_landmarks.py`: validates every reference kernel against the
  characteristic, analytically-derivable property of the *published model* it
  implements — Gompertz inflection at `Vmax/e`, the Simeoni tumor-static
  concentration `c* = λ0/k2`, the Norton-Simon stationary condition `E = g/k`,
  the bi-exponential nadir time, the Weibull median and proportional-hazards
  identity, Emax half-maximal effect at `EC50`, IO immune homeostasis `s/δ`,
  and more (15 checks across all 15 kernels).
- This is a second validation axis independent of the export round-trip: the
  round-trip proves exports agree with the kernel; landmarks prove the kernel
  *is* the model it names. A kernel can be self-consistent yet wrong — landmarks
  catch that.
- `docs/validation-landmarks.md`: catalogues each landmark, its closed-form
  condition, and its source. README documents the two-axis validation strategy.
- This is the honest reading of spec §9's "compare against published example
  simulations": the landmark is the published property, derived from the model's
  own equations — no digitized data is fabricated.

## [0.19.0] — Architecture as a tested contract + refreshed system design

- `tests/test_architecture.py`: pins the system structure so docs and code cannot
  silently drift — every declared subsystem has records (would have caught the
  empty `drug_effect`), the CLI export formats match the builders *and* the CI
  export sweep (would have caught `so`/`jsonld` missing from the loop), kernel
  kinds are known, no kernel is dead/orphaned, and the public API surface is stable.
- README architecture refreshed to the current 18-version system: a layered
  diagram (data → core → kernels → analyses → exports → presentation), a kernel
  taxonomy (ODE / survival / transform), and an updated round-trip data-flow.

## [0.18.0] — Norton-Simon kill model (fills the drug_effect subsystem)

- `norton_simon` kernel: the Norton-Simon hypothesis on a Gompertz growth law,
  `dV/dt = (g - k·E)·V·ln(Vmax/V)` — drug kill proportional to the GROWTH rate, so
  a smaller, faster-growing tumor is more chemo-sensitive. Closed form + rhs +
  per-state SBML round-trip; the scale-robust analytic-vs-ODE metric handles the
  Gompertz collapse toward zero.
- `drug_effect.norton_simon.nsclc` record fills the previously-empty `drug_effect`
  subsystem (spec §3) and joins the NSCLC divergence view as a third, distinct
  model — so the assumed *kill mechanism* (log-kill vs growth-proportional), not
  just the parameters, is now a visible model-selection axis. Norton 2005 citation.
- Kill-mechanism figure + notebook 13. Zero tier inflations after the add.

## [0.17.0] — PK composability bridge (dose → exposure → tumor → survival)

- `onkos.pk`: a small, illustrative PK bridge that realizes the spec's headline
  Hypnos-composability claim self-contained. `steady_state_metrics` returns the
  standard exposure metrics (C_avg, C_max, C_min, AUC_tau) for a one-compartment
  oral regimen (cornerstone `C_avg = F·Dose/(CL·tau)`); `concentration_profile`
  builds the multiple-dose Bateman curve; `from_profile` ingests an external
  (Hypnos) concentration-time profile onto the simulation grid.
- The full dose → C_avg → exposure-response → kill → tumor → OS/PFS chain now runs
  end to end: higher dose → more exposure → deeper response → longer OS.
- Onkos still does NOT model PK (that is Hypnos); the generators are clearly
  labelled illustrative, and `from_profile` is the real ingestion path.
- Composability-chain figure + notebook 12.

## [0.16.0] — Line of therapy + line-aware survival matching

- Fix: `_find_survival_links` matched only on tumor type, so a second-line context
  silently borrowed first-line survival models. It now matches on
  `(tumor_type, line)`; an unsupported line yields no survival curve, not a wrong
  one (mirroring the no-fallback rule for tumor type).
- NSCLC second-line context added: baseline (more advanced disease), OS + PFS
  survival links (shorter than first line), and a second-line Claret TGI (faster
  growth, weaker kill, faster resistance). The first-line-only Claret model is
  correctly excluded from the 2L divergence view; 2L survival is shorter than 1L.
- First-vs-second-line figure + notebook 11; zero tier inflations after the add.

## [0.15.0] — Evidence-based tier audit ("tiers are partly numeric")

- `onkos.audit` / `onkos audit`: derives the best confidence tier a clinical TGI /
  survival record's recorded evidence supports (external C-index → A/B; a
  poorly-identified kill/resistance term with IIV CV ≥ 70% → tier-C ceiling) and
  flags any record whose assigned tier *exceeds* its evidence (tier inflation).
- The inflation check runs inside `onkos validate`, so an over-claimed tier fails
  CI and cannot regress. The shipped dataset has zero inflations.
- The health report gains an "Evidence-based tier audit" section (inflations +
  conservative/upgrade-candidate counts).
- Operationalizes spec §5 ("tier assignment is partly numeric") and §9
  (predictive_performance feeds tier) — the honesty thesis applied to the tier
  field itself.

## [0.14.0] — Serializable virtual-trial result + rebuilt dashboard

- `Comparison.to_dict()` / `.to_json(include_curves=...)`: the virtual-trial
  divergence view is now a JSON-serializable result (per-model OS/PFS medians and
  TGI metrics, excluded models + reasons, OS/PFS divergence, optional curves) for
  dashboard / external-simulator ingestion. `onkos simulate --compare --json`.
- The Streamlit dashboard is rebuilt as a thin renderer over the tested package
  API: divergence view (tumor/OS/PFS), an analyze-a-model tab (uncertainty bands +
  sensitivity tornado), and a dataset browser. Its data logic is unit-tested and
  CI lints + compiles `dashboard/app.py` against the current API.
- Housekeeping: silence the external rdflib `ConjunctiveGraph` DeprecationWarning
  in pytest; lint scope now includes `dashboard/`.

## [0.13.0] — Cox survival link + survival-model-choice divergence

- `survival_cox_ph` kernel: Cox proportional hazards with a NONPARAMETRIC
  tabulated baseline survival `S0(t)` (`S(t|x) = S0(t)^exp(beta*x)`), the feature
  that distinguishes Cox from the parametric Weibull form. Completes the spec's
  "Cox and parametric OS/PFS link models".
- The Cox link (`survival_link.nsclc_os_cox`) is marked `structure.default: false`
  so it never auto-collides with the Weibull OS link; selecting it via
  `survival_link=` enables a "Weibull vs Cox" survival-model-choice comparison —
  a third uncertainty axis (median OS ~91 vs ~107 wk on the same TGI metric).
- Kernel framework: `uses_baseline` flag + `structure.baseline_survival` table
  (schema), injected by `simulate`. The baseline rides along in vt-json / JSON-LD.
- Survival-model-choice figure + notebook 10.

## [0.12.0] — Progression-free survival (PFS) endpoint

- Second survival endpoint: a PFS link per tumor context (parametric Weibull-PH on
  the week-8 TGI metric), distinguished by a `structure.endpoint` tag (OS|PFS).
  `simulate` returns a curve per endpoint (`Trajectory.survival`, `median_pfs`,
  `pfs_curve`); PFS is shorter than OS by construction.
- `compare` reports `pfs_divergence` / `median_pfs_range`; `sensitivity` accepts
  `target="median_pfs_weeks"`; the CLI shows OS and PFS side by side.
- 5 PFS survival-link records (NSCLC, breast, CRC, HCC, melanoma); the existing OS
  links are tagged `endpoint: "OS"`; schema gains the `endpoint` enum.
- OS-vs-PFS figure + notebook 09. Completes the spec's "OS/PFS link models".

## [0.11.0] — Linked data (JSON-LD / RDF)

- JSON-LD export (`onkos export --format jsonld`, `to_jsonld`, `dataset_jsonld`):
  records render as RDF — confidence tier, clinical-use prohibition, derivation
  context, transportability, and `bqbiol:isDescribedBy` DOI/PMID links become
  resolvable triples. The single `@context` shipped in
  `dataset/schema/context.jsonld` (previously unused) is now the source of truth.
- The virtual-trial JSON gains an `@context` + `@id`, making it valid JSON-LD.
- The COMBINE `.omex` now includes the JSON-LD document.
- Validated for real: tests expand the output with rdflib and assert the expected
  triples appear (added `rdflib` to the dev/test extra).

## [0.10.0] — Export-layer completion: PharmML SO + cross-format integrity

- PharmML **Standard Output (SO)** exporter (`onkos export --format so`): MLE
  population estimates, inter-individual variability as random-effect variances
  (`omega = ln(1+CV²)`, never mislabeled as estimate precision), external-validation
  diagnostics, and the universal Onkos annotations. Completes spec §7's
  "PharmML (+ SO)" and the `.omex` "SBML + PharmML + SO + provenance" bundle.
- The COMBINE `.omex` now includes the SO with its MIME type in the manifest.
- Round-trip hardening: rxode2 and Pumas parameter vectors are now re-parsed and
  checked, and a cross-format consistency test asserts NONMEM, SBML, PharmML-SO,
  rxode2, and Pumas all agree on the parameter values — one source of truth.

## [0.9.0] — Variance-based parameter sensitivity

- `onkos.sensitivity` / `onkos sensitivity`: ranks a record's IIV-bearing
  parameters by how much their variability drives a target (median OS or any
  metric). Independent lognormal sampling makes the standardized regression
  coefficient equal the input-target correlation, so squared SRCs give a
  first-order variance decomposition (with direction).
- Turns uncertainty bands into curation triage: it shows *which* parameter to
  verify first (for the Claret NSCLC model, the kill rate kD drives ~90% of the
  median-OS variance — more than the higher-CV resistance term).
- Tornado figure + notebook 08.

## [0.8.0] — TGI-metric extraction (Stein/Bruno panel)

- `onkos.metrics.extract_tgi_metrics`: model-agnostic extraction of the derived
  TGI metrics from any trajectory — depth of response, nadir / time-to-growth,
  the tumor growth-rate constant **k_g**, the shrinkage-rate constant **k_s**
  (via `k_s = k_g − s0`), and the RECIST-style **duration of response**.
- Recovers the biexponential kernel's generating k_g and k_s to within ~10% and
  the Claret growth constant k_L — a built-in correctness check; inapplicable
  metrics are returned as `nan`, never fabricated.
- Wired into `simulate()` (all prior metric keys preserved); the uncertainty
  ensemble now summarizes each metric over its finite samples.
- Annotated TGI-metric figure + notebook 07.

## [0.7.0] — Parameter-uncertainty propagation

- `onkos.simulate_ensemble` / `onkos uncertainty`: Monte-Carlo propagation of the
  stored inter-individual variability (`iiv_cv_percent`). Parameters with an IIV
  CV are sampled lognormally (median-preserving) and tumor-size, TGI-metric, and
  population-OS distributions are returned as median + percentile bands — closing
  the gap where IIV was stored but never used in simulation.
- `simulate(..., param_overrides=...)` for per-run parameter substitution.
- Uncertainty figure + notebook 06; design-decision and README sections framing
  this as the second uncertainty axis (parameter variability) alongside the
  model-selection divergence view.

## [0.6.0] — Phase F: hardening (releasable, self-reporting)

- `onkos report` — machine-generated dataset health & validation report
  (`docs/dataset-health.md`), kept in sync with the data by a CI/test gate.
- External-validation backfill: every eligible clinical TGI / survival model now
  carries an external-validation metric (`predictive_performance`) — coverage
  15/15.
- Releasability proven in CI: the wheel is built with the bundled dataset,
  installed into a clean environment, and exercised from outside the repository.
- Fixed dataset resolution precedence so the source `dataset/` always shadows a
  stale synced `_dataset/` during development (bundled copy is the wheel-only
  fallback).
- Release metadata: `.zenodo.json`, `CHANGELOG.md`, `py.typed`.
- Coverage figure (tier × subsystem heatmap + external-validation gauge).

## [0.5.0] — Phase E: immuno-oncology (hypothesis-tier)

- Kuznetsov 1994 tumor-immune QSP kernel (`io_tumor_immune`); control / escape /
  checkpoint-rescue dynamics.
- Two hypothesis-tier records (tier D), enforced by the validator; every export
  carries `DO NOT USE FOR PREDICTION`; excluded from the clinical view.
- Multi-initial-condition support (`init_inputs`); per-state initial amounts in
  SBML/NONMEM exports.

## [0.4.0] — Phase D: preclinical translation

- Multi-state ODE framework (optional analytic, `observable`, system integrator).
- Simeoni 2004 xenograft model (`simeoni_exp_linear`, `simeoni_tgi` with the
  signal-distribution transit chain) and in-vitro → in-vivo potency translation.
- Per-state SBML/NONMEM export + round-trip; preclinical excluded from the
  clinical view; xenograft-on-human floors to D.

## [0.3.0] — Phase C: survival + baselines

- `tumor_type_baselines` library and per-context Weibull-PH survival links across
  NSCLC, breast, CRC, HCC, melanoma; ≥2 eligible TGI models per context.
- Cross-context divergence; no-orphan-record invariant enforced in CI.

## [0.2.0] — Phase B: resistance + exposure-response

- Emax / sigmoid-Emax / power exposure-response kernels driving the kill term.
- Scalar and time-varying (PK-profile) drug effect; Hypnos composability.
- ER tier + transportability propagation; centralized version constant.

## [0.1.0] — Phase A: TGI spine

- Dataset-as-source-of-truth (JSON Schema + JSON-LD), the Claret 2009 TGI model,
  growth laws, NSCLC context, TGI→OS survival link, the virtual-trial divergence
  view, NONMEM + SBML export, and round-trip validation.
