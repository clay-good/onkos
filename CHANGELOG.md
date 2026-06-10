# Changelog

All notable changes to Onkos are documented here. Versions follow the phased
roadmap (spec §11). All parameter values are illustrative and `unverified` by
design; the infrastructure is real and tested.

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
