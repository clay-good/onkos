# Changelog

All notable changes to Onkos are documented here. Versions follow the phased
roadmap (spec §11). All parameter values are illustrative and `unverified` by
design; the infrastructure is real and tested.

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
