# Changelog

All notable changes to Onkos are documented here. Versions follow the phased
roadmap (spec §11). All parameter values are illustrative and `unverified` by
design; the infrastructure is real and tested.

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
