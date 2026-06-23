# Spec completion status

**Every Onkos spec is implemented, tested, and shipped.** This index is the
human-readable companion to [`tests/test_specs_complete.py`](../../tests/test_specs_complete.py),
which enforces the spec ↔ implementation mapping below as a CI invariant: a
research spec written without an implementation (or an implementation removed out
from under a spec) fails the suite. "All specs are complete" is therefore a
checked condition, not a claim.

Two spec layers exist:

- **`docs/specs/v0.1/spec.md`** — the foundational design spec (the envelope,
  subsystems, kernels, exports, tiers, safety line, and the phased roadmap).
- **`docs/specs/research/*.md`** — one self-contained forward-looking spec per
  research direction, each shipped as a tagged version.

## v0.1 design spec — phased roadmap (§11)

| Phase | Content | Status |
| --- | --- | --- |
| A — TGI spine | growth laws + Claret TGI + NSCLC context + TGI→OS link + divergence view; NONMEM/SBML; round-trip | ✅ v0.2 |
| B — Resistance + exposure-response | λ resistance, exposure-response, PharmML/rxode2/Pumas, IIV-CV | ✅ v0.2 |
| C — Survival + baselines | TGI-metric→survival breadth + `tumor_type_baselines` library (5 contexts) | ✅ v0.3 |
| D — Preclinical translation | Simeoni xenograft + in-vitro→in-vivo potency | ✅ v0.4 |
| E — Immuno-oncology (hypothesis-tier) | tumor–immune QSP, shipped non-predictive | ✅ v0.5 |
| F — Hardening | external-validation scaffold, COMBINE `.omex`, Zenodo, CITATION.cff | ✅ v0.6 |

**Enumerated requirement sets** (pinned by `test_specs_complete.py`):

- §2 growth-law family — exponential, logistic, Gompertz, Simeoni exp→linear,
  von Bertalanffy, power-law — **6/6 kernels shipped**.
- §6 reference kernels — growth, kill/transit, resistance, exposure-response,
  TGI-metric extraction, survival — all present and round-trip validated.
- §7 export layer — NONMEM · SBML · PharmML · SO · rxode2 · Pumas · virtual-trial
  JSON · JSON-LD · COMBINE `.omex` · CSV · BibTeX — **11/11 formats shipped**.

## Research specs (`docs/specs/research/`)

All 21 shipped. Each row is enforced by `SPEC_INDEX` in `test_specs_complete.py`.

| Spec | Version | Primary artifact | Landmark test |
| --- | --- | --- | --- |
| model-selection-uncertainty | v0.21 | `onkos.combine` | `test_combine.py` |
| practical-identifiability | v0.22 | `onkos.identify` | `test_identifiability.py` |
| combination-interaction | v0.23 | `onkos.interaction` | `test_interaction.py` |
| mechanistic-resistance | v0.24 | kernel `two_population_resistance` | `test_two_population.py` |
| survival-metric-choice | v0.25 | `onkos.simulate` (`link_metric`) | `test_survival_metric.py` |
| model-selection-budget | v0.26 | `onkos.budget` | `test_budget.py` |
| recist-orr-surrogate | v0.27 | `onkos.response` | `test_response.py` |
| duration-of-response | v0.28 | `onkos.response` | `test_duration.py` |
| cross-context-generalization | v0.29 | `onkos.compare` (5 contexts) | `test_context_library.py` |
| pfs-endpoint | v0.30 | `onkos.response` | `test_pfs_routes.py` |
| optimal-design | v0.31 | `onkos.design` | `test_design.py` |
| acquired-resistance | v0.32 | kernel `acquired_resistance` | `test_acquired.py` |
| burden-auc-bridge-metric | v0.33 | `onkos.metrics` (`log_burden_auc`) | `test_burden_auc.py` |
| joint-survival | v0.34 | `onkos.joint` | `test_joint.py` |
| loewe-additivity | v0.35 | `onkos.interaction` (Loewe) | `test_loewe.py` |
| exposure-response-extrapolation | v0.36 | `onkos.dose_response` | `test_dose_response.py` |
| early-surrogate-timing | v0.37 | `onkos.early_surrogate` | `test_early_surrogate.py` |
| model-discriminability | v0.38 | `onkos.discriminability` | `test_discriminability.py` |
| model-selection-atlas | v0.39 | `onkos.atlas` | `test_atlas.py` |
| von-bertalanffy-growth | v0.40 | kernel `growth_von_bertalanffy` | `test_landmarks.py` |
| power-law-growth | v0.41 | kernel `growth_power_law` | `test_landmarks.py` |

## Out of scope (declared, not pending)

Per spec §2/§10 these are **deliberate exclusions**, not unfinished specs:
per-patient prognosis / treatment recommendation (the hard safety line),
mechanistic intracellular signaling, hematologic-malignancy dynamics (a later
enumerated extension, not a v0.x item), and radiotherapy/surgery dynamics.

## The one genuinely open track: data verification (not a spec gap)

The *infrastructure* specs are complete; what remains is **scientific
verification of the dataset's parameter values** — a curation activity, not an
unimplemented spec. Its mechanism (`review_status: pending_human_review`), its
evidence dossier ([`docs/verification/workhorse-models.md`](../verification/workhorse-models.md)),
and its queue (`onkos review-queue`) are themselves built and shipped; promoting
records to `verified` is gated on human PDF sign-off per
[`CONTRIBUTING.md`](../../CONTRIBUTING.md).
