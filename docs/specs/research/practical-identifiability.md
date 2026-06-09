# Onkos — research spec: practical identifiability & trial-design information

**Status:** implemented in v0.22.0 (`onkos.identify`; steps 1–4 of §10). This is the
design-of-record; the methodological source anchors of §7 are documented but their
Crossref-verified citation curation is still pending, honest by design. Written in
the v0.1 house style; every parameter value is illustrative and `unverified` by
design — the infrastructure is the contribution.

**The question under the uncertainty.** Onkos already surfaces that the kill and
resistance terms carry coefficients of variation near 90% (spec §1, §4). The stated
reason is that *resistance is poorly identifiable from short trials.* But the project
asserts that reason without measuring it. This spec measures it: given a model and a
realistic clinical observation schedule, compute — from the **Fisher information of the
design** — the best precision with which each structural parameter *could* be estimated,
and flag the parameters a trial of that shape cannot pin down at all. The headline output
is the per-parameter predicted relative standard error (RSE) and a single
`practically_identifiable` verdict, with the resistance term's near-non-identifiability
made a number instead of an anecdote.

> The IIV CV that travels on a parameter answers *"how much does this quantity vary across
> patients?"* Practical identifiability answers a prior and different question: *"could the
> trial that produced this estimate have located the parameter at all, or is the 90% CV a
> symptom of a flat likelihood?"* A parameter can be both biologically variable **and**
> structurally unidentifiable; conflating the two is exactly how a fitted-but-meaningless
> resistance rate gets transported into the next program's control stream. This is the
> design-level companion to the parameter-level sensitivity triage (`onkos.sensitivity`)
> and the model-level model-selection triage (`onkos.combine`).

---

## 1. The problem this extends

Onkos quantifies how uncertainty *propagates* and how it *splits*. It does not yet
quantify whether the uncertainty was ever *resolvable* by the trial that reported it.

| Axis | Question | Today |
| --- | --- | --- |
| **Parameter (within-model)** | How far does the forecast move under the reported IIV CV? | ✅ `onkos.uncertainty` — lognormal Monte-Carlo bands. |
| **Sensitivity** | *Which* parameter's uncertainty drives the survival prediction (verify it first)? | ✅ `onkos.sensitivity` — variance attribution. |
| **Model-selection (between-model)** | How much does the answer depend on which eligible model you pick? | ✅ `onkos.combine` — law-of-total-variance split + model average. |
| **Practical identifiability (design)** | Given a realistic sampling schedule, could the parameter have been **estimated** at all — or is its CV a flat-likelihood artifact? | ⚠️ **Not measured.** The spec *asserts* resistance is poorly identifiable; nothing computes it. |

Sensitivity asks which uncertainty *matters*; identifiability asks which uncertainty is
*reducible by better data versus structurally stuck*. They are complementary: a parameter
can be high-sensitivity and well-identifiable (collect more of the same data and the
forecast tightens), or — the dangerous case the resistance term embodies — high-sensitivity
and **un**-identifiable (more data of the same shape will never tighten it; only a different
design, or an external constraint, can).

**Why this is the right deepening (and the right scope).** It (1) turns a load-bearing
*claim* of the project into a *measured, tested quantity*; (2) is pure post-processing over
the existing reference kernels — finite-difference sensitivities re-using
`simulate(..., param_overrides=...)`, no new ODE kernel, no dataset subsystem, near-zero
schema change; (3) is the canonical tool of pharmacometric **optimal design** (the Fisher
information matrix; PFIM/PopED), so it has direct regulatory-science precedent (§7); (4) is
*safe by construction* — it characterizes a *model-and-design*, never a patient, and emits
no individual prediction and no therapy ranking; and (5) sharpens the honest message: the
output most worth printing is *"this trial design cannot estimate this parameter,"* which is
the opposite of false precision.

---

## 2. The statistical framework

Fix a record's TGI/growth kernel with structural parameters `θ = (θ_1, …, θ_p)` at their
curated central values `θ*`, a clinical **observation schedule** `t_obs = (t_1, …, t_n)`
(the scan times), and a **residual-error model** for the tumor-size measurement. The
prediction at the central values is `y_i = f(t_i; θ*)` (the reference-kernel tumor-size
trajectory).

**Residual error.** Tumor-size assays are dominated by proportional error; Onkos uses the
standard combined model `σ_i = sqrt(σ_add² + (σ_prop · y_i)²)`, default `σ_prop = 0.2`
(20% CV), `σ_add = 0`. This is a *declared* measurement model, printed with every result —
not a fitted quantity.

**Local sensitivity matrix** `S ∈ ℝ^{n×p}`, `S_ij = ∂f(t_i; θ)/∂θ_j |_{θ*}`, computed by
central finite difference (relative step `h`, default `1e-4`) over the reference kernel —
the same perturb-and-resimulate machinery `onkos.sensitivity` uses, here read structurally
rather than stochastically.

**Fisher information matrix (the design information):**

```
M(θ*, t_obs) = Sᵀ W S,        W = diag(1/σ_i²)            (p × p)
```

This is the **individual, fixed-effects** design FIM — the information a single
representative tumor-size profile, sampled on `t_obs` with this residual error, carries about
the structural parameters under independent Gaussian residuals. (It is deliberately *not* the
full population/NLME FIM, which would require linearizing the random-effects structure; that
is a heavier, separable extension noted in §8. The individual FIM is the standard local
practical-identifiability diagnostic and is exactly what answers "can this short trial locate
λ?".)

**Cramér–Rao precision bound.** Where `M` is invertible, the parameter covariance is bounded
below by `C = M⁻¹`, giving the best achievable standard errors and the headline relative
standard error:

```
SE_j = sqrt(C_jj),     RSE_j = SE_j / |θ_j*|        (the predicted precision, in %)
```

A parameter with `RSE_j` above a declared ceiling (default 50%) is **practically
unidentifiable** under this design: no estimator can do better, so the reported point value
is not supported by a trial of that shape.

**Collinearity index (structural near-non-identifiability).** RSE inflates both when a
parameter is individually flat *and* when two parameters are confounded (their sensitivity
columns nearly parallel — e.g. a growth rate and a kill rate that trade off over a short
horizon). The standard scale-free diagnostic (Brun et al. 2001) normalizes each weighted
sensitivity column to unit length and takes

```
γ_K = 1 / sqrt( λ_min( S̃ᵀ S̃ ) ),     S̃_j = (W^{1/2} S_j) / ‖W^{1/2} S_j‖
```

`γ_K = 1` ⇔ orthogonal (perfectly separable directions); `γ_K → ∞` ⇔ collinear
(non-identifiable *combination*). The default ceiling is `γ_K < 15`. Reporting both RSE and
`γ_K` distinguishes *"this parameter is individually flat"* from *"this parameter is
confounded with another"* — different design fixes (denser early sampling vs. an external
constraint).

**The verdict** is the conjunction, so neither failure mode hides the other:

```
practically_identifiable  =  (max_j RSE_j < rse_ceiling)  AND  (γ_K < collinearity_ceiling)
```

---

## 3. The thesis tie — CV is not identifiability

The single most important output is the **side-by-side of stored `iiv_cv_percent` and
predicted `rse_percent`** for each parameter. The dataset's defining honesty move is
surfacing the ~90% CV on kill/resistance terms; this spec explains *where that CV comes
from*. A parameter that is **both high-CV and high-RSE** earns an explicit
`cv_is_identifiability_artifact` warning: its large reported variability is, at least in
part, a flat-likelihood artifact of the originating trial design, not a clean estimate of
biological spread — precisely the parameter whose out-of-context transport (Onkos's
load-bearing risk) is least defensible. This is the curation-triage payload: it ranks *which
parameter most needs a richer design or an external constraint before its estimate should be
trusted or reused.*

---

## 4. Tier & guardrail propagation (unchanged invariants)

Identifiability analysis is a diagnostic *about* a record under a design; it never changes
the record:

- **It cannot raise (or lower) a tier.** The result carries the record's own propagated tier
  for context; a well-identified parameter does not upgrade a C-tier record, and a flat one
  does not downgrade it. Identifiability and validation are orthogonal axes (a model can be
  externally validated yet have one locally-flat nuisance parameter, and vice versa).
- **It is design-relative and says so.** Every result names its `t_obs` and residual-error
  model; "unidentifiable" always means "under *this* design," never "in principle."
- **Population/design level only.** `f(t; θ*)` is the reference trajectory of a *published
  model*, not a patient's tumor; RSE is a property of the *estimator under a design*, not a
  prediction for anyone. No individual output is added.
- **No false precision in the other direction either.** A singular FIM yields `RSE = ∞` and
  `γ_K = ∞` (reported as such), never a silently regularized finite number — the honest
  statement is "the design does not identify this," not a fabricated bound.

---

## 5. Reference kernel & validation landmarks

No new ODE kernel — the analyzer is finite-difference sensitivity plus linear algebra over
existing trajectories. It gets its own **landmark suite** (in the spirit of
`tests/test_landmarks.py` and `tests/test_combine.py`): closed-form properties of the
information algebra itself, so the implementation is provably the estimator it claims to be.

| Landmark | Closed-form condition |
| --- | --- |
| **Exponential closed form** | For `y = V0·e^{kt}` with proportional error, `RSE_k = σ_prop / (|k|·sqrt(Σ_i t_i²))` exactly (the single-parameter CRLB). |
| **Information additivity** | `M(t_A ∪ t_B) = M(t_A) + M(t_B)` — Fisher information adds over independent observations. |
| **Monotonic information** | Adding an observation never raises any `RSE_j` (added information is positive-semidefinite). |
| **Residual-error scaling** | Scaling every `σ_i` by `c` scales every `RSE_j` by `c` (FIM ∝ 1/σ², cov ∝ σ²). |
| **Structural non-identifiability** | Two identical (parallel) sensitivity columns ⇒ singular FIM ⇒ `RSE = ∞`, `γ_K = ∞`, `practically_identifiable = False`. |
| **Orthogonal design** | Orthogonal weighted sensitivity columns ⇒ `γ_K = 1` (the lower bound). |
| **Collinearity scale-invariance** | Multiplying one parameter's sensitivity column by a constant leaves `γ_K` unchanged (it is the normalized diagnostic). |
| **Collinearity floor** | `γ_K ≥ 1` for any design. |
| **CRLB consistency** | `RSE_j` equals `sqrt(diag((SᵀWS)⁻¹))_j / |θ_j|` computed independently. |
| **Tier passthrough** | The result's tier equals the record's propagated tier (analysis cannot move it). |

This mirrors the project's two existing validation axes: round-trip proves exports match the
kernel; landmarks prove a kernel *is* the model it names; here the landmarks prove the
analyzer *is* the Fisher-information CRLB and the Brun collinearity index — not an
unconstrained precision guess.

---

## 6. API, CLI, and surface

**Python.** A new `onkos.identify` module with a pure, landmark-tested core
(`fisher_information`, `crlb_rse`, `collinearity_index`) and a record-facing
`identifiability(...)` returning an `Identifiability` dataclass:

```python
res = onkos.identifiability(
    ds, "resistance.claret_2009.tgi",
    context=dict(tumor_type="NSCLC", line="first"),
    schedule=[0, 6, 12, 18, 24, 36, 48],   # weeks — a realistic RECIST scan cadence
    sigma_prop=0.2,                         # declared residual error (20% CV)
)

res.practically_identifiable      # bool — the headline verdict (under this design)
res.collinearity_index            # γ_K
res.params                        # per-parameter: symbol, rse_percent, iiv_cv_percent, identifiable
res.worst                         # the least-identifiable parameter (curation triage)
res.tier                          # the record's propagated tier (unchanged)
res.warnings                      # cv_is_identifiability_artifact / singular_fim / …
res.to_dict()                     # carries the clinical-use prohibition, like every result
```

**CLI.**

```bash
onkos identify resistance.claret_2009.tgi --schedule 0,6,12,18,24,36,48 --sigma-prop 0.2
onkos identify resistance.claret_2009.tgi --json
```

prints, per parameter, the predicted RSE next to the stored IIV CV, the collinearity index,
the verdict, and any flat-likelihood-artifact flags.

**Report.** `onkos report` gains a per-model **practical-identifiability** section under a
fixed reference design: each clinical TGI model's least-identifiable parameter and a binned
verdict (identifiable / borderline / unidentifiable), so contributors can see *where a richer
trial design or an external constraint is needed before an estimate should be reused* — the
design-level analog of the sensitivity and model-selection triage. (Binned qualitatively,
like the model-selection-uncertainty summary, so the CI report-in-sync diff stays
byte-stable.)

**No new export model.** Identifiability is an analysis *over* a model under a design, not a
model; it adds no NONMEM/SBML/PharmML surface. (A future option, noted in §8, is an
`onkos:practicalIdentifiability` RDF predicate carrying the verdict under the reference
design — deferred until the reference design itself is something the project wants to pin in
the dataset rather than in code.)

---

## 7. Source anchors (methodological; DOIs added at curation time)

Well-established methods, not Onkos parameters; added to `dataset/citations/` through the
normal Crossref/PubMed-verified pipeline, honest until a human confirms each.

- **Fisher information for nonlinear design.** Mentré, Mallet & Baccar (1997), *Optimal
  design in random-effects regression models*, Biometrika — the population-FIM foundation the
  individual FIM here specializes.
- **The pharmacometric optimal-design tools.** Bazzoli, Retout & Mentré (2010), *PFIM*; and
  the PopED line — the established software embodiment of FIM-based design evaluation Onkos's
  analyzer mirrors one layer up (TGI/survival).
- **Practical vs. structural identifiability.** Raue et al. (2009), *Structural and practical
  identifiability analysis … using the profile likelihood*, Bioinformatics — the distinction
  (a parameter can be structurally identifiable yet practically flat under finite, noisy data)
  that this spec operationalizes via the CRLB.
- **Collinearity index.** Brun, Reichert & Künsch (2001), *Practical identifiability analysis
  of large environmental simulation models*, Water Resources Research — the source of the
  normalized-sensitivity collinearity index `γ_K` used here.
- **Cramér–Rao bound.** Standard; the lower bound on estimator variance that makes the
  inverse FIM a *best-case* precision, so "unidentifiable here" is a conservative verdict.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not the full population/NLME FIM.** Onkos computes the *individual* (fixed-effects) design
  FIM. The population FIM (information on the fixed effects *and* the IIV variances jointly,
  via a first-order linearization) is a real and separable extension; conflating the two would
  overclaim. The individual FIM is sufficient — and correct — for the question "can a short
  trial locate λ?".
- **Not optimal design.** The analyzer *evaluates* a given design; it does not *optimize* the
  schedule. Recommending sampling times edges toward prescription and is out of scope for v0.x
  (and the optimum is design-objective-dependent). Evaluation is descriptive; optimization is
  a future, clearly-bounded direction.
- **No individual identifiability.** "Identifiable" is always a property of the
  model-and-design, never of a patient.

---

## 9. Safety & scope (unchanged hard line)

- **Population / design level only.** Everything here is a property of a *published model under
  a stated trial design*. It is **not** an estimate of any person's parameters, tumor, or
  survival, and adds no individual-level output.
- **No therapy ranking.** The analyzer characterizes one model under one design; it never
  compares drugs or regimens.
- **Design-relative honesty.** "Unidentifiable" means "under this schedule and residual-error
  model," and the result says so; it is never read as an in-principle claim.
- **No false precision in either direction.** A flat or singular design yields `∞`, reported
  as `∞`; the analyzer never regularizes a fabricated finite bound, and it never tightens a
  tier.
- **The line, restated.** Any feature that takes a real patient's measurement and returns an
  estimate or a therapy choice **does not get built.** Measuring whether a trial design could
  identify a model parameter changes none of this; it makes the project's own
  poorly-identifiable claim legible, and stops there.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Information core** | `onkos.identify` with the FIM, the CRLB RSE, and the collinearity index over the existing kernels; the `Identifiability` result with per-parameter RSE-vs-CV. Landmark suite (§5). | RSE computes for the Claret λ term and every landmark passes; the resistance term is measurably the least identifiable. |
| **2 — Verdict & artifact flagging** | `practically_identifiable`, the `worst` parameter, and the `cv_is_identifiability_artifact` flag tying high CV to high RSE; singular-FIM honesty (`∞`, not a fabricated bound). | A parameter that is both high-CV and flat is flagged; a singular design reports `∞` and `practically_identifiable = False`. |
| **3 — CLI, figure, notebook** | `onkos identify`; an identifiability figure (RSE-vs-CV bars + an information-vs-schedule-length curve showing λ stays flat as the schedule lengthens); a CI-executed notebook; the clinical-use flag on `to_dict`. | The diagnostic is reachable from the CLI and visualized; the notebook runs in CI. |
| **4 — Report & curation triage** | `onkos report` ranks clinical TGI models by their least-identifiable parameter under a fixed reference design (binned), the design-level analog of the sensitivity/model-selection triage. | `onkos report` flags the models whose estimates a realistic trial cannot support, byte-stably. |

Step 1 alone is a self-contained, citable contribution: an open, validated tool that takes a
curated oncology TGI model and a realistic scan schedule and reports, with a Cramér–Rao
guarantee, which parameters that trial could not have estimated — turning the field's
qualitative *"resistance is poorly identifiable"* into a tested number, is something nobody
ships openly as an artifact.
