# Onkos — research spec: D-optimal trial design (can a better-designed trial rescue the flat parameter?)

**Status:** implemented in v0.31.0 (`onkos.design`). This is the design-of-record; the
methodological source anchors of §7 are documented but their Crossref-verified citation
curation is still pending, honest by design. Written in the v0.1 house style; every value is
illustrative and `unverified` by design — the infrastructure is the contribution.

**v0.22 asked whether a given trial *can* estimate a parameter; this asks for the *best* trial
that could, and whether even it is enough.** `onkos.identify` (v0.22) computes, from the design
Fisher information of a *fixed* sampling schedule, the Cramér-Rao precision (RSE) each
structural parameter could be estimated to — and found that for the Claret NSCLC model the
kill rate `kD` is well-identified (~9% RSE) while the growth rate `kL` and the resistance term
`λ` are flat (≈228% / 54% RSE) and confounded (collinearity γ≈25). That left the obvious next
question open: **was the schedule just badly chosen?** A pharmacometrician designing a trial
does not accept a fixed grid; they *choose* the sampling times. This spec adds that choice —
the **D-optimal sampling schedule** under a fixed budget of N measurements — and answers the
question v0.22 could not: it separates the parameters a *better trial* could rescue (the
flatness was circumstantial) from the parameters *no trial of this budget* can pin down (the
flatness is structural).

> The distinction is the whole point. "Your resistance term has a 96% CV" can mean two very
> different things: *you ran the wrong trial* (fixable — design a better one) or *no trial of
> this size can pin it down* (structural — the parameter is not estimable and the CV is a
> flat-likelihood artifact, exactly as v0.22 flagged). Optimal design is what separates them,
> and for the Claret model it separates them *within one model*: the D-optimal schedule rescues
> the borderline resistance term but leaves the deeply flat growth rate unidentifiable. Onkos
> computes the best schedule a fixed budget allows and reports which parameters it moves across
> the identifiability line — turning "badly designed trial" vs "structurally unidentifiable
> parameter" into a decided question rather than an excuse.

---

## 1. The problem this extends

The design Fisher information `M = S̃ᵀS̃` (v0.22) is a **sum of rank-one contributions, one per
sampling time**: `M = Σᵢ s̃ᵢ s̃ᵢᵀ`, where `s̃ᵢ` is the residual-weighted parameter-sensitivity
row at observation time `tᵢ`. Two consequences drive this spec:

| Consequence | What it enables |
| --- | --- |
| The information is **additive over timepoints** | The sensitivity rows can be computed once on a dense candidate grid; then *any* schedule's Fisher information is just the sum of its chosen rows — schedule optimization is pure linear algebra, **no re-simulation per candidate**. |
| The information depends on **which** times are chosen | Different schedules give different `M`, hence different parameter precision. The schedule is a *design variable*, and there is a best one. |

v0.22 evaluated one schedule. This spec searches the schedule space for the **D-optimal**
design — the N sampling times that maximize `det(M)`, i.e. minimize the volume of the joint
confidence ellipsoid of the parameters (the standard optimal-design criterion). It reuses the
v0.22 information core unchanged (`fisher_information`, `crlb_rse`, `collinearity_index`); the
only new machinery is the **selection of rows**.

**Why this is the right deepening (and the right scope).** It (1) completes the v0.22 arc —
v0.22 *evaluates* a design, this *optimizes* it, the move the reserve list flagged as "optimal
trial design (v0.22 stopped at evaluation)"; (2) is pure post-processing — it reuses the
existing kernels and the existing Fisher core, no new record, kernel, or export (mirrors
`identify`/`sensitivity`/`budget`); (3) is rigorous and standard — D-optimality is textbook
optimal-design theory, the FIM is the same one v0.22 ships; (4) is *safe by construction* — a
design-level quantity over published model structures, no patient data, no individual
prediction, no therapy ranking; (5) sharpens the honest message: the best schedule a fixed
budget allows still cannot rescue the resistance term, so its huge CV is **structural**, not a
design artifact — a strictly stronger, safer claim than v0.22's.

---

## 2. The D-optimal schedule, defined

For a dynamic (ODE) TGI record, a context, and a measurement budget of `N` samples over a
horizon `[0, H]`:

```
candidate grid   = a dense set of admissible measurement times in [0, H]
s̃ᵢ              = residual-weighted sensitivity row at grid time tᵢ   (v0.22 error model)
D-optimal design = argmax over N-subsets 𝒮 of   log det( Σ_{i∈𝒮} s̃ᵢ s̃ᵢᵀ )
```

The N-subset search is combinatorial, so Onkos uses **greedy forward selection** (start from the
mandatory baseline `t=0`, repeatedly add the time that most increases `log det M`) — a standard,
deterministic D-optimal heuristic. The reported optimal is then taken as the **better of the
greedy design and the uniform design**, which *guarantees the reported optimal is never worse
than uniform* (`D-efficiency ≥ 1` by construction). The baseline measurement is mandatory
(`include_baseline=True`), as in any real trial.

The design is scored against a **uniform** schedule of the same budget `N` over the same horizon
— the design a non-optimizing trial would use — with the standard summary:

```
D-efficiency = ( det M_optimal / det M_uniform )^(1/p)        (p = # parameters; ≥ 1)
per-parameter RSE, collinearity index γ, identifiable set     (reused from v0.22)
```

`D-efficiency = 1.3` means the optimal design extracts the information of a uniform design 1.3×
its size — a real, free precision gain from sampling times alone.

---

## 3. The result — design helps, but cannot rescue a flat parameter

For the Claret NSCLC model, a budget of `N = 7` over a 48-week horizon (illustrative, RSE %):

| Parameter | uniform RSE | D-optimal RSE | what the design did |
| --- | --- | --- | --- |
| `kD` (kill rate) | 9 | 9 | identifiable under both — already pinned |
| `λ` (resistance/regrowth) | 54 | **48** | **rescued** — the better design crosses the 50% line |
| `kL` (growth rate) | 228 | 199 | tightened but **still flat** — structurally unidentifiable |

The D-optimal design concentrates samples at the **kill phase** (≈8–13 wk) and the **regrowth
onset** (≈30 wk and the horizon tail) — the times where the parameters are most separable — and
it improves *every* parameter (D-efficiency ≈ 1.14, collinearity γ 25→22). Crucially, the two
flat parameters separate: the **borderline** resistance term `λ` (54%, just over the line) is
**rescued** by the better schedule (48%, under the line), while the **deeply** flat growth rate
`kL` (228%) is only tightened to 199% — it stays unidentifiable under the *optimal* schedule.
That is the finding, and it is sharper than v0.22's: optimal design is not a uniform "yes" or
"no" — it rescues the flatness that was *circumstantial* (a fixable design problem) and exposes
the flatness that is *structural* (`kL` cannot be pinned down by any schedule of this budget, so
its huge CV is a flat-likelihood artifact, not biological spread). The `λ` rescue sits near the
ceiling and is sensitive to the budget and grid; the `kL` verdict is robust across them — the
landmark suite pins the robust one.

Contrast a model the design *can* serve — the 2-parameter Wang biexponential (NSCLC):

| Parameter | uniform RSE | D-optimal RSE | identifiable? |
| --- | --- | --- | --- |
| `kg` (growth) | 12 | 10 | ✅ both |
| `ks` (shrinkage) | 34 | 28 | ✅ both |

Here both parameters are identifiable and D-optimal design tightens them further (D-efficiency
≈ 1.3). The method works when the model allows it — which is exactly why its *failure* on the
3-parameter resistance model is informative rather than a tooling artifact: **optimal design is
the control that proves the resistance term is structurally, not circumstantially,
unidentifiable.** This is the rigorous capstone to v0.22: it ranks which clinical TGI models a
realistic, *optimally designed* trial can support (the 2-parameter biexponential, yes; the
3-parameter resistance models, no), separating a design problem from a model problem.

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **Cannot move a tier.** A design analysis carries the record's propagated tier and never
  changes it (identical to `identify`/`sensitivity`/`budget`).
- **Design / population level only.** The output is a sampling schedule and a precision
  prediction for *estimating a published model's parameters* — never an individual prediction,
  never a dosing or treatment choice.
- **No therapy ranking.** It ranks *designs* (and, derivatively, which model a design can
  support), never treatments.
- **Honest about the heuristic.** Greedy forward selection is a heuristic, not a proof of global
  optimality; the reported design is guaranteed only `≥ uniform`. This is stated, and the
  D-efficiency makes the realized gain explicit.

---

## 5. Validation landmarks

No new kernel — the design core is row-selection over the v0.22 Fisher information. The landmark
suite (`tests/test_design.py`) pins the behavior:

| Landmark | Condition |
| --- | --- |
| **Closed-form selection** | on a hand-built scaled-sensitivity matrix whose D-optimal `N`-subset is known by construction (orthogonal high-magnitude rows), `d_optimal_rows` returns it. |
| **Additivity** | the Fisher information of a schedule equals the sum of its rows' outer products — the property the optimizer relies on. |
| **D-efficiency ≥ 1** | the reported optimal design's `det M` is never below the uniform design's (the guarantee from taking the better of greedy/uniform). |
| **Budget & baseline** | the returned schedule has exactly `N` times, all within `[0, H]`, and includes the baseline when `include_baseline=True`. |
| **Structural flat survives the best design** | for Claret NSCLC the optimal design lowers every RSE and the collinearity index, yet the deeply flat growth rate `kL` stays unidentifiable (`'kL' in structurally_flat`) — a structural flat the best schedule cannot remove (the robust verdict, not the budget-sensitive `λ` rescue). |
| **Method works when the model allows** | for the 2-parameter Wang biexponential both parameters are identifiable under the optimal design (`structurally_flat == []`) — the control that the `kL` failure above is the model's, not the tool's. |
| **Tier passthrough & guardrails** | the propagated tier rides through unchanged; the clinical-use flag is on every payload; a non-ODE kernel raises (a survival/transform record is not a trajectory to design for). |

---

## 6. API, CLI, and surface

```python
od = onkos.optimal_schedule(ds, "resistance.claret_2009.tgi",
                            context={"tumor_type": "NSCLC", "line": "first"},
                            n_samples=7, horizon=48.0)
od.optimal.schedule            # the D-optimal sampling times (weeks), baseline-anchored
od.optimal.rse_percent         # {"kL": 199, "kD": 9, "lambda": 48} — Cramér-Rao precision
od.d_efficiency                # ~1.14 — how much more informative than uniform
od.rescues_any                 # True — the better design rescued the borderline λ
od.structurally_flat           # ["kL"] — deeply flat even under the best schedule

# the pure, landmark-tested core: D-optimal row selection over a scaled-sensitivity matrix
from onkos.design import d_optimal_rows
rows = d_optimal_rows(scaled_sens, n=7, seed_rows=(0,))   # baseline-anchored greedy
```

**CLI.** `onkos design <id>` prints the uniform-vs-D-optimal RSE table, the optimal schedule,
the D-efficiency, and the structural-flat verdict.

**No new export model** — a design is an analysis over a record's kernel, not a model.

---

## 7. Source anchors (methodological; DOIs added at curation time)

- **Optimal design of experiments.** The standard treatment of D-/A-/E-optimality and the
  Fisher-information design criterion (e.g. Atkinson, Donev & Tobias, *Optimum Experimental
  Designs*; Fedorov's exchange algorithm) — the basis for maximizing `det(M)` over sampling
  schedules.
- **Optimal design in pharmacometrics.** The population-optimal-design literature (e.g.
  PopED / PFIM tooling and the Mentré school) applying the FIM to PK/PD sampling-time selection
  — the field this mirrors at the design (not population-random-effects) level.
- **Practical identifiability.** Already cited for v0.22 (`identify`): the structural-vs-
  practical identifiability distinction that optimal design here makes operational.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not population-optimal design.** The FIM here is the individual/design FIM v0.22 ships, not
  the population FIM with random-effects blocks (the PopED/PFIM object). A population extension
  is a clean later step.
- **Not a global optimum.** Greedy forward selection is a heuristic guaranteed only `≥ uniform`;
  a Fedorov exchange refinement or an exact search is a noted follow-on.
- **Not dose / regimen optimization.** It optimizes *measurement times* for parameter precision,
  never doses or treatments.
- **Not A-/E-optimality (yet).** D-optimality (det) is the shipped criterion; A-optimality
  (trace of `M⁻¹`) and E-optimality are clean additions on the same core.
- **No new record.** Design is computed over the existing TGI records.

---

## 9. Safety & scope (unchanged hard line)

- **Design / population level only.** The output is a sampling schedule and a parameter-precision
  prediction for a *published model*, never an individual prognosis.
- **No therapy ranking, no dosing.** It ranks designs, never treatments or doses.
- **Cannot raise a tier.** A design analysis carries the worst tier of the record unchanged.
- **Honest about its limits.** The greedy heuristic, the fixed budget, and the design-FIM scope
  are stated; the realized gain is reported as D-efficiency, never overclaimed.
- **The line, restated.** Any feature that takes a real patient's data and returns a sampling
  schedule *for that patient* or a treatment choice **does not get built.** Making the best trial
  design under a fixed budget computable changes none of this.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Selection core** | `d_optimal_rows(scaled_sens, n, seed_rows)` — greedy forward D-optimal row selection over a scaled-sensitivity matrix. | The closed-form and additivity landmarks pass. |
| **2 — Bind to a record** | `optimal_schedule(...)` builds the dense-grid sensitivities (reusing the v0.22 machinery), runs the core, and scores optimal vs uniform (RSE, collinearity, D-efficiency). | The Claret design table and the `D-efficiency ≥ 1` guarantee hold. |
| **3 — The structural verdict** | `rescues_any` / `structurally_flat` — does the optimal design move any parameter across the identifiability line? | The "design helps but cannot rescue" finding and the biexponential control are shown. |
| **4 — Surfaces** | `onkos design` CLI; a uniform-vs-optimal figure + a CI-executed notebook; README section. | The design is reachable, visualized, and documented with the project's guardrails. |

Step 1–2 alone is a self-contained contribution: an open, landmark-tested D-optimal sampling
designer over published TGI models, reusing the existing Fisher-information core. Step 3 is the
payload: showing that the *best* schedule a fixed budget allows still cannot identify the
resistance term — that its huge CV is structural, not a design failure — is the rigorous
capstone to the v0.22 identifiability work, shipped as a tested artifact.
