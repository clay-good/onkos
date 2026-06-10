# Onkos — research spec: survival-metric choice (which TGI metric predicts OS)

**Status:** implemented in v0.25.0 (`link_metric`-configurable survival links + a
growth-rate-constant OS link). This is the design-of-record; the methodological source
anchors of §7 are documented but their Crossref-verified citation curation is still
pending, honest by design. Written in the v0.1 house style; every value is illustrative
and `unverified` by design — the infrastructure is the contribution.

**Completing the v0.24 finding.** The mechanistic-resistance work (v0.24) ended on an
uncomfortable note: the phenomenological and mechanistic resistance models diverge ~5× in
the tumor *tail* yet agree on overall survival — because the survival link reads only the
**week-8 change**, a shrinkage surrogate blind to the regrowth tail. That is not a quirk of
those two models; it is a property of the *metric the survival link consumes*. This spec
makes that metric an explicit, swappable choice and adds the **growth-rate-constant (k_g)**
link — the tail-sensitive, more-prognostic Stein/Bruno quantity — turning *"which
on-treatment metric predicts survival"* into a first-class model-selection axis. The result
is sharp: switching the metric can **re-rank, even invert, which model looks better.**

> Onkos already treats *which TGI model* (the divergence view) and *which survival-model
> structure* (Weibull vs Cox) as model-selection axes. The metric that bridges them — the
> single on-treatment number fed into the hazard — has been a silent constant (week-8
> change). It is arguably the most consequential choice of all: the entire go/no-go
> apparatus of early oncology rests on the claim that an early tumor-size endpoint predicts
> survival, and the field actively debates *which* endpoint. Onkos makes the choice visible
> and shows that it changes the answer.

---

## 1. The problem this extends

A composed Onkos forecast is `TGI model → on-treatment metric → survival link → OS`. The
first and third arrows are already model-selection axes; the **middle arrow was fixed**:

| Bridge metric | Captures | Status before this spec |
| --- | --- | --- |
| **week-8 relative change** | early shrinkage depth | ✅ the hardcoded covariate of every link |
| **growth-rate constant k_g** | post-nadir regrowth rate (the tail) | ⚠️ extracted (`tgi_metrics`) but never used to drive survival |
| (landmark size, AUC, TTP, …) | other facets of the trajectory | future |

The week-8 surrogate is the standard early endpoint, but it is, by construction, blind to
what happens after week 8 — exactly the regrowth dynamics that resistance models disagree
about and that the Stein/Bruno literature finds *most* prognostic for OS. Fixing the bridge
metric silently imports that blindness into every survival forecast.

**Why this is the right deepening (and the right scope).** It (1) completes the project's
own freshly-surfaced v0.24 finding rather than adding breadth; (2) is a near-zero-code change
plus one dataset record — the metric becomes a declared field (`structure.link_metric`),
defaulting to the existing week-8 behavior, so nothing regresses; (3) has direct precedent —
the tumor growth-rate constant as the dominant OS predictor is the central result of the
Stein/Bruno program (§7); (4) is *safe by construction* — population/trial-level survival
under a declared metric, no individual prognosis, no therapy ranking; and (5) sharpens the
honest message to its strongest form: *the choice of surrogate metric can invert which
treatment-model looks better, and the field has been assuming one.*

---

## 2. The mechanism

A survival link declares which on-treatment metric drives its hazard:

```
x = metrics[ link.structure.link_metric ]          # default: "week8_relative_change"
S(t) = exp( −(t/scale)^shape · exp(beta · x) )      # Weibull proportional hazards
```

The default is unchanged (`week8_relative_change`), so every existing link and curve is
byte-identical. A link may instead declare `link_metric: "tumor_growth_rate_kg"` — the
post-nadir log-linear regrowth rate — with a positive `beta` (faster regrowth ⇒ higher
hazard). A metric that **did not occur** (e.g. `k_g` for a tumor that never regrows) is
`nan`, which maps to the **no-effect covariate** `x = 0` (the baseline hazard) — so a
complete responder gets the best, not an undefined, survival.

The k_g link ships **non-default** (`structure.default = false`), reached via the explicit
`survival_link=` argument exactly like the Cox alternative, so it never collides with the
default week-8 link on the same endpoint. *Which metric* and *which baseline structure* are
two orthogonal, opt-in survival-model choices.

---

## 3. The result — metric choice re-ranks, and inverts

For the NSCLC first-line models, the two metrics tell different stories (illustrative
medians, weeks):

| Model | week-8 change | OS (week-8 link) | k_g | OS (k_g link) |
| --- | --- | --- | --- | --- |
| Claret (phenomenological resistance) | −0.82 | 91 | 0.021 | 39 |
| Two-population (mechanistic resistance) | −0.87 | **94** | 0.025 | **32** |
| Norton-Simon (complete responder) | −0.29 | 58 | — (none) | **102** |
| Wang biexponential | −0.20 | 54 | 0.020 | 44 |

Two headline inversions fall out:

1. **The resistance-model ranking flips.** Under week-8 the *mechanistic* model looks
   better (deeper early shrinkage, OS 94 > 91); under k_g it looks worse (faster regrowth,
   OS 32 < 39). The v0.24 tail divergence, invisible to week-8, becomes the deciding factor —
   *and points the other way.* Which survival metric you assume determines which resistance
   model wins.
2. **The complete responder is undervalued by the surrogate.** Norton-Simon eradicates the
   tumor (slow early shrinkage, no regrowth). The week-8 surrogate scores it mediocre (58,
   below both resistance models); the k_g metric — seeing no regrowth — correctly makes it
   the longest survivor (102). An early-shrinkage endpoint systematically penalizes a
   slow-but-complete responder, which is precisely the molecule class an early go/no-go gate
   should not discard.

This is the model-selection thesis at its most consequential: not "the models disagree," but
"the *endpoint you chose to compare them* disagrees, and it can reverse your decision."

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **The link is a contributing component; worst-input-wins still governs.** The k_g link is
  tier C; a forecast through it inherits the worst tier of its chain, and an out-of-context
  transport floors to D exactly as for the week-8 link.
- **Higher discrimination is recorded, not assumed.** The k_g link carries a higher external
  C-index (illustratively 0.69 vs 0.64) because k_g out-discriminates early shrinkage in the
  Stein/Bruno data — but it is still `unverified`, and a better metric does **not** raise the
  composed tier.
- **No individual output, no ranking.** The curves are trial-level. The re-ranking is a
  statement about *models under a metric*, never a recommendation between treatments for a
  patient.
- **The default is sacred.** Making the metric configurable changes no existing curve; the
  week-8 behavior is the default and is covered by a backward-compatibility landmark.

---

## 5. Validation landmarks

No new kernel — the survival kernel is unchanged; the addition is a declared covariate
source. The landmark suite (`tests/test_survival_metric.py`) pins the behavior:

| Landmark | Condition |
| --- | --- |
| **Backward compatibility** | a link with no `link_metric` reads `week8_relative_change`; default OS curves are unchanged. |
| **Metric actually moves OS** | the k_g link gives a materially different median OS than the week-8 link. |
| **Resistance-ranking inversion** | `OS(two-pop) > OS(claret)` under week-8 **and** `OS(two-pop) < OS(claret)` under k_g. |
| **Complete-responder re-ranking** | `OS(norton) < OS(claret)` under week-8 **and** `OS(norton) > OS(claret)` under k_g. |
| **Undefined-metric floor** | `k_g = nan` (no regrowth) ⇒ covariate 0 ⇒ a finite, baseline (best-case) survival, never a nan curve. |
| **Non-default + opt-in** | the k_g link has `default = false` and is absent from the default OS divergence; reached only via `survival_link=`. |
| **Transport** | used out of context, the k_g link floors the composed tier to D with a warning. |

---

## 6. API, CLI, and surface

No new module or CLI verb — the metric choice flows through the existing `survival_link=`
argument and the dataset:

```python
# Default (week-8 surrogate) vs the tail-sensitive growth-rate metric.
default = onkos.simulate(ds, "resistance.nsclc_first_line.two_population", context=ctx)
kg = onkos.simulate(ds, "resistance.nsclc_first_line.two_population", context=ctx,
                    survival_link="survival_link.nsclc_os_growth_rate")
default.median_os, kg.median_os        # the metric choice re-ranks the resistance models
```

```bash
onkos simulate resistance.nsclc_first_line.two_population        # week-8 OS link (default)
# the k_g link is reachable through the Python survival_link= argument and the exports
```

The new link exports like any survival record (vt-JSON / JSON-LD / CSV / BibTeX) and carries
the universal clinical-use prohibition. (A `--survival-link` CLI flag is a trivial, deferred
add; the Python and export surfaces are the substantive ones.)

---

## 7. Source anchors (methodological; DOIs added at curation time)

- **The growth-rate constant as the dominant OS predictor.** Stein, W.D. et al. (2008),
  *Tumor growth rates derived from data for patients in a clinical trial correlate strongly
  with survival* (and the Stein/Bruno program) — the result that k_g out-discriminates early
  shrinkage for OS, the basis for the k_g link (already in `dataset/citations` as
  `stein-2008-grc`).
- **Model-based OS prediction from on-treatment metrics.** Claret, Bruno et al. — the
  framework that links a TGI metric to OS (the week-8 default), already cited.
- **Surrogate-endpoint validity.** The broad literature on early tumor-size endpoints as OS
  surrogates — the debate this feature makes computable rather than rhetorical.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not a claim that k_g is the right metric.** Onkos ships *both* and shows they disagree;
  it does not declare a winner. The contribution is making the choice legible.
- **Not a new metric.** k_g is already extracted; this spec only lets it drive survival. New
  bridge metrics (landmark size, tumor AUC, TTP, ctDNA) are clean follow-ons.
- **Not joint longitudinal-survival modeling.** The link remains a two-stage (metric →
  hazard) surrogate, not a joint model; the joint formulation is a separate, larger spec.
- **No individual-level survival.** Trial-level only.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only.** Every curve is a trial-level survival function of a
  published model under a declared metric, never a person's prognosis.
- **No therapy ranking.** The re-ranking is of *models under a metric*, a methodological
  statement; it never recommends a treatment.
- **The metric is declared, not fitted.** `link_metric`, `beta`, and the recorded C-index are
  illustrative, `unverified`, and tier-bounded; a more prognostic metric cannot raise a tier.
- **The line, restated.** Any feature that takes a real patient's tumor measurement and
  returns a survival estimate or a therapy choice **does not get built.** Making the
  surrogate-metric choice explicit and its consequences computable changes none of this.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Configurable metric** | `structure.link_metric` read by `simulate` (default week-8); nan ⇒ baseline covariate; backward-compatibility landmark. | Existing curves unchanged; a link can declare its driving metric. |
| **2 — The k_g link + the inversion** | a non-default growth-rate-constant OS link for NSCLC; the re-ranking / inversion landmarks. | The k_g link computes, and the resistance-model ranking inverts vs week-8. |
| **3 — Presentation** | a figure (week-8 vs k_g OS, with the inversion) + a CI-executed notebook; the higher recorded C-index. | The metric-choice axis is visualized and runs in CI. |
| **4 — Documentation** | README section framing the bridge metric as a model-selection axis; roadmap, layout, and the report counts updated. | The survival-metric axis is documented with the project's rigor and guardrails. |

Step 1 alone is a self-contained, honest contribution: making the on-treatment metric that
bridges tumor dynamics to survival a *declared, swappable* field — instead of a hidden
constant — is the smallest change that turns the field's surrogate-endpoint debate into
something a tool can compute. Step 2 shows it matters: the metric choice inverts the answer.
