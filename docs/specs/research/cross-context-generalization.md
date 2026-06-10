# Onkos — research spec: cross-context generalization (the findings are not NSCLC artifacts)

**Status:** implemented in v0.29.0 (mechanistic resistance models + tail-sensitive survival
links for breast, CRC, HCC, melanoma). This is a *breadth* release: it adds no new
methodology — it tests whether the project's headline findings, all previously demonstrated
on NSCLC first line, **generalize across solid-tumor contexts**. The methodological source
anchors are those of the underlying specs (v0.24–v0.28). Written in the v0.1 house style;
every value is illustrative and `unverified` by design — the infrastructure and the
generalization are the contribution.

**The credibility question.** By v0.28 the project had eight analytical capabilities — the
resistance-mechanism divergence axis, the model-selection budget, the ORR → OS surrogate,
and the depth-vs-durability finding among them — but each headline *result* lived in a
single context: NSCLC first line. A reviewer's first question is therefore unavoidable: *do
these findings generalize, or are they artifacts of one tumor type's illustrative
parameters?* This release answers it by giving the four other curated solid-tumor contexts
(breast, CRC, HCC, melanoma) the **same two pieces** NSCLC had — a mechanistic
two-population resistance model and a tail-sensitive (growth-rate-constant) OS link — and
showing that every NSCLC-only result reappears, unchanged in direction, in all four.

> The corrected recipe (from the v0.28 post-mortem): the ORR-surrogate inversion and the
> budget's survival-link axis need *a model that decouples response depth from durability* in
> the context, not merely a second survival link. A 2-model context where the resistance
> model wins on both early shrinkage and tail growth is concordant by construction. Adding
> the two-population model — universal response, but a fast resistant regrowth — supplies the
> decoupling, and the k_g link supplies the tail-sensitive endpoint that sees it.

---

## 1. What was NSCLC-only, and why

| Finding | Spec | Why it was NSCLC-only |
| --- | --- | --- |
| Resistance *mechanism* as a divergence axis (phenomenological vs mechanistic) | v0.24 | only NSCLC had a two-population model |
| Budget **survival-link axis** `V_link` (and the model×link interaction) | v0.26 | only NSCLC had ≥2 OS survival links |
| ORR → OS surrogate **conditional on the survival mechanism** | v0.27 | needs the tail-sensitive `k_g` link *and* a depth≠durability model |
| **Depth ≠ durability** (highest ORR = shortest DoR) | v0.28 | needs the two-population model's broad-but-brief responses |

All four share a root cause: the non-NSCLC contexts carried only a phenomenological (Claret)
resistance model + a biexponential model + a single week-8 OS link. That set is *concordant
by construction* — the Claret model has both the deepest early shrinkage (highest ORR) and
the slowest tail growth (best k_g-OS), so ORR and tail-driven OS agree, the budget's
survival-link axis is empty, and there is no resistance-mechanism contrast.

---

## 2. The two additions per context

For each of breast, CRC, HCC, and melanoma first line, v0.29 adds:

1. **A mechanistic two-population resistance model** (`resistance.<tt>_first_line.two_population`,
   the `two_population_resistance` kernel): a sensitive clone killed by the drug + a
   pre-existing resistant clone that outgrows. The kill potency `kd` is matched to that
   context's Claret model (so the contrast isolates the resistance *mechanism*), and the
   resistant growth `kgr` is set *faster* than the sensitive growth — so the model responds
   in nearly everyone (deep early shrink, high ORR) but its responses are **brief** (a fast
   resistant regrowth → short DoR → a high tail growth-rate constant).

2. **A tail-sensitive `k_g` OS link** (`survival_link.<tt>_os_growth_rate`, non-default,
   `link_metric: tumor_growth_rate_kg`, calibrated per context): the survival endpoint that
   reads the regrowth tail rather than the early shrinkage — the dimension on which the
   broad-but-brief responder is penalized.

Both are tier C, illustrative, and `unverified`; citations are the existing `foo-michor-2014`
/ `goldie-coldman-1979` (resistance) and `stein-2008-grc` (growth-rate survival).

---

## 3. The result — every finding generalizes

With the two additions, each context goes from a 2-model / 1-link layout to a 3-model /
2-link one, and every NSCLC-only result reappears (illustrative, 0–312 wk horizon):

| Context | ORR→OS under week-8 | ORR→OS under k_g | budget survival-link share | dominant axis |
| --- | --- | --- | --- | --- |
| **NSCLC** 1L | concordant (0/6) | **discordant** (4/6) | 24% | model×link interaction |
| **breast** 1L | concordant (0/3) | **discordant** (2/3) | 72% | survival-link |
| **CRC** 1L | concordant (0/3) | **discordant** (2/3) | 52% | survival-link |
| **HCC** 1L | concordant (0/3) | **discordant** (3/3) | 21% | parameter |
| **melanoma** 1L | concordant (0/3) | **discordant** (2/3) | 54% | survival-link |

In **every** context: (a) the two-population model has the **highest ORR** but the
**shortest DoR** (depth ≠ durability); (b) under the tail-sensitive `k_g` link that same
model has the **worst** OS, so ORR — faithful under the week-8 surrogate — **mis-ranks** OS
(the conditional surrogacy of v0.27); and (c) the budget's survival-link axis is now real,
making the survival-model choice the **dominant** structural axis in three of the four new
contexts. The model-selection budget summary now flags **5 of 6 contexts as
structure-dominated** (was 4/6, with the survival-link axis empty for the non-NSCLC ones).

The four headline findings are therefore **not NSCLC artifacts** — they reproduce, unchanged
in direction, across five solid-tumor contexts. That generalization is the contribution.

---

## 4. Tier & guardrail propagation (unchanged invariants)

- **Worst-input-wins, transport floors, population-only** — every new record obeys the same
  invariants as the rest of the dataset (tier C, out-of-context → D, no individual output).
- **The k_g links are non-default and opt-in**, exactly like the NSCLC one, so they never
  perturb the default divergence view; they enter only through `survival_link=` and the
  budget's eligible-link enumeration.
- **Calibration is declared, not fitted.** The per-context `kgr` and the k_g-link
  `scale`/`beta` are illustrative values chosen so the contexts are on a realistic OS scale;
  they are `unverified` and tier-bounded, and no finding depends on their precise values —
  only on the structural fact that the resistant clone outgrows.

---

## 5. Validation landmarks

No new kernel or analysis — the validation is that the existing landmark-tested machinery
produces the generalized result. Two cross-context tests pin it:

| Landmark | Condition |
| --- | --- |
| **Survival-link axis populated** | every first-line solid-tumor context has ≥2 OS links and `v_link > 0` (`tests/test_budget.py`). |
| **ORR-surrogate generalizes** | in every first-line solid-tumor context, ORR is concordant with OS under the week-8 link and **discordant** under the k_g link (`tests/test_response.py`). |
| **Single-link control** | NSCLC second line still has one OS link, so the v0.21 collapse + `single_survival_link` flag remain tested. |
| **Round-trip** | the new two-population records are multi-state ODE kernels; they round-trip to SBML/NONMEM through the existing `two_population_resistance` binding (covered by the kernel's round-trip, the records share it). |

---

## 6. API, CLI, and surface

No new surface — the additions are dataset records, reached through every existing tool:

```python
# The resistance-mechanism divergence, ORR-surrogate, and budget now work in any context.
ctx = dict(tumor_type="breast", line="first")
onkos.compare(ds, purpose="tgi", context=ctx)                 # 3 models incl. two-population
onkos.response_vs_survival(ds, context=ctx,
                           survival_link="survival_link.breast_os_growth_rate")  # discordant
onkos.model_selection_budget(ds, context=ctx, endpoint="OS")  # survival-link axis populated
```

```bash
onkos response --surrogate --tumor-type breast --survival-link survival_link.breast_os_growth_rate
onkos budget --tumor-type CRC --line first
```

---

## 7. Source anchors

The methodological anchors are those of the underlying specs: Goldie-Coldman 1979 /
Foo-Michor 2014 (two-population resistance, v0.24), Stein 2008 (the growth-rate-constant OS
link, v0.25), RECIST 1.1 / Eisenhauer 2009 (response, v0.27–v0.28). This release adds curated
*content*, not methods, so it introduces no new methodological citations.

---

## 8. Deliberate non-goals (so the scope stays honest)

- **Not real per-context parameters.** The new models are illustrative and tier C; this
  release demonstrates *generalizability of the findings*, not validated tumor-type
  parameters. Promoting them is the standard PDF-verification work.
- **Not second-line breadth.** First-line solid-tumor contexts only; second line and the
  preclinical/IO subsystems are unchanged.
- **Not a Cox link per context.** Only the tail-sensitive k_g link is added (the one the
  findings need); a per-context Cox link (a second survival-structure axis) is a clean,
  separable extension.
- **No new endpoint or analysis.** Pure curated content over the v0.24–v0.28 machinery.

---

## 9. Safety & scope (unchanged hard line)

- **Population / trial level only**, no individual output, no therapy ranking — the new
  records obey every guardrail of the subsystems they join.
- **Illustrative and unverified by design**, tier-bounded; the generalization is of the
  *structural findings*, not a claim about any real drug or tumor type.
- **The line, restated.** Adding curated content that makes the headline findings
  reproducible across contexts changes nothing about the population-only, no-prognosis,
  no-recommendation boundary.

---

## 10. Phased steps

| Step | Content | Done = |
| --- | --- | --- |
| **1 — Records** | a two-population model + a k_g OS link for breast, CRC, HCC, melanoma (8 records), calibrated per context. | The dataset validates; each context has a 3-model / 2-link layout. |
| **2 — Generalization tests** | the survival-link-axis-populated and ORR-surrogate-generalizes cross-context landmarks; the single-link control moved to NSCLC 2L. | The findings reproduce in all five contexts and CI enforces it. |
| **3 — Presentation** | a cross-context figure (the week-8 vs k_g discordance and the structure-dominance across tumor types) + a CI-executed notebook; refreshed report/figure numbers. | The generalization is visualized and the report shows 5/6 structure-dominated. |
| **4 — Documentation** | README updated so the headline findings are stated as dataset-wide, with the per-context table; roadmap. | The project no longer reads as an NSCLC demo. |

Step 1–2 alone is the contribution: showing, with CI-enforced tests, that the resistance-
mechanism divergence, the model-selection budget's survival-link axis, the conditional
ORR → OS surrogacy, and the depth-vs-durability finding all reproduce across five solid-tumor
contexts — turning four single-context demonstrations into a general claim about how these
oncology models behave.
