# Workhorse-model verification dossier

**Status: literature-sourced, human PDF sign-off pending.** This document is the
evidence package for promoting the dataset's most-used records from
`review_status: unverified` toward `verified`. It was assembled by automated
literature review against accessible sources (open-access full text, PubMed
Central, regulatory reviews, and peer-reviewed reproductions of the original
tables). Per [`CONTRIBUTING.md`](../../CONTRIBUTING.md), **an LLM may assemble
this evidence but may not set `review_status: verified` on its own authority** —
a human must confirm each value against the source PDF and make the flip. Every
number below carries its source and an access status so that confirmation is a
lookup, not a re-derivation.

> **Why this matters.** The dataset ships every parameter value as *illustrative*
> by design (see the CHANGELOG preamble). That is honest, but it also means no
> value here can be cited in a real analysis yet. This dossier is the first step
> in closing that gap: it separates what is **structurally** correct (the model
> form matches the cited paper) from what is **numerically** grounded (the value
> matches a published estimate), and it flags the places where the current
> illustrative numbers are actively misleading and should be corrected before any
> promotion.

## How to read the confidence column

- **Structure: confirmed** — the model's equations match the cited paper, verified
  against an accessible source.
- **Values: confirmed (open)** — the published estimate was read from an
  open-access source (or an open reproduction of the original table). High enough
  to act on after a human cross-check.
- **Values: not located** — the original estimate sits behind a paywall and is not
  quoted in any accessible secondary source. Do **not** invent it; promotion of
  these records is blocked on institutional access to the source PDF.
- **Values: illustrative** — the current dataset number is a placeholder with no
  literature grounding (and in one case traces to a known-illustrative third-party
  default — see Claret below).

## Summary

| Record | Structure | Values | Status / action |
| --- | --- | --- | --- |
| `preclinical_translation.simeoni_2004.xenograft` | confirmed | confirmed (open reproduction) | ✅ **Applied** — real paclitaxel/A2780 set adopted; `review_status` stays `unverified` pending PDF (see §3) |
| `immuno_oncology.kuznetsov_1994.tumor_immune` | confirmed | **confirmed (open)** — already matches the published set | ✅ **Applied** — provenance upgraded; stays tier D / hypothesis (see §5) |
| `tgi_metrics.wang_2009.biexponential` | biexponential = Stein-Fojo form; Wang's own model was exp-decay + *linear* | not located | ✅ **Applied** — description corrected to fix the attribution (see §2) |
| `resistance.claret_2009.tgi` | confirmed | illustrative (CRC table not located; NSCLC variant available) | Re-anchor or relabel; C-index is illustrative-by-design, see note in §1 |
| `tgi_metrics.*.biexp` / `tgi_metrics.bruno_2020.*` | confirmed (biexponential is a real Stein-Fojo form) | illustrative | Ground magnitudes against Stein 2008 ranges (see §2) |
| `drug_effect.norton_simon.nsclc` | confirmed | illustrative (human Gompertz constants are a distribution) | Ground Vmax / growth constant against Norton 1988 (see §4) |

---

## §1 — Claret 2009 clinical TGI model

**Record:** `resistance.claret_2009.tgi`
**Source:** Claret L, et al. *J Clin Oncol* 2009;27(25):4103-4108. DOI
[10.1200/JCO.2008.21.0807](https://doi.org/10.1200/JCO.2008.21.0807), PMID 19636014.

**Structure — confirmed (high).** The model is a single ODE on tumor size:
`dTS/dt = KL·TS − KD(t)·Exposure·TS` with `KD(t) = KD0·exp(−λ·t)` — first-order
growth `KL`, exposure-driven kill `KD0`, and an exponentially decaying kill effect
(rate `λ`) encoding emergent resistance. This matches the record's kernel.
Confirmed via the MonolixSuite TGI library and BioModels MODEL1708310001 (both
accessible), which cite Claret 2009.

**Provenance of the paper's numbers — colorectal, not NSCLC.** Claret 2009 was
derived on **capecitabine in metastatic colorectal cancer (CRC)**, with tumor size
as RECIST sum of longest diameters and a **per-week** time base. The record is
parameterized for an illustrative **dacomitinib / NSCLC** context — i.e. it borrows
the *structure* but not the paper's *data*.

**Values — not located + one actively misleading number.**
- The CRC Table 1 estimates (KL, KD0, λ, baseline, and their IIV CV%) are **not
  quoted in any accessible source**; the JCO original and ResearchGate full text
  are paywalled. Promotion of the CRC numbers is blocked on institutional access.
- ⚠️ **Do not trust the circulating "Claret values" `KL=0.021, KD=0.025,
  λ=0.053`** (BioModels MODEL1708310001). That entry is non-curated, is in
  **per-day** units (implausible for the per-week CRC model — 0.021/day ≈ 0.15/week
  growth is far too fast for CRC), carries **no IIV** although the paper reports it,
  and uses a unitless baseline of 200. These are third-party illustrative defaults.
  **The record's own `kL=0.021` / `λ=0.061` appear to descend from this same
  illustrative source** and must not be presented as Claret's estimates.
- ⚠️ **The record's `predictive_performance` C-index is the wrong metric type for
  this paper.** It lists `OS_C_index_external = 0.62`. Claret 2009 reported **no
  C-index**; it reported a **median-OS prediction**: predicted median OS **431 days
  (90% PI 362–514)** vs observed **401 days** on the phase-III capecitabine arm, and
  a predicted survival improvement of **39 days (90% PI −21 to 110)** vs observed
  35 days (parametric lognormal survival model; baseline tumor size and week-7
  tumor-size change were significant at P < .00001). Source: PMID 19636014 abstract
  (accessible).

  **Important caveat — do not delete this field in isolation.** An external
  C-index labeled *(illustrative)* is carried on **37 records** as the dataset's
  uniform validation currency, and it is **load-bearing**: `onkos.audit` uses its
  presence as the external-validation gate for tier A/B, `onkos.combine` reads it
  for evidence-weighted model averaging, and the discriminability/survival-metric
  views compare the values across links. It is illustrative *by the same design as
  every parameter value*, not a one-off fabrication. The correct fix is therefore
  **not** to strip Claret's C-index (that would be inconsistent and would break
  tested audit/averaging behavior) but, when this record is promoted, to record the
  metric Claret actually reported (median-OS prediction) and decide dataset-wide
  whether the illustrative C-index scheme should be replaced with real, per-record
  validation metrics. That is a maintainer-level design decision, flagged here.

**A real NSCLC grounding does exist (biexponential Claret-family).** An open-access
external-validation study in **ALK+ NSCLC** (alectinib vs crizotinib;
[PMC10363035](https://pmc.ncbi.nlm.nih.gov/articles/PMC10363035/)) fits a
biexponential TGI model "implemented as in Claret et al." (growth `KG` + shrinkage
`KS`, **no λ term**), per-week:

| Arm | KG (1/wk) | KS (1/wk) | TS0 (mm) |
| --- | --- | --- | --- |
| Alectinib | 0.00196 (RSE 13.2%, IIV CV 138%) | 0.0342 (RSE 9.8%, IIV CV 97.1%) | 57.0 (IIV CV 66.1%) |
| Crizotinib | 0.00438 (RSE 12.1%, IIV CV 105%) | 0.0373 (RSE 10.5%, IIV CV 91.3%) | 57.0 |

**Recommended actions (maintainer):**
1. Immediately correct or remove the fabricated `OS_C_index_external = 0.62`; if a
   validation field is wanted, record the real median-OS prediction (431 vs 401 d)
   and cite it, or leave the field empty until the paywalled value is confirmed.
2. Decide the record's intent: either **(a)** re-anchor it to the actual CRC
   capecitabine context (and obtain the Table 1 values from the PDF), or **(b)**
   keep the NSCLC context but relabel it honestly as a biexponential Claret-family
   record grounded on PMC10363035 (note: that variant has no λ resistance term, so
   the resistance axis would need a separate, explicitly-illustrative record).
3. The IIV CV magnitudes the record carries (kD 89%, λ 96%) are *qualitatively*
   well-supported — the NSCLC fits above show 90–140% CV on the kill/growth terms —
   so the "resistance terms carry ~90% CV" thesis is real even though these specific
   numbers are illustrative.

---

## §2 — Stein-Fojo growth-rate constants & Wang 2009 NSCLC

**Records:** `tgi_metrics.wang_2009.biexponential`, the `*.biexp` family, and
`tgi_metrics.bruno_2020.breast_biexponential`.

**Stein 2008 — confirmed (open), with a context correction.**
[Stein WD, et al. *The Oncologist* 2008;13(10):1046-1054](https://pmc.ncbi.nlm.nih.gov/articles/PMC3313464/)
(open access). The growth-rate-constant model is `f(t) = exp(−d·t) + exp(g·t) − 1`
(`f` = marker normalized to day 0, `t` in **days**, `g` growth-rate constant, `d`
regression-rate constant). ⚠️ **Context correction:** the original cohort is
**metastatic castration-resistant prostate cancer measured by PSA** (112 patients,
two NCI trials; validation 42 patients on ixabepilone), *not* renal cell carcinoma
/ sorafenib (those are companion papers). Published magnitudes:
- `g`: median **10^−2.5 day⁻¹ ≈ 0.0032/day (≈ 0.022/week)**, ~1,500-fold spread.
- `d`: median **10^−1.7 day⁻¹ ≈ 0.020/day (≈ 0.14/week)**, ~50-fold spread.
- Prognostic strength (this is the headline the dataset leans on): survival
  correlated with **log(g)** at Pearson **r = −0.72 (p < 0.0001)** and only weakly
  with log(d) (r = −0.22); HR for above- vs below-median log(g) was **5.14 (95% CI
  3.10–8.52)**. This strongly supports the dataset's "growth rate is the
  well-identified, prognostic term" tier rationale.

**Wang 2009 — structural mismatch + values not located.**
[Wang Y, et al. *Clin Pharmacol Ther* 2009;86(2):167-174](https://pubmed.ncbi.nlm.nih.gov/19440187/)
(abstract accessible; full text paywalled). ⚠️ **The record models a *biexponential*
form, but Wang 2009 used "a mixed exponential decay and linear growth model"**:
`TS(t) = BSL·exp(−SR·t) + PR·t`, where `SR` is an exponential shrinkage-rate
constant and **`PR` is a *linear* progression rate (mm/time), not an exponential
rate constant**. The OS model is parametric with three predictors: **ECOG, baseline
tumor size, and week-8 tumor-size change** (advanced NSCLC, four FDA registration
trials). The numeric SR/PR estimates and the week-8 hazard ratio are **not quoted
in any accessible source** — blocked on the paywall.

**Recommended actions (maintainer):**
1. ✅ **Applied (partial).** The record's description now states explicitly that the
   biexponential form is the **Stein-Fojo growth-rate-constant model** and that Wang
   2009's own NSCLC model was a distinct exp-decay-plus-linear parameterization; the
   record id and citation are unchanged to avoid breaking the many references to it.
   A *fuller* fix — renaming the record id to drop the Wang attribution, or adding a
   real exp-decay-plus-linear kernel that matches Wang 2009 — is left to the
   maintainer because the id is referenced across the README, tests, and notebooks.
2. ✅ **Applied** to `tgi_metrics.wang_2009.biexponential`: its `kg` (0.018/week) and
   `ks` (0.055/week) `source_locator`s now cite the open-access Stein 2008 ranges
   (g median ≈ 0.022/week, ~1500-fold spread; d median ≈ 0.14/week, ~50-fold spread)
   and confirm both values sit within them. This is a deliberately *weaker* grounding
   than Simeoni's: it is a **cross-tumor magnitude anchor** (Stein's cohort was
   prostate/PSA), not a same-context fit, and the values were left unchanged because
   they are already in-range and the record is referenced across 9 test files. The
   other in-context `*.biexp` records still carry illustrative magnitudes.
3. `bruno_2020.breast_biexponential` cites a *review* (Bruno 2020 CCR) for specific
   values; reviews summarize rather than estimate, so these should stay
   illustrative or be re-anchored to a primary breast-cancer TGI fit.

---

## §3 — Simeoni 2004 preclinical xenograft model

**Record:** `preclinical_translation.simeoni_2004.xenograft`
**Source:** Simeoni M, et al. *Cancer Research* 2004;64(3):1094-1101. DOI
[10.1158/0008-5472.CAN-03-2524](https://doi.org/10.1158/0008-5472.CAN-03-2524),
PMID 14871843.

**Structure — confirmed (high).** Exponential-then-linear unperturbed growth
`dx1/dt = λ0·x1 / [1 + (λ0/λ1·w)^ψ]^(1/ψ)`; drug damages proliferating cells at
rate `k2·c(t)·x1`; damaged cells traverse a transit chain `x2→x3→x4` at rate `k1`
before dying; observed weight `w = x1+x2+x3+x4`; `ψ` fixed at 20. Confirmed via the
PubMed abstract, BMC Cancer [PMC7076937](https://pmc.ncbi.nlm.nih.gov/articles/PMC7076937/),
and the MonolixSuite TGI library.

**Values — confirmed (open reproduction).** The original paper is paywalled (AACR
403), but its parameter table is reproduced in
[Strömberg & Hooker, PMC5660732](https://pmc.ncbi.nlm.nih.gov/articles/PMC5660732/)
(open access), Table I, attributed to Simeoni 2004. Units: λ0 (1/day), λ1 (g/day),
w0 (g), k1 (1/day), k2 (ng⁻¹·mL·day⁻¹), ψ fixed = 20.

| Compound (cell line) | λ0 (1/day) | λ1 (g/day) | w0 (g) | k1 (1/day) | k2 (ng⁻¹·mL·day⁻¹) |
| --- | --- | --- | --- | --- | --- |
| Paclitaxel (A2780) exp.3 | 0.311 | 0.656 | 0.033 | 0.968 | 6.29e-4 |
| Paclitaxel (A2780) exp.4 | 0.273 | 0.814 | 0.055 | 0.968 | 6.29e-4 |
| 5-FU | 0.215 | 0.412 | 0.065 | 0.056 | 2.021e-3 |
| Compound A | 0.349 | 0.363 | 0.010 | 0.405 | 3.45e-4 |
| Compound B exp.7 | 0.309 | 0.796 | 0.034 | 0.517 | 2.89e-4 |
| Compound B exp.8 | 0.369 | 0.511 | 0.016 | 0.615 | 2.93e-4 |

Untreated A2780 growth is characterized by λ0 ≈ 0.27–0.31/day, λ1 ≈ 0.66–0.81 g/day.
The threshold ("tumor-static") concentration concept `Ct ≈ λ0/k2` is part of the
framework (a worked value of ~1100 ng/mL is quoted in the literature); treat the
closed form as standard-but-unverified-at-source.

**✅ Applied.** The record now carries the **paclitaxel/A2780 exp.3** values
(λ0=0.311/day, λ1=0.656 g/day, k1=0.968/day, k2=6.29e-4 mL·ng⁻¹·day⁻¹, ψ=20 fixed),
with each parameter's `source_locator` citing PMC5660732 Table I and this dossier.
The implied tumor-static concentration is now `Ct = λ0/k2 ≈ 494 ng/mL`; the
qualitative tests in `tests/test_preclinical.py` were re-dosed to physical
concentrations around/above that threshold (the old illustrative doses were matched
to a potency ~19× too high) and the dataset-health report was regenerated.
**`review_status` deliberately remains `unverified`:** the source is a high-quality
*secondary reproduction*, so the `verified` flip still requires a human reading
Simeoni 2004 Table I directly — the values are high-confidence, the sign-off is not
an LLM's to give.

---

## §4 — Norton-Simon kill on Gompertz growth

**Record:** `drug_effect.norton_simon.nsclc`
**Sources:** Norton L. *The Oncologist* 2005;10(6):370-381 (open); Norton L. "A
Gompertzian Model of Human Breast Cancer Growth." *Cancer Res* 1988;48:7067-7071.

**Structure — confirmed (high).** Kill proportional to the unperturbed Gompertzian
growth rate: `dV/dt = g·V·ln(Vmax/V) − k·E·V·ln(Vmax/V)`. Confirmed via
[Mistry et al., PMC3228251](https://pmc.ncbi.nlm.nih.gov/articles/PMC3228251/)
(open), which writes the treated form `dN/dt = {αN − βN·ln N}·[1 − D(t)]` — exactly
the kill-proportional-to-growth-rate structure.

**Values — illustrative; the human growth constant is a distribution, not a point.**
Norton 1988 (read from an open PDF mirror) reports, for **human breast cancer**:
carrying capacity **N(∞) = 3.1 × 10¹² cells**, lethal burden **10¹² cells**, and a
growth-rate constant **b that is lognormal** with mean `log_e(b) = −2.9`, SD 0.71,
with `t` in **months** (median b ≈ 0.055/month). The record's `Vmax = 200 mm` and
`g = 0.08/week` are illustrative NSCLC-scale numbers, not these. For experimental
solid tumors, [Vaghi et al., PLOS Comput Biol 2020, PMC7059968](https://pmc.ncbi.nlm.nih.gov/articles/PMC7059968/)
(open) gives Gompertz `dV/dt = αV·ln(K/V)` with breast α = 0.58/day, K = 2,600 mm³
(mouse) — note these are murine, not human.

**Recommended action (maintainer):** the Norton-Simon record is best kept explicitly
illustrative — there is no single published human NSCLC Gompertz constant to anchor
it, and Norton's own message is *kinetic heterogeneity* (a distribution). Document
the growth-constant IIV honestly rather than implying a point estimate. The
record's purpose (contrast the kill *mechanism* against log-kill in the divergence
view) does not require a verified value.

---

## §5 — Kuznetsov 1994 tumor-immune model

**Record:** `immuno_oncology.kuznetsov_1994.tumor_immune`
**Source:** Kuznetsov VA, Makalkin IA, Taylor MA, Perelson AS. *Bull Math Biol*
1994;56(2):295-321, PMID 8186756.

**Structure & values — confirmed (open).** This is the one workhorse record whose
numbers are **already the genuine published values.** The original is paywalled, but
the canonical **nondimensional** parameter set is reproduced verbatim in
[d'Onofrio, arXiv 1309.3337](https://arxiv.org/abs/1309.3337) (peer-reviewed in
*Math Comput Modelling* 2008) and independently corroborated by
[Marchant et al., arXiv 1909.05203](https://arxiv.org/abs/1909.05203) (*J Theor Biol*
2020) and EBI BioModels BIOMD0000000762. The dataset's values match:

| Symbol | Record | Published | Match |
| --- | --- | --- | --- |
| σ (s) source | 0.118 | 0.1181 | ✓ |
| ρ (rho) recruitment | 1.131 | 1.131 | ✓ |
| η (eta) half-saturation | 20.19 | 20.19 | ✓ |
| μ (mu) inactivation | 0.00311 | 0.00311 | ✓ |
| δ (delta) effector death | 0.374 | 0.3743 | ✓ |
| α (alpha) tumor growth | 1.636 | 1.636 | ✓ |
| β (beta) inverse capacity | 0.002 | 0.002 | ✓ |

Scaling (d'Onofrio): `t_true = 9.9·t_nondimensional` days; one nondimensional cell
unit = 10⁶ cells.

**Action taken:** the `extraction.source_locator` for these seven dynamics
parameters has been upgraded from "illustrative" to cite the confirming open
reproductions — **no value changed.** The record correctly **stays tier D /
hypothesis / non-predictive**: the parameters reproduce immune control, dormancy,
and escape *qualitatively* (they were fit to BCL₁ B-lymphoma in mouse spleen, and
the immunotherapy effect `E` is a hypothesis-tier augmentation). Verified parameter
*values* do not upgrade a model that is non-predictive *by construction* — exactly
the distinction the tier system exists to make.

---

## Maintainer checklist to set `review_status: verified`

For each record, per [`CONTRIBUTING.md`](../../CONTRIBUTING.md):

1. Open the source PDF (this dossier gives DOI/PMID and the specific table).
2. Confirm structure, every value + units + IIV, the derivation context, and the
   transportability boundary.
3. Where this dossier says "not located," obtain the value from the PDF — do not
   import third-party defaults (see the Claret BioModels warning).
4. Flip `review_status` and update `source_locator` to the confirmed table/page.
5. Re-run `onkos validate`, `onkos report` (and `diff` against
   `docs/dataset-health.md`), `pytest -q`, and the affected notebooks.
