# Onkos

**A curated, citation-backed, tier-annotated dataset of tumor-growth-inhibition
(TGI) models, exposure-response links, and TGI-metric → survival models — the
machinery oncology drug development runs on — exported into the standard
pharmacometric and systems-biology formats (NONMEM, SBML, PharmML,
nlmixr2/rxode2, Pumas).**

> ⚠️ **NOT a clinical decision tool. NOT a prognostic calculator. NOT a
> treatment recommender.** Population/trial-level forward simulation only, for
> drug-development methodology, simulation, and education. Every export carries
> `onkos:clinicalUse = "PROHIBITED — research / drug-development / education only"`.

*Onkos* (Greek *ὄγκος*, "mass, swelling") is the literal root of *onco-*. It is
the third in a family with **Nidus** (gestational physiology) and **Hypnos**
(anesthetic PK/PD), sharing one thesis: **a model is only as trustworthy as its
weakest, least-validated input — so make that a first-class, machine-readable
field.**

[![CI](https://github.com/clay-good/onkos/actions/workflows/ci.yml/badge.svg)](https://github.com/clay-good/onkos/actions/workflows/ci.yml)
&nbsp;v0.10 · Code: MIT · Data: CC-BY-4.0 · Python ≥ 3.9

---

## The problem

Oncology has the highest drug attrition of any therapeutic area, and the field's
response is model-informed drug development: link drug exposure to tumor-size
dynamics, and tumor dynamics to overall survival (OS), so early data forecast
late outcomes and gate go/no-go decisions. The workhorse models (Gompertz,
Simeoni, **Claret**, Stein/Bruno growth-rate-constant) live in per-drug,
per-trial papers, carry enormous and under-communicated uncertainty (resistance
terms with ~90% CV are routine), and are **derived in one context then silently
transported to another**, where their predictive validity is unknown.

Onkos is the missing curated layer: it says, honestly, *which TGI model and which
parameters, for which tumor type and line, derived from which trial, validated
how far beyond it, with what confidence — and how much the survival prediction
changes if you'd picked a different model.*

---

## The headline feature: virtual-trial divergence

Pick a tumor type, line, and drug-effect size. Onkos overlays the simulated
tumor-size and **population OS** curves across *every eligible TGI model*, greys
out the models whose `transportability` envelope the context violates (with the
reason), and quantifies the divergence in the survival prediction. **This makes
model-selection risk in go/no-go decisions measurable** — the exact risk that,
unquantified, sends drugs into doomed phase-3 trials.

![Virtual-trial divergence](docs/images/divergence.png)

In the figure above (NSCLC, first line, E = 1.0), two NSCLC-validated models that
fit early tumor data comparably imply median OS anywhere from ~54 to ~91 weeks.
Every model validated only on another tumor type is **greyed out automatically**
because applying it to NSCLC leaves its validated envelope (tier → D + warning).
That spread *is* the model-selection risk.

```text
$ onkos simulate --compare --tumor-type NSCLC --line first --drug-effect 1.0

  [C] resistance.claret_2009.tgi              median OS 90.8
  [C] tgi_metrics.wang_2009.biexponential     median OS 53.7
  [-] resistance.crc_first_line.claret        EXCLUDED
        (tumor_type 'NSCLC' is outside validated ['CRC'] -> tier_down_to_D and warn)
  [-] ... 6 more excluded for out-of-context transport (breast, HCC, melanoma)

  OS divergence (max pointwise): 0.247
  Median OS range: (53.7, 90.8)
```

### The second uncertainty axis: parameter variability

The divergence view quantifies *model-selection* uncertainty. The other axis is
*parameter* uncertainty: the dataset records inter-individual variability
(`iiv_cv_percent`) on its high-uncertainty kill/resistance terms specifically so
they cannot pose as point estimates — and `onkos.simulate_ensemble` makes that
stored variability flow into the prediction. Parameters with an IIV CV are
sampled lognormally (the standard pharmacometric convention; the median is
preserved) and the tumor-size, TGI-metric, and population-OS distributions are
returned as bands.

![Parameter-uncertainty bands](docs/images/uncertainty.png)

For the Claret NSCLC model, whose resistance and kill terms carry ~90% CV, this
turns a deceptively precise "80% week-8 shrinkage" into an honest
`[-100%, -30%]` band and a median-OS interval of roughly 58–103 weeks — the
uncertainty was always in the data; now it is in the answer.

### TGI metrics — the Stein/Bruno panel

Every simulated trajectory is summarized into the derived metrics oncology
pharmacometrics actually reports (spec §3, §6): depth of response, nadir and
**time-to-growth**, the **tumor growth-rate constant k_g** (the strongly
prognostic Stein/Bruno quantity), the shrinkage-rate constant k_s, and the
RECIST-style **duration of response** (partial response → progression). They feed
both the survival link (via the week-8 change) and the reports.

![TGI-metric extraction](docs/images/tgi_metrics.png)

The extractor is **model-agnostic** — it estimates k_g / k_s from the trajectory
the way the Stein method estimates them from RECIST data, so the metrics are
comparable across the Claret, biexponential, and Simeoni kernels. It is also
self-checking: run on the biexponential kernel it recovers that kernel's
*generating* k_g and k_s to within ~10%, and on the Claret model the extracted
k_g recovers the model's growth constant k_L. Metrics that don't apply (no
regrowth, no RECIST response) are returned as `nan`, never fabricated.

```python
m = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential", context=ctx).metrics
m["tumor_growth_rate_kg"]       # late-phase log-linear regrowth rate (≈ generating kg)
m["tumor_shrinkage_rate_ks"]    # initial shrink rate via k_s = k_g − s0
m["time_to_growth_weeks"]       # nadir time when genuine regrowth follows
m["duration_of_response_weeks"] # RECIST PR (−30%) → PD (+20% from nadir); nan if no PR
```

### Which uncertainty to verify first: sensitivity analysis

Propagation gives a band; `onkos.sensitivity` attributes that band's variance to
individual parameters. Because each IIV-bearing parameter is sampled
independently, a parameter's standardized regression coefficient equals its
correlation with the target and the squared coefficients partition the explained
variance (a first-order Sobol decomposition). This is **curation triage**: the
spec's highest-leverage contribution is verifying records against the source PDF
(§9), and this says *which* parameter's uncertainty actually moves the survival
prediction.

![Parameter sensitivity tornado](docs/images/sensitivity.png)

A genuinely useful nuance falls out: for the Claret NSCLC model the kill rate
`kD` (CV 89%) drives ~90% of the median-OS variance — *more* than the resistance
term `lambda` (CV 96%), because influence is variability × effect-strength, not
CV alone. So `kD` is where verification pays off first.

```python
res = onkos.sensitivity(ds, "resistance.claret_2009.tgi", context=ctx, target="median_os_weeks")
res.r_squared                  # first-order variance explained
res.dominant.symbol            # "kD" — verify this parameter first
[(p.symbol, round(p.contribution, 2), p.src) for p in res.indices]
```

---

## Install & quick start

```bash
git clone https://github.com/clay-good/onkos
cd onkos
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"        # or: pip install -e .   (runtime only)

onkos validate                 # JSON-Schema-validate the dataset
onkos info                     # counts by subsystem / tier / review status
onkos simulate --compare       # the divergence view (NSCLC, 1L by default)
```

### Python API (cheat sheet)

```python
import numpy as np
import onkos

ds = onkos.load()
m = ds["resistance.claret_2009.tgi"]
m.tier                                     # "C"
m.derivation_context.tumor_type            # "NSCLC"
m.transportability.validated_tumor_types   # ("NSCLC",)
m["lambda"].iiv_cv_percent                 # 96   -> uncertainty is first-class
m.review_status                            # "unverified"
m.primary_citation.doi                     # "10.1200/JCO.2008.21.0807"

# Population-level forward simulation (NO individual prognosis, NO therapy ranking)
ctx = dict(tumor_type="NSCLC", line="first")
traj = onkos.simulate(ds, "resistance.claret_2009.tgi",
                      context=ctx, drug_effect=1.0, t=np.linspace(0, 104, 209))
traj.tumor_size, traj.os_curve             # tumor-size + population OS trajectory
traj.tier, traj.warnings                   # propagated tier + transport warnings
traj.metrics["week8_relative_change"]      # the TGI metric feeding the survival link

# Virtual-trial comparison — the headline feature
cmp = onkos.compare(ds, purpose="tgi", context=ctx, drug_effect=1.0)
cmp.os_divergence                          # model-choice dependence of the OS prediction
cmp.median_os_range                        # (lo, hi) median OS across models
cmp.excluded                               # models greyed out for out-of-context transport

# Parameter uncertainty — propagate the stored IIV CVs (Monte-Carlo bands)
ens = onkos.simulate_ensemble(ds, "resistance.claret_2009.tgi", context=ctx, n=400, seed=0)
ens.tumor_size.median, ens.tumor_size.lo, ens.tumor_size.hi   # 5–95% band arrays
ens.os_curve.lo, ens.os_curve.hi                               # population-OS band
ens.metrics["median_os_weeks"]             # {"median", "lo", "hi"}

# Sensitivity — attribute the OS-prediction variance to parameters (verify first)
res = onkos.sensitivity(ds, "resistance.claret_2009.tgi", context=ctx, target="median_os_weeks")
res.dominant.symbol                        # "kD"
res.indices                                # ranked [ParamSensitivity(symbol, src, contribution), …]
```

### CLI (cheat sheet)

| Command | Does |
| --- | --- |
| `onkos version` | print version |
| `onkos validate` | JSON-Schema + referential-integrity check of the dataset |
| `onkos info` | counts by subsystem / tier / review status |
| `onkos report [--output FILE]` | dataset health & external-validation report (Markdown) |
| `onkos simulate <id> [--tumor-type --line --drug-effect]` | one model's trajectory + metrics |
| `onkos simulate --compare` | virtual-trial divergence across eligible models |
| `onkos uncertainty <id> [--n --seed]` | Monte-Carlo parameter-uncertainty bands (propagates IIV CV) |
| `onkos sensitivity <id> [--target --n]` | rank parameters by how much their IIV drives a target metric |
| `onkos export --format <fmt> --output <dir>` | generate artifacts |

Export formats: `nonmem`, `sbml`, `pharmml`, `so` (PharmML Standard Output),
`rxode2`, `pumas`, `vt-json`, `omex`, `csv`, `bibtex`. The COMBINE `.omex` bundles
SBML + PharmML + the SO + virtual-trial JSON + provenance into one citable
archive.

### Dashboard

```bash
pip install -e ".[dashboard]"
streamlit run dashboard/app.py
```

---

## The record — the unit of curation

A record is a structured object, not a scalar. Two kinds share one schema: a
**model** record (e.g. the Claret 2009 TGI model) and a **context-baseline**
record (e.g. NSCLC first-line baseline growth). The fields that carry the
project:

- **`derivation_context`** — the exact drug, tumor type, line, trial, and
  measurement basis a parameter came from. Machine-readable, mandatory.
- **`transportability`** — how far beyond that origin it has actually been
  validated. Crossing this boundary forces a tier penalty.
- **`iiv_cv_percent`** — inter-individual variability on the high-uncertainty
  kill/resistance terms, so a 90%-CV term cannot masquerade as a point estimate.

```jsonc
{
  "id": "resistance.claret_2009.tgi",
  "kind": "model", "purpose": "tgi", "subsystem": "resistance",
  "kernel": "claret_tgi",
  "structure": { "growth_law": "exponential",
                 "kill_model": "first_order_exposure_driven",
                 "resistance": "exponential_decay_of_kill" },
  "parameters": [
    { "symbol": "kL",     "tier": "B", "value": {"central": 0.021, "units": "1/week"} },
    { "symbol": "kD",     "tier": "C", "iiv_cv_percent": 89, "value": {"central": 0.30, "units": "1/week per effect-unit"} },
    { "symbol": "lambda", "tier": "C", "iiv_cv_percent": 96, "value": {"central": 0.061, "units": "1/week"} }
  ],
  "derivation_context": { "drug": "dacomitinib", "drug_class": "EGFR_TKI",
                          "tumor_type": "NSCLC", "line_of_therapy": "first" },
  "transportability": { "validated_tumor_types": ["NSCLC"],
                        "out_of_context_action": "tier_down_to_D and warn" },
  "tier": "C", "review_status": "unverified", "primary_citation": "claret-2009-tgi"
}
```

> **Honesty note.** v0.1 parameter values are *illustrative* and `unverified` by
> design — see the [verification checklist](CONTRIBUTING.md). The infrastructure
> (schema, kernels, tier propagation, round-trip-validated exports) is real and
> tested; promoting records to `verified` from source PDFs is the
> highest-leverage contribution.

---

## Confidence tiers and propagation

| Tier | Meaning |
| --- | --- |
| **A** | Model + parameters externally validated; TGI→survival link held in ≥1 *independent* trial; broad context. |
| **B** | One robust model from a well-powered trial with at least a partial external check. |
| **C** | Single trial, narrow tumor type/line; no external validation; high-CV kill/resistance terms. |
| **D** | Transported outside its validated context, **or** hypothesis-tier (e.g. immuno-oncology). **Not predictive.** |

Two rules are enforced in code (`onkos/tiers.py`, tested in `tests/`):

1. **Worst input wins.** A composed simulation (`growth + drug_effect +
   resistance + exposure_response + survival_link`) inherits the worst component
   tier.
2. **Out-of-context transport forces a tier floor of D + a warning.** You cannot
   get an A-looking forecast from a model validated only on a different tumor
   type. This is what greys models out in the divergence view.

![Tier distribution](docs/images/tiers.png)

---

## Models & kernels

Every model binds to a **pure-NumPy/SciPy reference kernel** in
`onkos/export/reference.py`, the single computational ground truth. `E` is the
drug-effect magnitude that scales the kill term — supplied directly or derived
from a PK exposure through an exposure-response kernel (below).

| Kernel | Kind | Dynamics | Records |
| --- | --- | --- | --- |
| `growth_exponential` | ODE | `dV/dt = kg·V` | `growth_laws.exponential` |
| `growth_logistic` | ODE | `dV/dt = kg·V·(1 − V/Vmax)` | `growth_laws.logistic` |
| `growth_gompertz` | ODE | `dV/dt = kg·V·ln(Vmax/V)` | `growth_laws.gompertz` |
| `claret_tgi` | ODE | `dy/dt = kL·y − kD·E·e^(−λt)·y` (resistance = exp-decay of kill) | `resistance.claret_2009.tgi` |
| `biexp_tgi` | ODE | `y = y0·(e^(−ks·E·t) + e^(kg·t) − 1)` (shrink + regrowth) | `tgi_metrics.wang_2009.*`, `tgi_metrics.bruno_2020.*` |
| `survival_weibull_ph` | survival | `S(t) = exp(−(t/scale)^shape · e^(β·x))`, `x` = week-8 change | `survival_link.nsclc_os_week8` |
| `er_emax` | exposure-response | `E = Emax·C/(EC50+C)` | `exposure_response.emax_generic`, `…dacomitinib_egfr.emax` |
| `er_sigmoid_emax` | exposure-response | `E = Emax·C^γ/(EC50^γ+C^γ)` | `exposure_response.sigmoid_emax_generic` |
| `er_power` | exposure-response | `E = slope·C^θ` | `exposure_response.power_generic` |
| `simeoni_exp_linear` | ODE | `dw/dt = λ0·w / (1+(λ0·w/λ1)^ψ)^(1/ψ)` (exp→linear) | `growth_laws.simeoni_exp_linear` |
| `simeoni_tgi` | ODE (4-state) | transit-chain TGI; observe `w = x1+x2+x3+x4` | `preclinical_translation.simeoni_2004.xenograft` |
| `ivive_power` | exposure-response | `potency = scale·IC50^power` | `preclinical_translation.ivive_potency` |
| `io_tumor_immune` | ODE (2-state) | Kuznetsov tumor–immune predator-prey (**hypothesis-tier**) | `immuno_oncology.kuznetsov_1994.tumor_immune`, `…poorly_immunogenic.hypothesis` |

---

## Exposure-response & PK composability (Phase B)

The exposure-response (ER) layer maps a PK exposure metric `C` (C_avg, AUC,
C_max) to the drug-effect magnitude `E` that drives a TGI model's kill term. This
makes the **potency** of a regimen first-class (with its own tier and IIV) and
completes the chain **PK → exposure → tumor dynamics → survival** — the seam
where a [Hypnos](#licensing--citation) PK record composes with an Onkos TGI
model. A *time-varying* exposure (a full PK profile aligned to `t`) yields a
time-varying `E(t)`, and the tumor ODE is integrated numerically; a scalar
exposure uses the fast closed form.

![Exposure-response and PK-driven tumor dynamics](docs/images/exposure_response.png)

```python
import numpy as np, onkos
ds = onkos.load()
ctx = dict(tumor_type="NSCLC", line="first")

# Scalar exposure -> Emax transform -> drug effect -> Claret TGI -> OS
traj = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                      exposure=200.0,                                   # C_avg in µg/L
                      exposure_response="exposure_response.dacomitinib_egfr.emax")

# Time-varying PK profile (e.g. piped from Hypnos) -> E(t) -> ODE integration
t = np.linspace(0, 104, 209)
C = 300.0 * np.exp(-0.02 * t)                                          # declining exposure
traj = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                      exposure=C, exposure_response="exposure_response.emax_generic", t=t)
```

```text
$ onkos simulate resistance.claret_2009.tgi \
    --exposure 200 --exposure-response exposure_response.dacomitinib_egfr.emax
resistance.claret_2009.tgi  tier=C  (exposure=200.0 via exposure_response.dacomitinib_egfr.emax)
```

The ER record's tier and `transportability` propagate like any other component:
an ER model validated only on NSCLC/EGFR-TKI floors an out-of-context simulation
to **D** with a warning, exactly as the TGI and survival components do.

---

## Tumor-context library (Phase C)

The divergence view is only broadly useful if it has a context to run in. Phase C
builds the `tumor_type_baselines` library and the matching per-context survival
links, so every supported tumor type carries:

- a **baseline** (`tumor_type_baselines.*`) — baseline SLD `y0` and unperturbed
  growth, supplying the simulation's initial conditions;
- a **survival link** (`survival_link.*_os_week8`) — a tumor-specific Weibull-PH
  OS model whose scale reflects that indication's baseline prognosis;
- **≥2 eligible TGI models** (a Claret resistance form + a biexponential form),
  so model-selection risk is measurable rather than hypothetical.

| Context (1L) | baseline SLD | OS scale (wk) | eligible TGI models | OS divergence |
| --- | --- | --- | --- | --- |
| NSCLC | 80 mm | 60 | Claret 2009 · Wang 2009 biexp | 0.25 |
| breast | 55 mm | 130 | breast Claret · Bruno 2020 biexp | 0.16 |
| CRC | 90 mm | 95 | CRC Claret (capecitabine) · CRC biexp | 0.26 |
| HCC | 110 mm | 48 | HCC Claret · HCC biexp | 0.33 |
| melanoma | 60 mm | 85 | melanoma Claret · melanoma biexp | 0.25 |

![Tumor-context library](docs/images/context_library.png)

Each context resolves its own baseline and survival link automatically; a model
from one tumor type applied to another is greyed out (floored to **D**) by the
same transportability rule. Values are illustrative and `unverified` by design.

```python
import onkos
ds = onkos.load()
for tt in ["NSCLC", "breast", "CRC", "HCC", "melanoma"]:
    cmp = onkos.compare(ds, purpose="tgi", context=dict(tumor_type=tt, line="first"))
    print(tt, len(cmp.included), round(cmp.os_divergence, 2))
```

---

## Preclinical translation (Phase D)

The discovery-to-clinic bridge. Onkos implements the canonical **Simeoni 2004**
xenograft PK/PD model — the project's first *multi-state* ODE system. Unperturbed
growth is exponential then linear; drug at concentration `E` damages
proliferating cells (`x1`) at rate `k2·E`, and damaged cells traverse a
**signal-distribution transit chain** `x2→x3→x4` (rate `k1`) before dying, which
produces the characteristic *delayed* cell death. The observed tumor weight is
`w = x1+x2+x3+x4`.

![Preclinical Simeoni model](docs/images/preclinical.png)

```
dx1/dt = λ0·x1 / (1+(λ0·w/λ1)^ψ)^(1/ψ) − k2·E·x1      (proliferating)
dx2/dt = k2·E·x1 − k1·x2                               (damaged, transit)
dx3/dt = k1·x2  − k1·x3
dx4/dt = k1·x3  − k1·x4                  w = x1+x2+x3+x4 (observed weight)
```

Multi-state systems have no closed form, so the kernel framework integrates them
numerically and validates exports state-by-state: the SBML round-trip re-parses
**each** rate rule's MathML and checks it against the reference `rhs`, and the
NONMEM stream emits one `$DES` compartment per state. A concentration profile can
drive the kill term directly (`exposure=...`), so a Hypnos PK curve composes here
too.

```python
# Dose-dependent xenograft TGI (concentration drives the kill term directly)
tr = onkos.simulate(ds, "preclinical_translation.simeoni_2004.xenograft",
                    context=dict(tumor_type="ovarian_xenograft"), drug_effect=120.0)
tr.tumor_size                      # total tumor weight w(t)
tr.os_curve                        # None — preclinical models carry no survival link
```

**In-vitro → in-vivo translation.** `preclinical_translation.ivive_potency` maps
an in-vitro potency (e.g. IC50) to an in-vivo potency parameter
(`potency = scale·IC50^power`). The assumption that in-vitro potency predicts
in-vivo activity is itself what must be validated (Rocchetti 2007), so the record
is tiered and annotated accordingly. Preclinical records are **excluded from the
clinical divergence view** and applying xenograft parameters to a human tumor
floors the result to **D** — the translation gap, made explicit.

---

## Immuno-oncology (Phase E) — represented honestly, not predictively

> 🛑 **HYPOTHESIS-TIER. NOT FOR PREDICTION.** The immuno-oncology subsystem ships
> **tier D by construction** because the quantitative validation to do otherwise
> honestly does not yet exist (spec §2, §3, §5, §10).

Onkos includes the **Kuznetsov 1994** tumor–immune QSP model — a 2-state
predator-prey system (tumor + effector cells) that reproduces the field's
qualitative regimes: immune **control / dormancy**, immune **escape**, and the
bistable **rescue** when an immunotherapy effect (e.g. checkpoint blockade)
pushes the system across its threshold.

```
d tumor/dt    = α·tumor·(1 − β·tumor) − (1+E)·effector·tumor
d effector/dt = s + ρ·effector·tumor/(η+tumor) − μ·effector·tumor − δ·effector
```

![Tumor–immune QSP (hypothesis-tier)](docs/images/immuno_oncology.png)

The non-predictive stance is enforced in code, not just documented:

- the **validator rejects** any immuno-oncology record (or parameter) that is not
  tier D — `onkos validate` fails otherwise;
- every export carries `onkos:predictionStatus = "DO NOT USE FOR PREDICTION
  (hypothesis-tier)"` and the virtual-trial JSON sets `DO_NOT_USE_FOR_PREDICTION: true`;
- IO models are **excluded from the clinical divergence view** and never receive
  a survival link (no OS curve).

This is the Nidus "Phase-C" convention: the frontier is represented so it can be
explored and exported, but it can never masquerade as validated.

---

## Dataset health & releasing (Phase F)

`onkos report` turns the dataset's own honesty fields into a machine-generated
health report ([`docs/dataset-health.md`](docs/dataset-health.md)), kept in sync
with the data by a CI gate. It surfaces tier and review-status coverage, the
external-validation backlog, and the hypothesis-tier records.

![Dataset health coverage](docs/images/coverage.png)

```text
$ onkos report | head
# Onkos dataset health report
- Records: 33  ·  Citations: 8
- Verified (PDF-checked): 0 / 33
- External-validation coverage (clinical TGI + survival models): 15 / 15  100%
```

**Releasable, proven in CI.** The dataset is the source of truth at the repo
root; `scripts/sync_dataset_into_package.py` copies it into the package as
`_dataset/` for packaging. The CI `release` job builds the wheel, installs it
into a clean environment, and runs `onkos validate / info / report` and a
simulation **from outside the repository** — proving the bundled dataset ships
and resolves without a source checkout. Resolution order puts the source
`dataset/` first so a stale synced copy can never shadow live edits during
development; the bundled `_dataset/` is the wheel-only fallback.

Release metadata: [`CHANGELOG.md`](CHANGELOG.md), [`.zenodo.json`](.zenodo.json)
(Zenodo concept DOI on first deposit), [`CITATION.cff`](CITATION.cff), and a
`py.typed` marker so downstream type-checkers see Onkos's annotations.

---

## Architecture

The **dataset is the single source of truth**; everything else is a
deterministic projection.

```mermaid
flowchart TD
    DS["<b>dataset/</b> — source of truth<br/>JSON records + JSON Schema + JSON-LD context<br/>models · params · derivation context · transportability · tiers · citations"]
    DS -->|sync_dataset_into_package.py| PKG["<b>onkos</b> package<br/>load · filter · validate · simulate · compare"]
    PKG --> CLI["<b>onkos</b> CLI"]
    PKG --> DASH["Streamlit dashboard<br/>browse + virtual-trial divergence"]
    PKG --> NB["Notebooks<br/>executed in CI (nbmake)"]
    PKG --> EXP["<b>onkos.export</b>"]
    EXP --> NM["NONMEM"]
    EXP --> SBML["SBML L3v2"]
    EXP --> PHARMML["PharmML"]
    EXP --> RX["rxode2 / Pumas"]
    EXP --> VT["virtual-trial JSON"]
    EXP --> OMEX["COMBINE .omex"]
```

```mermaid
flowchart LR
    P["dataset/records/*.json"] --> REG["registry.py<br/>bind record → kernel"]
    REG --> REF["reference.py<br/>NumPy/SciPy kernels"]
    REG --> B["nonmem / sbml / pharmml / rxode2 / pumas"]
    ANN["annotate.py<br/>clinicalUse=PROHIBITED · tier · DOI RDF"] --> B
    REF -. "round-trip validates (1e-6 algebraic / 1e-4 ODE)" .-> B
```

### Round-trip validation — why exports cannot lie

Each ODE kernel declares three *independent* expressions of the same dynamics: a
closed-form `analytic` solution, a hand-written `rhs`, and an `rhs_infix` string.
CI checks ([`tests/test_roundtrip.py`](tests/test_roundtrip.py)):

- **analytic vs. SciPy ODE integration** → agreement to ~1e-4 (single-state
  closed forms; validates the rhs);
- **SBML re-parsed**: the generated MathML rate law is converted back to an
  expression and evaluated against `rhs` → ~1e-6, **per state** (so the
  multi-state Simeoni system is checked compartment-by-compartment);
- **NONMEM re-parsed**: `$THETA` initial estimates must equal the dataset values,
  and one `$DES` compartment is emitted per state;
- **rxode2 / Pumas / SO re-parsed**: the parameter vector is read back from each
  and must equal the dataset values;
- **cross-format consistency**: NONMEM, SBML, PharmML-SO, rxode2, and Pumas must
  all agree on the parameter values — one source of truth, five renderings.

Multi-state kernels (Simeoni) have no closed form, so the analytic check is
skipped for them and the rhs is instead pinned by the per-state SBML round-trip
plus behavioral tests (exp→linear growth, dose-dependent shrinkage, transit
delay). An export bug therefore cannot ship silently.

### Design decisions

| Decision | Rationale |
| --- | --- |
| Pure Python (NumPy/SciPy); R/Julia only as export targets | Nothing is compute-bound; R/Julia models are generated artifacts, not runtime deps. |
| Dataset is the centerpiece; everything else is presentation | The durable contribution is the curated, tiered, context-annotated parameters. |
| `derivation_context` + `transportability` are first-class | Out-of-context transport is the dominant silent error; machine-enforcing it is the load-bearing idea. |
| IIV CV surfaced on kill/resistance terms | A ~90%-CV term must not present as a point estimate. |
| Tiers + transport warnings propagate; worst input wins | A forecast is only as trustworthy as its least-validated component. |
| Population-level forward simulation only | The line between research tool and clinical tool is exactly individual prediction. Onkos stays on the safe side by construction. |
| Exposure-response is a separate, tiered kernel (not baked into the TGI model) | Potency/uncertainty are drug-specific and reusable; decoupling them lets one ER record drive many TGI models and keeps the PK→effect seam explicit and tier-propagating. |
| Scalar exposure uses the closed form; time-varying PK integrates the ODE | Exactness and speed for the common case; correctness for a full PK profile, where the constant-E closed form would be wrong. |
| Multi-state kernels keep `analytic` optional; an `observable` maps states → the measured quantity | The Simeoni transit model has no closed form. Numerical integration + a per-state SBML round-trip preserve export-correctness guarantees without forcing a closed form; the observable (total weight = Σ compartments) decouples the measured signal from the latent states. |
| Hypothesis-tier (immuno-oncology) is enforced in code, not just documented | A "do not predict" note in prose is easy to ignore. The validator fails the build if an IO record is not tier D, exports carry a machine-readable `predictionStatus`, and IO is excluded from the clinical view — so the frontier can be explored but never masquerade as validated. |
| Parameter uncertainty is propagated, not just stored | Storing `iiv_cv_percent` but simulating on central values would let a ~90%-CV term pose as a point estimate. `simulate_ensemble` samples IIV lognormally (median-preserving) so the reported variability flows into tumor/OS bands — the second uncertainty axis alongside model-selection divergence. |
| TGI metrics are extracted model-agnostically (Stein method), not read from params | Reading k_g/k_s off a record only works for the biexponential; the Claret/Simeoni structures have no such params. Extracting them from the trajectory the way Stein extracts them from RECIST data makes the metric panel uniform across kernels — and recovers the generating rates as a built-in correctness check. |
| Sensitivity uses independent sampling so first-order indices are correlations | Sampling each IIV parameter independently makes the standardized regression coefficient equal the input-target correlation and the squared SRCs partition the variance — a first-order Sobol decomposition with no extra design. It also exposes that CV alone ≠ influence (influence is CV × effect-strength), pointing verification at the parameter that actually moves the prediction. |
| PharmML SO carries IIV as random-effect variance, never as estimate precision (RSE) | The SO's job is to report results, but the dataset curates inter-individual variability, not the precision of the population estimate. Encoding IIV as `omega = ln(1+CV²)` is faithful; fabricating an RSE we don't have would not be — so the precision block is deliberately omitted. |
| Composable with Hypnos | A shared export/annotation convention lets a Hypnos PK record drive an Onkos TGI model end to end via an exposure-response record. |

---

## Repository layout

```
onkos/
├── dataset/                     # SOURCE OF TRUTH
│   ├── schema/                  # JSON Schema + JSON-LD context
│   ├── records/                 # one JSON per model / context-baseline
│   └── citations/               # Crossref/PubMed citation records
├── python/onkos/
│   ├── load · filter · validate · tiers · simulate · metrics · compare · uncertainty · sensitivity · report · cli
│   ├── py.typed                 # PEP 561 typing marker
│   └── export/                  # registry · reference · nonmem · sbml · pharmml · pharmml_so
│       · rxode2 · pumas · virtual_trial_json · combine · annotate
├── dashboard/app.py             # Streamlit: browse + divergence view
├── notebooks/                   # executed in CI (nbmake)
├── scripts/                     # sync_dataset_into_package · make_figures
├── tests/                       # schema · simulate · round-trip · CLI · report · …
├── docs/                        # essay · specs/v0.1/spec.md · dataset-health.md · images/
├── CHANGELOG.md · CITATION.cff · .zenodo.json   # release metadata
└── .github/workflows/ci.yml     # lint · test (3.9–3.12) · exports · releasable wheel
```

---

## Scope & safety

**In scope:** unperturbed growth laws; drug-effect/kill models; resistance/
regrowth (the λ term); exposure-response links; TGI-derived metrics; TGI-metric →
survival models; tumor-type/line baselines; a separated preclinical-translation
subsystem; immuno-oncology *only* as a hypothesis-tier, non-predictive subsystem.

**Out of scope (hard line, not a roadmap item):** any per-patient prognosis,
survival estimate for a real person, treatment recommendation, or therapy
ranking. The tell that the project has crossed its line is any feature that takes
a real patient's tumor measurement and returns a prognosis or a therapy choice.
**That feature does not get built.** See [spec §10](docs/specs/v0.1/spec.md).

---

## Roadmap

| Phase | Content | Status |
| --- | --- | --- |
| **A — TGI spine** | Growth laws + Claret TGI + NSCLC context + TGI→OS link + divergence view; NONMEM + SBML; round-trip validation. | ✅ v0.1 |
| **B — Resistance + exposure-response** | Emax / sigmoid-Emax / power ER kernels driving the kill term; scalar **and** time-varying PK-driven simulation (Hypnos composability); ER tier + transportability propagation; PharmML + rxode2/Pumas; IIV-CV surfaced. | ✅ v0.2 |
| **C — Survival + baselines** | `tumor_type_baselines` library + per-context Weibull-PH survival links across NSCLC, breast, CRC, HCC, melanoma; ≥2 eligible TGI models per context; cross-context divergence; orphan-record invariant enforced in CI. | ✅ v0.3 |
| **D — Preclinical translation** | Multi-state ODE framework; Simeoni 2004 xenograft model (exp→linear growth + signal-distribution transit chain); in-vitro→in-vivo potency translation; per-state SBML/NONMEM export + round-trip. | ✅ v0.4 (this release) |
| **E — Immuno-oncology** | Kuznetsov tumor–immune QSP, hypothesis-tier (tier D), non-predictive; tier-D enforced by the validator; DO-NOT-PREDICT annotation on every export; excluded from the clinical view. | ✅ v0.5 (this release) |
| **F — Hardening** | External-validation backfill (coverage 15/15); `onkos report` health report with CI sync gate; wheel-build releasability proven in CI; `.omex`, `.zenodo.json`, `CHANGELOG.md`, `py.typed`, `CITATION.cff`. | ✅ v0.6 (this release) |

The phased roadmap (spec §11, Phases A–F) is now fully implemented. Remaining work
is **breadth and verification**: promoting `unverified` records to `verified` from
source PDFs, and adding more drugs / tumor types / lines (see CONTRIBUTING.md).

---

## Licensing & citation

- **Code:** MIT ([LICENSE](LICENSE)).
- **Dataset:** CC-BY-4.0 ([LICENSE-DATASET](LICENSE-DATASET)).
- **Citation:** [`CITATION.cff`](CITATION.cff). When you use a record, cite Onkos
  **and** the original source via `record.primary_citation.doi`.

Sibling projects: **Nidus** (gestational physiology, per-parameter tier) and
**Hypnos** (anesthetic PK/PD, applicability envelope). Hypnos and Onkos compose:
a Hypnos PK record can drive the exposure-response of an Onkos TGI model, giving
an open, tier-annotated PK → exposure → tumor-dynamics → survival chain.
