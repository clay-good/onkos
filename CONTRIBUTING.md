# Contributing to Onkos

Thank you for helping build an honest ground-truth layer for oncology TGI
modeling. The single highest-leverage contribution is **promoting `unverified`
records to `verified` by reading the source PDF** — with the oncology-specific
twist that the **derivation context and transportability claims** deserve the
most scrutiny, because that is where over-broad reuse originates.

## The dataset is the source of truth

`dataset/records/*.json` is canonical. Everything else (the Python package, CLI,
dashboard, exports) is a deterministic projection. Never hand-edit an export;
regenerate it (`onkos export ...`, `scripts/make_figures.py`).

## Adding or editing a record

1. Add a JSON file under `dataset/records/`. The **filename must equal the
   record `id`** (e.g. `resistance.claret_2009.tgi.json`).
2. Conform to `dataset/schema/record.schema.json`. Run `onkos validate`.
3. The **record-level tier is the worst contributing parameter tier.** A test
   enforces this (`tests/test_dataset.py`).
4. Every citation key must resolve to a file in `dataset/citations/`.
5. If the record binds a reference kernel, set `kernel` to a name in
   `onkos.export.reference.KERNELS`. New dynamics need a new kernel with an
   `analytic` solution, a hand-written `rhs`, and an `rhs_infix` string so the
   round-trip tests can validate exports.

## Confidence tiers (A/B/C/D)

| Tier | Meaning |
| --- | --- |
| A | Model + parameters externally validated; TGI to survival link held in >=1 independent trial; broad context. |
| B | One robust model from a well-powered trial with at least a partial external check. |
| C | Single trial, narrow tumor type/line; no external validation; high-CV kill/resistance terms. |
| D | Transported outside its validated context, or hypothesis-tier (e.g. immuno-oncology mechanistic params). **Not predictive.** |

Rules that propagate automatically (do not bypass them):

- **Worst input wins.** A composed simulation inherits the worst component tier.
- **Out-of-context transport forces a tier floor of D + a warning.** Set each
  record's `transportability.validated_tumor_types` / `validated_lines` honestly.
- **Surface uncertainty.** Kill and resistance terms carry `iiv_cv_percent`; a
  term known only to ~90% CV must not present as a point estimate.

## PDF-verification checklist (to set `review_status: verified`)

Open the source and confirm, field by field:

1. **Structure & parameterization** — the model form and the exact form of the
   resistance / exposure-response terms.
2. **Every parameter value, its units, and its IIV** — this is where
   transcription errors hide.
3. **Derivation context** — drug, tumor type, line, trial, n, measurement basis.
4. **Validated transportability boundary** — how far beyond its origin the fit
   was actually checked.

LLMs may assist but never promote a record on their own authority. The verified
count is reported honestly by `onkos info`.

## Development

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
ruff check python/ tests/ scripts/
pytest -q
pytest --nbmake notebooks/ -q
```

## Safety boundary (non-negotiable)

Onkos does population/trial-level forward simulation only. Any feature that takes
a real patient's tumor measurement and returns a prognosis or a therapy choice
is out of scope and will not be merged. See spec §10.
