"""Dataset health & validation report — the Phase-F hardening surface.

Turns the dataset's own honesty fields (tier, review status, external-validation
performance, hypothesis-tier flags) into a machine-generated Markdown report so
coverage is visible and tracked rather than asserted.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from .load import Dataset

# Records for which an external OS-validation metric is meaningful: clinical TGI
# and survival-link models. Growth laws, baselines, transforms, preclinical, and
# hypothesis-tier IO are excluded by construction.
_VALIDATION_ELIGIBLE_PURPOSES = {"tgi", "survival_link"}
_NON_CLINICAL_SUBSYSTEMS = {"immuno_oncology", "preclinical_translation"}


def validation_eligible(ds: Dataset):
    return [
        r
        for r in ds
        if r.purpose in _VALIDATION_ELIGIBLE_PURPOSES
        and r.subsystem not in _NON_CLINICAL_SUBSYSTEMS
    ]


def external_validation_coverage(ds: Dataset):
    """(n_validated, n_eligible, fraction) for clinical TGI / survival models."""
    eligible = validation_eligible(ds)
    validated = [r for r in eligible if r.predictive_performance]
    n = len(eligible)
    return len(validated), n, (len(validated) / n if n else 0.0)


def stats(ds: Dataset) -> dict:
    n_val, n_elig, frac = external_validation_coverage(ds)
    return {
        "version": ds.version,
        "n_records": len(ds),
        "n_citations": len(ds.citations),
        "by_subsystem": dict(sorted(ds.by_subsystem().items())),
        "by_tier": {k: ds.by_tier().get(k, 0) for k in ("A", "B", "C", "D")},
        "by_review_status": dict(sorted(ds.by_review_status().items())),
        "tier_by_subsystem": _tier_by_subsystem(ds),
        "external_validation": {"validated": n_val, "eligible": n_elig, "fraction": frac},
        "hypothesis_tier": [r.id for r in ds if r.subsystem == "immuno_oncology"],
        "verified": sum(1 for r in ds if r.review_status == "verified"),
    }


# Model-selection-uncertainty bins. Reported qualitatively (not as raw Monte-Carlo
# floats) so the CI report-in-sync diff stays byte-stable across platforms; the
# context fractions sit comfortably away from these bin edges.
_MSU_HORIZON_WEEKS = 156.0
_MSU_N = 120
_MSU_LOW, _MSU_HIGH = 0.25, 0.45


def _msu_label(fraction: float) -> str:
    if fraction >= _MSU_HIGH:
        return "high"
    if fraction >= _MSU_LOW:
        return "moderate"
    return "low"


def model_selection_summary(ds: Dataset) -> list[dict]:
    """Per-clinical-context model-selection-uncertainty summary (research spec §10,
    step 4). For each tumor-type/line context with ≥2 eligible TGI models, the
    irreducible model-choice fraction is binned (low/moderate/high) so contexts
    can be ranked by where adding a better-validated model has the most value.

    Deferred import of :mod:`onkos.compare` avoids an import cycle (compare →
    combine → uncertainty → simulate, none of which import report)."""
    from .compare import compare

    contexts = sorted(
        {
            (r.derivation_context.tumor_type, r.derivation_context.line_of_therapy)
            for r in ds
            if r.kind == "context_baseline"
            and r.derivation_context
            and r.derivation_context.tumor_type
            and r.derivation_context.line_of_therapy
        }
    )
    t = np.linspace(0.0, _MSU_HORIZON_WEEKS, int(2 * _MSU_HORIZON_WEEKS) + 1)
    out: list[dict] = []
    for tumor_type, line in contexts:
        cmp = compare(
            ds, purpose="tgi", context={"tumor_type": tumor_type, "line": line},
            drug_effect=1.0, t=t,
        )
        if len(cmp.included) < 2:  # a single model has no cross-model disagreement
            continue
        ma = cmp.model_average(target="median_os_weeks", endpoint="OS", weights="equal",
                               n=_MSU_N, seed=0)
        out.append({
            "tumor_type": tumor_type,
            "line": line,
            "n_models": len(cmp.included),
            "tier": ma.tier,
            "risk": _msu_label(ma.model_selection_fraction),
        })
    # Rank high → moderate → low, then alphabetically — a deterministic order
    # that does not depend on the raw (platform-sensitive) fraction.
    rank = {"high": 0, "moderate": 1, "low": 2}
    out.sort(key=lambda d: (rank[d["risk"]], d["tumor_type"], d["line"]))
    return out


# Practical-identifiability reference design (research spec practical-identifiability
# §6): a realistic RECIST scan cadence + 20% proportional assay error. The verdict
# reported is the binary `practically_identifiable` bool (RSE < 50% on every
# parameter AND collinearity index < 15) plus the least-identifiable parameter — no
# raw floats, so the CI report-in-sync diff stays byte-stable. The eligible models
# bifurcate cleanly (worst RSE either <=47% or >=127%), well clear of the boundary.
_PI_SIGMA_PROP = 0.2


def practical_identifiability_summary(ds: Dataset) -> list[dict]:
    """Per-clinical-TGI-model practical-identifiability under a fixed reference
    design: can a realistic trial even estimate the model's structural parameters,
    or is its precision (and the stored IIV CV) a flat-likelihood artifact? Ranks the
    models whose estimates a realistic trial cannot support to the top — the
    design-level analog of the sensitivity / model-selection curation triage.

    Deferred import avoids loading the simulation stack at report import time."""
    from .export.registry import get_kernel
    from .identify import identifiability

    out: list[dict] = []
    for r in ds:
        if r.purpose not in ("tgi", "metric") or r.subsystem == "immuno_oncology":
            continue
        if r.kernel is None or get_kernel(r).kind != "ode":
            continue
        dc = r.derivation_context
        ctx = {
            "tumor_type": dc.tumor_type if dc else None,
            "line": dc.line_of_therapy if dc else None,
        }
        res = identifiability(ds, r.id, context=ctx, sigma_prop=_PI_SIGMA_PROP)
        worst = res.worst
        out.append({
            "record_id": r.id,
            "tier": res.tier,
            "n_params": len(res.params),
            "identifiable": res.practically_identifiable,
            "worst": worst.symbol if worst else "-",
        })
    # Not-identifiable first (the triage priority), then alphabetically — a
    # deterministic order independent of the raw (platform-sensitive) RSE values.
    out.sort(key=lambda d: (d["identifiable"], d["record_id"]))
    return out


# Model-selection-budget reference depth. The per-context structural fraction is
# reported as a binary (structure- vs parameter-dominated) split at 0.5 — the eligible
# contexts cluster near ~0.40 and ~0.65, well clear of the boundary, so the CI
# report-in-sync diff stays byte-stable.
_BUDGET_N = 80


def model_selection_budget_summary(ds: Dataset) -> list[dict]:
    """Per-context model-selection budget (research spec model-selection-budget §6,
    step 4): the share of a composed OS forecast's variance that is irreducible
    structural-choice risk (TGI model + survival link + their interaction) vs
    parameter noise a bigger trial would shrink, plus whether the survival-model axis
    is even cross-checked (≥2 eligible survival links). Ranks where standardizing an
    assumption has the most leverage.

    Deferred import avoids loading the budget stack at report import time."""
    from .budget import eligible_survival_links, model_selection_budget

    contexts = sorted(
        {
            (r.derivation_context.tumor_type, r.derivation_context.line_of_therapy)
            for r in ds
            if r.kind == "context_baseline"
            and r.derivation_context
            and r.derivation_context.tumor_type
            and r.derivation_context.line_of_therapy
        }
    )
    out: list[dict] = []
    for tumor_type, line in contexts:
        ctx = {"tumor_type": tumor_type, "line": line}
        try:
            b = model_selection_budget(ds, context=ctx, endpoint="OS", n=_BUDGET_N, seed=0)
        except ValueError:
            continue  # no eligible TGI models / links for this context
        out.append({
            "tumor_type": tumor_type,
            "line": line,
            "n_models": len(b.models),
            "n_links": len(eligible_survival_links(ds, ctx, "OS")),
            "tier": b.tier,
            # Binary, edge-safe at 0.5 (contexts cluster near 0.40 / 0.65).
            "dominated_by": "structure" if b.structural_fraction >= 0.5 else "parameter",
        })
    # Structure-dominated first (where standardization buys the most), then alphabetical.
    out.sort(key=lambda d: (d["dominated_by"] != "structure", d["tumor_type"], d["line"]))
    return out


def _tier_by_subsystem(ds: Dataset) -> dict:
    table: dict = {}
    for r in ds:
        table.setdefault(r.subsystem, Counter())[r.tier] += 1
    return {sub: dict(c) for sub, c in sorted(table.items())}


def _bar(frac: float, width: int = 20) -> str:
    filled = round(frac * width)
    return "█" * filled + "░" * (width - filled)


def build_report(ds: Dataset) -> str:
    s = stats(ds)
    ev = s["external_validation"]
    lines = [
        "# Onkos dataset health report",
        "",
        "> GENERATED by `onkos report` — do not hand-edit. NOT FOR CLINICAL USE.",
        "",
        f"- **Dataset version:** {s['version']}",
        f"- **Records:** {s['n_records']}  ·  **Citations:** {s['n_citations']}",
        f"- **Verified (PDF-checked):** {s['verified']} / {s['n_records']}",
        f"- **External-validation coverage** (clinical TGI + survival models): "
        f"{ev['validated']} / {ev['eligible']}  `{_bar(ev['fraction'])}` {ev['fraction'] * 100:.0f}%",
        "",
        "## Records by subsystem",
        "",
        "| subsystem | records |",
        "| --- | --- |",
    ]
    lines += [f"| {k} | {v} |" for k, v in s["by_subsystem"].items()]

    lines += ["", "## Confidence tiers", "", "| tier | records | meaning |", "| --- | --- | --- |"]
    meaning = {
        "A": "externally validated, broad context",
        "B": "robust single trial, partial external check",
        "C": "single trial, no external validation",
        "D": "out-of-context transport, or hypothesis-tier (not predictive)",
    }
    lines += [f"| {t} | {s['by_tier'][t]} | {meaning[t]} |" for t in ("A", "B", "C", "D")]

    lines += ["", "## Tier × subsystem", "", "| subsystem | A | B | C | D |", "| --- | --- | --- | --- | --- |"]
    for sub, counts in s["tier_by_subsystem"].items():
        lines.append(
            f"| {sub} | {counts.get('A', 0)} | {counts.get('B', 0)} | "
            f"{counts.get('C', 0)} | {counts.get('D', 0)} |"
        )

    lines += ["", "## Review status", "", "| status | records |", "| --- | --- |"]
    lines += [f"| {k} | {v} |" for k, v in s["by_review_status"].items()]

    missing = [r.id for r in validation_eligible(ds) if not r.predictive_performance]
    lines += ["", "## External-validation backlog", ""]
    if missing:
        lines.append("Clinical models without a recorded external-validation metric:")
        lines += [f"- `{rid}`" for rid in sorted(missing)]
    else:
        lines.append("All eligible clinical models carry an external-validation metric.")

    # Evidence-based tier audit (spec §5: tiers are partly numeric).
    from .audit import audit_tiers

    findings = audit_tiers(ds)
    inflated = [f for f in findings if f.status == "inflated"]
    conservative = [f for f in findings if f.status == "conservative"]
    lines += [
        "",
        "## Evidence-based tier audit",
        "",
        f"- **Tier inflation** (assigned tier exceeds recorded evidence): {len(inflated)}"
        + (" ✅" if not inflated else ""),
        f"- **Conservative** (external validation recorded; could upgrade if trusted): "
        f"{len(conservative)} / {len(findings)} clinical TGI / survival records",
    ]
    if inflated:
        lines.append("\nInflated records (MUST be corrected):")
        lines += [f"- `{f.record_id}`: tier {f.assigned} > ceiling {f.ceiling}" for f in inflated]

    msu = model_selection_summary(ds)
    if msu:
        lines += [
            "",
            "## Model-selection uncertainty by context",
            "",
            "Of everything uncertain in a composed OS forecast, how much is *irreducible "
            "model-choice risk* (between-model disagreement that more data on any one model "
            "cannot resolve) versus estimable parameter noise? High-risk contexts are where "
            "adding a better-validated TGI model has the most value (curation triage). "
            "Equal-weight combination; population/trial level only. See `onkos.combine`.",
            "",
            "| context | line | eligible models | tier | model-selection risk |",
            "| --- | --- | --- | --- | --- |",
        ]
        lines += [
            f"| {m['tumor_type']} | {m['line']} | {m['n_models']} | {m['tier']} | {m['risk']} |"
            for m in msu
        ]

    pid = practical_identifiability_summary(ds)
    if pid:
        n_unident = sum(1 for m in pid if not m["identifiable"])
        lines += [
            "",
            "## Practical identifiability by model",
            "",
            "Under a realistic RECIST scan cadence (weeks 0, 6, 12, 18, 24, 36, 48; 20% "
            "proportional assay error), could a trial of that shape even *estimate* each "
            "model's structural parameters? A model flagged *not identifiable* has a "
            "parameter whose predicted relative standard error exceeds 50% (or a collinear, "
            "non-separable parameter combination) — so its reported point value, and often "
            "the large IIV CV travelling with it, are partly a flat-likelihood artifact of "
            "the design rather than clean estimates. These are the models where a richer "
            "trial design or an external constraint is needed before an estimate should be "
            "reused (curation triage). Design level only; identifiability cannot move a "
            "tier. See `onkos.identify`.",
            "",
            f"- **Not practically identifiable** under the reference design: {n_unident} / "
            f"{len(pid)} clinical TGI models.",
            "",
            "| model | params | tier | identifiable? | least-identifiable |",
            "| --- | --- | --- | --- | --- |",
        ]
        for m in pid:
            verdict = "yes" if m["identifiable"] else "**no**"
            worst = "—" if m["identifiable"] else f"`{m['worst']}`"
            lines.append(
                f"| `{m['record_id']}` | {m['n_params']} | {m['tier']} | {verdict} | {worst} |"
            )

    budget = model_selection_budget_summary(ds)
    if budget:
        n_struct = sum(1 for b in budget if b["dominated_by"] == "structure")
        lines += [
            "",
            "## Model-selection budget by context",
            "",
            "Splitting a composed OS forecast's total variance across its structural choices "
            "(TGI model + survival link + their interaction) versus parameter noise a bigger "
            "trial would shrink. A context *dominated by structure* is one where most of the "
            "forecast uncertainty is irreducible model-choice risk — standardizing an "
            "assumption (or adding a second survival link where there is only one, so the "
            "survival-model axis is even cross-checked) has more leverage there than more "
            "patients. Equal-weight balanced design; population/trial level only. See "
            "`onkos.budget`.",
            "",
            f"- **Structure-dominated** (structural share ≥ 50% of forecast variance): "
            f"{n_struct} / {len(budget)} contexts.",
            "",
            "| context | line | TGI models | survival links | tier | uncertainty dominated by |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for b in budget:
            links = b["n_links"] if b["n_links"] > 1 else f"{b['n_links']} ⚠"
            lines.append(
                f"| {b['tumor_type']} | {b['line']} | {b['n_models']} | {links} | {b['tier']} | "
                f"{b['dominated_by']} |"
            )

    if s["hypothesis_tier"]:
        lines += [
            "",
            "## Hypothesis-tier (NOT FOR PREDICTION)",
            "",
            "These records ship tier D by construction and are excluded from the clinical view:",
        ]
        lines += [f"- `{rid}`" for rid in s["hypothesis_tier"]]

    lines += [
        "",
        "## Honesty note",
        "",
        "All v0.x parameter values are illustrative and `unverified` by design. Promoting "
        "records to `verified` by reading the source PDF (especially the derivation-context "
        "and transportability claims) is the highest-leverage contribution — see CONTRIBUTING.md.",
        "",
    ]
    return "\n".join(lines)
