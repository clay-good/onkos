"""Onkos command-line interface."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__, compare, load, sensitivity, simulate, simulate_ensemble
from .export.combine import build_omex
from .export.jsonld import to_jsonld
from .export.nonmem import to_nonmem
from .export.pharmml import to_pharmml
from .export.pharmml_so import to_pharmml_so
from .export.pumas import to_pumas
from .export.registry import get_kernel
from .export.rxode2 import to_rxode2
from .export.sbml import to_sbml
from .export.virtual_trial_json import to_virtual_trial_json
from .report import build_report
from .validate import validate_dataset

_TEXT_EXPORTERS = {
    "nonmem": (to_nonmem, ".mod"),
    "sbml": (to_sbml, ".xml"),
    "pharmml": (to_pharmml, ".pharmml"),
    "so": (to_pharmml_so, ".so.xml"),
    "rxode2": (to_rxode2, ".R"),
    "pumas": (to_pumas, ".jl"),
    "vt-json": (to_virtual_trial_json, ".json"),
    "jsonld": (to_jsonld, ".jsonld"),
}


def _cmd_version(_args) -> int:
    print(f"onkos {__version__}")
    return 0


def _cmd_validate(_args) -> int:
    errors = validate_dataset()
    if errors:
        for e in errors:
            print(f"FAIL: {e}", file=sys.stderr)
        print(f"\n{len(errors)} validation error(s).", file=sys.stderr)
        return 1
    ds = load()
    print(f"OK: {len(ds)} record(s) valid against schema.")
    return 0


def _cmd_report(args) -> int:
    md = build_report(load())
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(md)
        print(f"Wrote {out}")
    else:
        print(md)
    return 0


def _cmd_info(_args) -> int:
    ds = load()
    print(f"Onkos dataset version {ds.version} — {len(ds)} records\n")
    print("By subsystem:")
    for k, v in sorted(ds.by_subsystem().items()):
        print(f"  {k:<26} {v}")
    print("\nBy tier:")
    for k in ("A", "B", "C", "D"):
        print(f"  {k}  {ds.by_tier().get(k, 0)}")
    print("\nBy review status:")
    for k, v in sorted(ds.by_review_status().items()):
        print(f"  {k:<12} {v}")
    return 0


def _cmd_simulate(args) -> int:
    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    # When an exposure-response record is supplied, derive E from the exposure
    # metric; otherwise use the scalar drug effect.
    er = args.exposure_response
    exposure = args.exposure
    drug_effect = None if (er and exposure is not None) else args.drug_effect
    effect_desc = (
        f"exposure={exposure} via {er}" if (er and exposure is not None)
        else f"drug_effect={args.drug_effect}"
    )

    if args.compare:
        cmp = compare(
            ds, purpose="tgi", context=ctx, drug_effect=drug_effect,
            exposure=exposure, exposure_response=er,
        )
        print(f"Virtual-trial comparison — {ctx}, {effect_desc}\n")
        for tr in cmp.included:
            print(f"  [{tr.tier}] {tr.record_id:<48} median OS {tr.median_os}")
        for rid, reason in cmp.excluded:
            print(f"  [-] {rid:<48} EXCLUDED ({reason})")
        print(f"\n  OS divergence (max pointwise): {cmp.os_divergence:.3f}")
        print(f"  Median OS range: {cmp.median_os_range}")
        return 0

    if not args.record:
        print("error: provide a RECORD id or --compare", file=sys.stderr)
        return 2
    tr = simulate(
        ds, args.record, context=ctx, drug_effect=drug_effect,
        exposure=exposure, exposure_response=er,
    )
    print(f"{args.record}  tier={tr.tier}  ({effect_desc})")
    for k, v in tr.metrics.items():
        print(f"  {k:<24} {v:.4f}")
    if tr.median_os is not None:
        print(f"  median_os_weeks          {tr.median_os:.2f}")
    for w in tr.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_uncertainty(args) -> int:
    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    er, exposure = args.exposure_response, args.exposure
    drug_effect = None if (er and exposure is not None) else args.drug_effect
    ens = simulate_ensemble(
        ds, args.record, context=ctx, drug_effect=drug_effect, exposure=exposure,
        exposure_response=er, n=args.n, seed=args.seed,
    )
    lo, hi = ens.ci
    print(f"{args.record}  tier={ens.tier}  (n={ens.n}, {lo:g}-{hi:g}% bands)\n")
    print(f"  {'metric':<24} {'median':>10} {'lo':>10} {'hi':>10}")
    for k, v in ens.metrics.items():
        print(f"  {k:<24} {v['median']:>10.3f} {v['lo']:>10.3f} {v['hi']:>10.3f}")
    for w in ens.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_sensitivity(args) -> int:
    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    res = sensitivity(ds, args.record, context=ctx, target=args.target, n=args.n, seed=args.seed)
    print(
        f"{args.record}  target={res.target}  "
        f"(n={res.n_used}/{res.n}, first-order R^2={res.r_squared:.2f})\n"
    )
    print(f"  {'parameter':<12} {'IIV CV%':>8} {'SRC':>8} {'contribution':>14}")
    for p in res.indices:
        bar = "█" * round(p.contribution * 20)
        print(f"  {p.symbol:<12} {p.iiv_cv_percent:>8.0f} {p.src:>+8.3f} "
              f"{p.contribution * 100:>12.0f}%  {bar}")
    if res.dominant:
        print(f"\n  -> verify '{res.dominant.symbol}' first (drives the prediction).")
    return 0


def _write_csv(ds, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["record_id", "symbol", "label", "central", "units", "iiv_cv_percent", "tier", "citation"]
        )
        for r in ds:
            for p in r.parameters:
                cit = p.primary_citation.key if p.primary_citation else ""
                w.writerow([r.id, p.symbol, p.label, p.value.central, p.value.units,
                            p.iiv_cv_percent, p.tier, cit])


def _write_bibtex(ds, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n\n".join(c.bibtex() for c in ds.citations.values()) + "\n")


def _cmd_export(args) -> int:
    ds = load()
    fmt = args.format
    out = Path(args.output)

    if fmt == "csv":
        target = out if out.suffix == ".csv" else out / "parameters.csv"
        _write_csv(ds, target)
        print(f"Wrote {target}")
        return 0
    if fmt == "bibtex":
        target = out if out.suffix == ".bib" else out / "citations.bib"
        _write_bibtex(ds, target)
        print(f"Wrote {target}")
        return 0
    if fmt == "omex":
        out.mkdir(parents=True, exist_ok=True)
        n = 0
        for r in ds:
            if r.kernel is None:
                continue
            build_omex(r, str(out / f"{r.id}.omex"))
            n += 1
        print(f"Wrote {n} .omex archive(s) to {out}/")
        return 0

    builder, ext = _TEXT_EXPORTERS[fmt]
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for r in ds:
        if r.kernel is None:
            continue
        spec = get_kernel(r)
        if fmt in ("nonmem", "sbml", "rxode2", "pumas") and spec.kind != "ode":
            continue
        (out / f"{r.id}{ext}").write_text(builder(r))
        n += 1
    print(f"Wrote {n} {fmt} file(s) to {out}/")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="onkos", description="Onkos — TGI/survival dataset toolkit.")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("version", help="print version").set_defaults(func=_cmd_version)
    sub.add_parser("validate", help="JSON-Schema-validate the dataset").set_defaults(func=_cmd_validate)
    sub.add_parser("info", help="counts by subsystem / tier / review status").set_defaults(func=_cmd_info)

    rp = sub.add_parser("report", help="dataset health & validation report (Markdown)")
    rp.add_argument("--output", default=None, help="write report to a file instead of stdout")
    rp.set_defaults(func=_cmd_report)

    sp = sub.add_parser("simulate", help="population-level forward simulation")
    sp.add_argument("record", nargs="?", help="record id (omit with --compare)")
    sp.add_argument("--tumor-type", default="NSCLC")
    sp.add_argument("--line", default="first")
    sp.add_argument("--drug-effect", type=float, default=1.0)
    sp.add_argument(
        "--exposure", type=float, default=None,
        help="PK exposure metric (e.g. C_avg in ug/L); requires --exposure-response",
    )
    sp.add_argument(
        "--exposure-response", default=None,
        help="exposure-response record id mapping exposure -> drug effect E",
    )
    sp.add_argument("--compare", action="store_true", help="virtual-trial divergence across models")
    sp.set_defaults(func=_cmd_simulate)

    up = sub.add_parser(
        "uncertainty", help="Monte-Carlo parameter-uncertainty bands (propagates IIV CV)"
    )
    up.add_argument("record", help="record id")
    up.add_argument("--tumor-type", default="NSCLC")
    up.add_argument("--line", default="first")
    up.add_argument("--drug-effect", type=float, default=1.0)
    up.add_argument("--exposure", type=float, default=None)
    up.add_argument("--exposure-response", default=None)
    up.add_argument("--n", type=int, default=200, help="number of Monte-Carlo samples")
    up.add_argument("--seed", type=int, default=0)
    up.set_defaults(func=_cmd_uncertainty)

    np_ = sub.add_parser(
        "sensitivity", help="rank parameters by how much their IIV drives a target metric"
    )
    np_.add_argument("record", help="record id")
    np_.add_argument("--tumor-type", default="NSCLC")
    np_.add_argument("--line", default="first")
    np_.add_argument("--target", default="median_os_weeks", help="metric key or median_os_weeks")
    np_.add_argument("--n", type=int, default=400, help="number of Monte-Carlo samples")
    np_.add_argument("--seed", type=int, default=0)
    np_.set_defaults(func=_cmd_sensitivity)

    ep = sub.add_parser("export", help="generate export artifacts")
    ep.add_argument(
        "--format",
        required=True,
        choices=list(_TEXT_EXPORTERS) + ["omex", "csv", "bibtex"],
    )
    ep.add_argument("--output", required=True)
    ep.set_defaults(func=_cmd_export)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
