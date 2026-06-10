"""Onkos command-line interface."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

from . import __version__, compare, load, sensitivity, simulate, simulate_ensemble
from .combine import SCHEMES, WEIGHTS_ARE_COMBINATION_NOT_POSTERIOR
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


def _cmd_audit(_args) -> int:
    from .audit import audit_tiers

    findings = audit_tiers(load())
    inflated = [f for f in findings if f.status == "inflated"]
    conservative = [f for f in findings if f.status == "conservative"]
    print(f"Evidence-based tier audit — {len(findings)} clinical TGI/survival records\n")
    print(f"  {'record':<44} {'tier':>4} {'ceiling':>8} {'status':>13}")
    for f in findings:
        mark = "!" if f.status == "inflated" else " "
        print(f"{mark} {f.record_id:<44} {f.assigned:>4} {f.ceiling:>8} {f.status:>13}")
    print(f"\n  inflated (tier exceeds evidence): {len(inflated)}")
    print(f"  conservative (could upgrade if evidence trusted): {len(conservative)}")
    if inflated:
        print("\n  FAIL: tier inflation detected.", file=sys.stderr)
        return 1
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
        if args.json:
            print(cmp.to_json(include_curves=args.include_curves))
            return 0
        print(f"Virtual-trial comparison — {ctx}, {effect_desc}\n")
        for tr in cmp.included:
            mos = f"{tr.median_os:.1f}" if tr.median_os else "n/r"
            mpfs = f"{tr.median_pfs:.1f}" if tr.median_pfs else "n/r"
            print(f"  [{tr.tier}] {tr.record_id:<48} median OS {mos:>6}  PFS {mpfs:>6}")
        for rid, reason in cmp.excluded:
            print(f"  [-] {rid:<48} EXCLUDED ({reason})")
        print(f"\n  OS  divergence {cmp.os_divergence:.3f}  | median OS range  {cmp.median_os_range}")
        print(f"  PFS divergence {cmp.pfs_divergence:.3f}  | median PFS range {cmp.median_pfs_range}")
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
    if tr.median_pfs is not None:
        print(f"  median_pfs_weeks         {tr.median_pfs:.2f}")
    for w in tr.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_compare(args) -> int:
    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    er, exposure = args.exposure_response, args.exposure
    drug_effect = None if (er and exposure is not None) else args.drug_effect
    cmp = compare(
        ds, purpose="tgi", context=ctx, drug_effect=drug_effect,
        exposure=exposure, exposure_response=er,
    )
    if not cmp.included:
        print("No eligible (in-context) models for this context.", file=sys.stderr)
        return 1

    ma = None
    if args.average:
        ma = cmp.model_average(
            target=args.target, endpoint=args.endpoint, weights=args.weights,
            n=args.n, seed=args.seed,
        )

    if args.json:
        print(cmp.to_json(include_curves=args.include_curves, model_average=ma))
        return 0

    print(f"Virtual-trial comparison — {ctx}, drug_effect={args.drug_effect}\n")
    for tr in cmp.included:
        mos = f"{tr.median_os:.1f}" if tr.median_os else "n/r"
        print(f"  [{tr.tier}] {tr.record_id:<48} median OS {mos:>6}")
    for rid, reason in cmp.excluded:
        print(f"  [-] {rid:<48} EXCLUDED ({reason})")
    print(f"\n  OS divergence {cmp.os_divergence:.3f}  | median OS range {cmp.median_os_range}")

    if ma is not None:
        print(
            f"\n  Model-averaged {ma.target} = {ma.point:.1f}  [{ma.endpoint}, scheme={ma.scheme}, "
            f"tier={ma.tier}]"
        )
        print(
            f"  Variance: within(parameter)={ma.within_var:.1f}  "
            f"between(model-selection)={ma.between_var:.1f}"
        )
        print(
            f"  >> model_selection_fraction = {ma.model_selection_fraction:.2f}  "
            f"(irreducible model-choice risk)"
        )
        print(f"  weights: {WEIGHTS_ARE_COMBINATION_NOT_POSTERIOR}")
        for rid, w in ma.weights.items():
            print(f"    {w:>6.3f}  {rid}")
        if args.decompose:
            dec = cmp.uncertainty_decomposition(target=args.target, endpoint=args.endpoint,
                                                n=args.n, seed=args.seed)
            print("\n  Per-scheme decomposition:")
            print(f"    {'scheme':<10} {'point':>8} {'within':>9} {'between':>9} {'frac':>6}")
            for scheme, row in dec.items():
                print(f"    {scheme:<10} {row['point']:>8.1f} {row['within_var']:>9.1f} "
                      f"{row['between_var']:>9.1f} {row['model_selection_fraction']:>6.2f}")
            print(f"\n  weight_sensitivity (point swing across schemes) = {ma.weight_sensitivity:.2f}")
        for w in ma.warnings:
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


def _cmd_identify(args) -> int:
    from .identify import identifiability

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    schedule = [float(x) for x in args.schedule.split(",")] if args.schedule else None
    kw = {"schedule": schedule} if schedule else {}
    res = identifiability(
        ds, args.record, context=ctx, drug_effect=args.drug_effect,
        sigma_prop=args.sigma_prop, sigma_add=args.sigma_add, **kw,
    )
    if args.json:
        print(res.to_json())
        return 0
    sched = ", ".join(f"{x:g}" for x in res.schedule)
    print(
        f"{args.record}  tier={res.tier}  (design: scans [{sched}] wk, "
        f"residual σ={res.sigma_prop:g}·y + {res.sigma_add:g})\n"
    )
    print(f"  {'parameter':<12} {'central':>10} {'pred. RSE':>10} {'IIV CV':>8}  identifiable?")
    for p in res.params:
        rse = "  inf" if not math.isfinite(p.rse_percent) else f"{p.rse_percent:>8.0f}%"
        cv = "   -" if p.iiv_cv_percent is None else f"{p.iiv_cv_percent:>6.0f}%"
        mark = "yes" if p.identifiable else "NO"
        print(f"  {p.symbol:<12} {p.central:>10.4g} {rse:>10} {cv:>8}  {mark}")
    gamma = "inf" if not math.isfinite(res.collinearity_index) else f"{res.collinearity_index:.1f}"
    verdict = "IDENTIFIABLE" if res.practically_identifiable else "NOT identifiable"
    print(
        f"\n  collinearity index γ_K = {gamma}  (ceiling {res.collinearity_ceiling:g})  "
        f"->  {verdict} under this design"
    )
    if res.worst and not res.worst.identifiable:
        print(f"  -> richer design / external constraint needed first for '{res.worst.symbol}'.")
    for w in res.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_design(args) -> int:
    from .design import optimal_schedule

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    od = optimal_schedule(
        ds, args.record, context=ctx, drug_effect=args.drug_effect,
        n_samples=args.n_samples, horizon=args.horizon,
        sigma_prop=args.sigma_prop, sigma_add=args.sigma_add,
    )
    if args.json:
        print(od.to_json())
        return 0
    print(
        f"{args.record}  tier={od.tier}  (D-optimal design: {od.n_samples} scans over "
        f"{od.horizon:g} wk, residual σ={args.sigma_prop:g}·y)\n"
    )
    syms = list(od.uniform.rse_percent)
    print(f"  {'parameter':<12} {'uniform RSE':>12} {'D-opt RSE':>12}  identifiable (D-opt)?")
    for s in syms:
        u = od.uniform.rse_percent[s]
        o = od.optimal.rse_percent[s]
        us = "  inf" if not math.isfinite(u) else f"{u:>10.0f}%"
        os_ = "  inf" if not math.isfinite(o) else f"{o:>10.0f}%"
        mark = "yes" if s in od.optimal.identifiable else "NO"
        print(f"  {s:<12} {us:>12} {os_:>12}  {mark}")
    sched = ", ".join(f"{x:g}" for x in od.optimal.schedule)
    print(f"\n  D-optimal schedule (wk) = [{sched}]")
    print(f"  D-efficiency over uniform = {od.d_efficiency:.2f}x  "
          f"(collinearity {od.uniform.collinearity_index:.1f} -> {od.optimal.collinearity_index:.1f})")
    if od.rescues_any:
        rescued = sorted(set(od.optimal.identifiable) - set(od.uniform.identifiable))
        print(f"  >> the better design RESCUED: {', '.join(rescued)} (crossed the identifiability line)")
    if od.structurally_flat:
        print(f"  >> STRUCTURALLY FLAT even under the best design: {', '.join(od.structurally_flat)}")
    for w in od.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_response(args) -> int:
    from .response import objective_response_rate, response_vs_survival

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    link = args.survival_link

    if args.surrogate:
        rs = response_vs_survival(ds, context=ctx, survival_link=link, n=args.n, seed=args.seed)
        if args.json:
            print(rs.to_json())
            return 0
        link_name = link.split(".")[-1] if link else "week-8 (default)"
        print(f"ORR -> OS surrogate — {ctx}, OS link = {link_name}\n")
        print(f"  {'model':<48} {'ORR':>5} {'DCR':>5} {'median OS':>10}")
        for r in sorted(rs.rows, key=lambda x: -x["orr"]):
            mos = f"{r['median_os_weeks']:.0f}" if r["median_os_weeks"] else "n/r"
            print(f"  [{r['tier']}] {r['record_id']:<44} {r['orr']:>5.2f} {r['dcr']:>5.2f} {mos:>10}")
        verdict = "ORR faithfully ranks OS" if rs.orr_predicts_os else "ORR MIS-RANKS OS"
        print(
            f"\n  >> discordant model pairs = {rs.discordant_pairs}/{rs.total_pairs} "
            f"(fraction {rs.discordant_fraction:.2f})  ->  {verdict}"
        )
        return 0

    if args.durability:
        rs = response_vs_survival(ds, context=ctx, survival_link=link, n=args.n, seed=args.seed)
        if args.json:
            print(rs.to_json())
            return 0
        link_name = link.split(".")[-1] if link else "week-8 (default)"
        print(f"Breadth vs durability — {ctx}, OS link = {link_name}\n")
        print(f"  {'model':<48} {'ORR':>5} {'DoR':>7} {'OS':>6}   (ORR=breadth, DoR=durability)")
        for r in sorted(rs.rows, key=lambda x: -x["orr"]):
            dor = "n/r" if r["median_dor_weeks"] is None else f"{r['median_dor_weeks']:.0f}"
            mos = f"{r['median_os_weeks']:.0f}" if r["median_os_weeks"] else "n/r"
            print(f"  [{r['tier']}] {r['record_id']:<44} {r['orr']:>5.2f} {dor:>7} {mos:>6}")
        print("\n  >> depth is not durability: the highest response rate can be the least durable")
        return 0

    if not args.record:
        print("error: provide a RECORD id (or --surrogate / --durability)", file=sys.stderr)
        return 2
    rr = objective_response_rate(ds, args.record, context=ctx, drug_effect=args.drug_effect,
                                 survival_link=link, n=args.n, seed=args.seed)
    if args.json:
        print(rr.to_json())
        return 0
    print(f"{args.record}  tier={rr.tier}  ({ctx}, n={rr.n})\n")
    print(f"  ORR (objective response rate) = {rr.orr:.2f}    DCR (disease control) = {rr.dcr:.2f}")
    if rr.median_dor_weeks is not None:
        cens = f"  ({rr.dor_censored_fraction:.0%} censored)" if rr.dor_censored_fraction else ""
        print(f"  median DoR (durability)       = {rr.median_dor_weeks:.0f} wk{cens}  "
              f"[{rr.n_responders} responders]")
    if rr.median_os_weeks is not None:
        print(f"  median OS (same trial)        = {rr.median_os_weeks:.0f} wk")
    print("\n  RECIST best-response distribution:")
    for cat in ("CR", "PR", "SD", "PD"):
        frac = rr.distribution[cat]
        print(f"    {cat}  {frac * 100:>4.0f}%  {'█' * round(frac * 30)}")
    for w in rr.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_pfs(args) -> int:
    from .response import pfs_route_divergence, progression_free_survival

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}

    if args.routes:
        div = pfs_route_divergence(ds, context=ctx, drug_effect=args.drug_effect,
                                   landmark_weeks=args.landmark, n=args.n, seed=args.seed)
        if args.json:
            print(div.to_json())
            return 0
        print(f"PFS two routes — {ctx}  (mechanistic RECIST TTP vs statistical week-8 link)\n")
        print(f"  {'model':<48} {'mech TTP':>9} {'stat PFS':>9} {'ratio':>6}")
        for r in sorted(div.rows, key=lambda x: -(x["median_ttp_weeks"] or 0)):
            mech = "n/r" if r["median_ttp_weeks"] is None else f"{r['median_ttp_weeks']:.0f}"
            stat = "n/r" if r["median_pfs_link_weeks"] is None else f"{r['median_pfs_link_weeks']:.0f}"
            ratio = "—" if r["route_ratio"] is None else f"{r['route_ratio']:.2f}"
            print(f"  [{r['tier']}] {r['record_id']:<44} {mech:>9} {stat:>9} {ratio:>6}")
        verdict = "the routes agree" if div.routes_agree else "the PFS ROUTE inverts the ranking"
        print(
            f"\n  >> route-discordant model pairs = {div.discordant_pairs}/{div.total_pairs} "
            f"(fraction {div.discordant_fraction:.2f})  ->  {verdict}"
        )
        return 0

    if not args.record:
        print("error: provide a RECORD id (or --routes)", file=sys.stderr)
        return 2
    pf = progression_free_survival(ds, args.record, context=ctx, drug_effect=args.drug_effect,
                                   landmark_weeks=args.landmark, n=args.n, seed=args.seed)
    if args.json:
        print(pf.to_json())
        return 0
    print(f"{args.record}  tier={pf.tier}  ({ctx}, n={pf.n})\n")
    mech = "n/r" if pf.median_ttp_weeks is None else f"{pf.median_ttp_weeks:.0f} wk"
    cens = f"  ({pf.ttp_censored_fraction:.0%} censored)" if pf.ttp_censored_fraction else ""
    print(f"  mechanistic PFS (RECIST TTP)  = {mech}{cens}")
    print(f"  progression-free @ {pf.landmark_weeks:.0f} wk    = {pf.mechanistic_pfs_rate:.2f}")
    if pf.has_pfs_link and pf.median_pfs_link_weeks is not None:
        print(f"  statistical PFS (week-8 link) = {pf.median_pfs_link_weeks:.0f} wk")
    if pf.route_ratio is not None:
        agree = "routes agree" if 0.85 <= pf.route_ratio <= 1.18 else "routes disagree"
        print(f"  route ratio (mech / stat)     = {pf.route_ratio:.2f}  ({agree})")
    for w in pf.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_joint(args) -> int:
    from .joint import compare_joint_vs_two_stage

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    cmp = compare_joint_vs_two_stage(ds, context=ctx, drug_effect=args.drug_effect, alpha=args.alpha)
    print(
        f"Two-stage (week-8, proportional hazards) vs joint (current-value link, alpha={args.alpha}) "
        f"— {ctx}\n"
    )
    print(f"  {'model':<48} {'2-stage':>8} {'joint':>7} {'HR(end)/HR(8wk)':>16}")
    for r in sorted(cmp.rows, key=lambda x: -(x["joint_median"] or 0)):
        ts = "n/r" if r["two_stage_median"] is None else f"{r['two_stage_median']:.0f}"
        jt = "n/r" if r["joint_median"] is None else f"{r['joint_median']:.0f}"
        phv = "—" if not _finite(r["ph_violation"]) else f"{r['ph_violation']:.1f}x"
        print(f"  [{r['tier']}] {r['record_id']:<44} {ts:>8} {jt:>7} {phv:>16}")
    print(
        f"\n  >> rank-discordant model pairs (two-stage vs joint) = {cmp.rank_discordant_pairs}"
        f"   max PH-violation = {cmp.max_ph_violation:.0f}x"
    )
    print(
        "  the two-stage links assume a CONSTANT hazard ratio; the joint link's HR rises as a "
        "resistant clone regrows — a non-proportional hazard they cannot represent."
    )
    return 0


def _finite(x) -> bool:
    import math

    return isinstance(x, (int, float)) and math.isfinite(x)


def _cmd_budget(args) -> int:
    from .budget import _AXIS_LABELS, model_selection_budget

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    b = model_selection_budget(ds, context=ctx, endpoint=args.endpoint, n=args.n, seed=args.seed)
    if args.json:
        print(b.to_json())
        return 0
    print(
        f"Model-selection budget — {ctx}, {args.endpoint}  "
        f"(grid: {len(b.models)} TGI models x {len(b.links)} survival links, tier {b.tier})\n"
    )
    print(f"  grand-mean {b.target} = {b.grand_mean:.1f}   total variance = {b.total:.1f}\n")
    print(f"  {'axis':<46} {'share':>7}")
    for axis, frac in sorted(b.fractions.items(), key=lambda kv: kv[1], reverse=True):
        bar = "█" * round(frac * 30)
        print(f"  {_AXIS_LABELS[axis]:<46} {frac * 100:>5.0f}%  {bar}")
    print(
        f"\n  >> structural-choice share = {b.structural_fraction * 100:.0f}% "
        f"(irreducible by a bigger trial); parameter share = "
        f"{b.fractions['parameter'] * 100:.0f}%"
    )
    print(f"  >> dominant axis: {_AXIS_LABELS[b.dominant]} — standardize / validate this first")
    for w in b.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_interactions(args) -> int:
    from .interaction import SYNERGY_IS_AN_ASSUMPTION, compare_interactions

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    cmp = compare_interactions(
        ds, args.record, context=ctx, effect_a=args.effect_a, effect_b=args.effect_b,
        psi=args.psi,
    )
    if args.json:
        print(cmp.to_json())
        return 0
    print(
        f"{args.record}  tier={cmp.tier}  (combination: E_A={args.effect_a:g} + "
        f"E_B={args.effect_b:g}, {ctx})\n"
    )
    print(f"  {'interaction':<14} {'combined E':>11} {'median OS':>11}")
    for label, tr in cmp.trajectories.items():
        e = cmp.combined_effects[label]
        mos = f"{tr.median_os:.1f}" if tr.median_os else "n/r"
        print(f"  {label:<14} {e:>11.3f} {mos:>11}")
    rng = cmp.median_os_range
    span = f"{rng[0]:.0f}-{rng[1]:.0f}" if rng else "n/r"
    print(
        f"\n  OS divergence across interaction models = {cmp.os_divergence:.3f}  "
        f"| median OS range {span} wk"
    )
    print("  >> the interaction model is a model-selection axis (not a measured quantity)")
    if cmp.warnings:
        print(f"  ! {SYNERGY_IS_AN_ASSUMPTION}")
    return 0


def _cmd_atlas(args) -> int:
    from .atlas import model_selection_atlas

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    a = model_selection_atlas(ds, context=ctx)
    if args.json:
        print(a.to_json())
        return 0
    print(
        f"Model-selection atlas — {ctx}  tier={a.tier}\n"
        f"  each axis in its OWN unit (a survey, not a decomposition; see `onkos budget`)\n"
    )
    print(f"  {'axis':<26} {'headline':>10}  unit / detail")
    for e in a.entries:
        h = "n/a" if e.headline is None else f"{e.headline}"
        print(f"  {e.label:<26} {h:>10}  {e.unit}")
        print(f"  {'':<26} {'':>10}  └ {e.detail}")
    print(
        "\n  >> a one-call map of where the model-selection risk lies for this context; "
        "use each axis's own command for the deep dive (onkos joint / dose-response / "
        "discriminability / …)."
    )
    return 0


def _cmd_discriminability(args) -> int:
    from .discriminability import model_discriminability

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    md = model_discriminability(ds, context=ctx, survival_link=args.survival_link,
                                power=args.power, alpha=args.alpha)
    if args.json:
        print(md.to_json())
        return 0
    link = args.survival_link.split(".")[-1] if args.survival_link else "default week-8"
    print(
        f"tier={md.tier}  ({ctx}, link={link}, power={args.power}, alpha={args.alpha})\n"
        f"  required events to distinguish each model pair's OS curves:\n"
    )
    for p in sorted(md.pairs, key=lambda x: x["required_events"]):
        a = p["record_a"].split(".")[-1]
        b = p["record_b"].split(".")[-1]
        d = p["required_events"]
        if not (d == d) or d == float("inf"):  # nan/inf
            tag, ds_ = "INDISTINGUISHABLE", "inf"
        else:
            tag = "feasible" if d < 500 else ("large" if d < 3000 else "INFEASIBLE")
            ds_ = f"{d:,.0f}"
        print(f"  HR={p['hazard_ratio']:>5.2f}  events={ds_:>10}  {tag:<17} {a} vs {b}")
    print(
        f"\n  >> {md.n_indistinguishable}/{len(md.pairs)} model pairs are practically "
        "indistinguishable (need an infeasible trial): the model choice cannot be resolved by "
        "the data — only assumed. The silent model-selection risk, quantified in events."
    )
    return 0


def _cmd_early_surrogate(args) -> int:
    from .early_surrogate import surrogate_timing_fidelity

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    st = surrogate_timing_fidelity(ds, context=ctx, reference_link=args.reference_link)
    if args.json:
        print(st.to_json())
        return 0
    ref = [r.split(".")[-1] for r in st.reference_ranking]
    print(
        f"tier={st.tier}  ({ctx})  durable-benefit reference: {st.reference_link.split('.')[-1]}\n"
        f"  tail-aware ranking (best->worst): {' > '.join(ref)}\n"
    )
    print(f"  {'landmark wk':>11} | {'discordant pairs vs durable benefit':>36}")
    for r in st.rows:
        bar = "#" * r["discordant_pairs"]
        print(f"  {r['week']:>11.0f} | {r['discordant_pairs']:>2}/{st.total_pairs}  {bar}")
    print(
        f"\n  >> earliest readout = {st.earliest_discordance}/{st.total_pairs} discordant, "
        f"latest = {st.latest_discordance}/{st.total_pairs}: reading the surrogate earlier "
        "over-rewards deep-but-doomed responders (the ctDNA-timing bias)"
    )
    return 0


def _cmd_dose_response(args) -> int:
    from .dose_response import compare_er_extrapolation

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    cmp = compare_er_extrapolation(ds, args.record, context=ctx, c_ref=args.c_ref, e_ref=args.e_ref)
    if args.json:
        print(cmp.to_json())
        return 0
    print(
        f"{args.record}  tier={cmp.tier}  (ER-shape extrapolation anchored at "
        f"C_ref={args.c_ref:g} -> E_ref={args.e_ref:g}, {ctx})\n"
    )
    shapes = [er.split(".")[-1] for er in cmp.er_ids]
    print(f"  {'dose':>7} | " + " ".join(f"{s:>7}" for s in shapes) + f" | {'OS spread':>9}")
    for r in cmp.rows:
        effs = " ".join(f"{r['effects'][er]:>7.2f}" for er in cmp.er_ids)
        tag = "  <- anchor" if abs(r["dose"] - args.c_ref) < 1e-9 else ""
        print(f"  {r['dose']:>7.0f} | {effs} | {r['os_divergence']:>8.0f}w{tag}")
    print(
        f"\n  OS spread at the studied dose = {cmp.reference_os_divergence:.1f}w (anchored); "
        f"max on extrapolation = {cmp.max_os_divergence:.0f}w"
    )
    print(
        "  >> the ER-model choice is invisible at the studied dose but a model-selection axis "
        "off it — a dose-extrapolation risk, sharpest on de-escalation"
    )
    return 0


def _cmd_loewe(args) -> int:
    from .interaction import compare_additivity_references

    ds = load()
    ctx = {"tumor_type": args.tumor_type, "line": args.line}
    cmp = compare_additivity_references(
        ds, args.record, context=ctx, dose_a=args.dose_a, dose_b=args.dose_b,
        er_a=args.er_a, er_b=args.er_b,
    )
    if args.json:
        print(cmp.to_json())
        return 0
    print(
        f"{args.record}  tier={cmp.tier}  (doses d_A={args.dose_a:g} via {args.er_a}, "
        f"d_B={args.dose_b:g} via {args.er_b}, {ctx})\n"
    )
    print(f"  {'reference':<8} {'combined E':>11} {'median OS':>11}")
    for ref, tr in cmp.trajectories.items():
        mos = f"{tr.median_os:.1f}" if tr.median_os else "n/r"
        print(f"  {ref:<8} {cmp.combined_effects[ref]:>11.3f} {mos:>11}")
    rng = cmp.median_os_range
    span = f"{rng[0]:.0f}-{rng[1]:.0f}" if rng else "n/r"
    print(
        f"\n  OS divergence across additivity references = {cmp.os_divergence:.3f}  "
        f"| median OS range {span} wk"
    )
    print(
        "  >> the additivity REFERENCE is a model-selection axis; Loewe alone passes the "
        "sham-combination test (a drug with itself is exactly additive)"
    )
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
    sub.add_parser("audit", help="evidence-based tier audit (flags tier inflation)").set_defaults(func=_cmd_audit)

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
    sp.add_argument("--json", action="store_true", help="emit the comparison as JSON (with --compare)")
    sp.add_argument("--include-curves", action="store_true", help="include tumor/OS/PFS arrays in --json")
    sp.set_defaults(func=_cmd_simulate)

    cp = sub.add_parser(
        "compare",
        help="virtual-trial divergence + model averaging (model-selection uncertainty)",
    )
    cp.add_argument("--tumor-type", default="NSCLC")
    cp.add_argument("--line", default="first")
    cp.add_argument("--drug-effect", type=float, default=1.0)
    cp.add_argument("--exposure", type=float, default=None)
    cp.add_argument("--exposure-response", default=None)
    cp.add_argument("--average", action="store_true",
                    help="model-average the eligible models + decompose the uncertainty")
    cp.add_argument("--weights", default="equal", choices=list(SCHEMES),
                    help="weighting scheme (combination weights, NOT model posteriors)")
    cp.add_argument("--target", default="median_os_weeks",
                    help="target: median_os_weeks | median_pfs_weeks | a metric key")
    cp.add_argument("--endpoint", default="OS", choices=["OS", "PFS"])
    cp.add_argument("--decompose", action="store_true", help="show the per-scheme decomposition")
    cp.add_argument("--n", type=int, default=200, help="within-model ensemble depth")
    cp.add_argument("--seed", type=int, default=0)
    cp.add_argument("--json", action="store_true", help="emit the result as JSON")
    cp.add_argument("--include-curves", action="store_true", help="include curve arrays in --json")
    cp.set_defaults(func=_cmd_compare)

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

    ip = sub.add_parser(
        "identify",
        help="practical identifiability: can a trial design even estimate the parameters?",
    )
    ip.add_argument("record", help="record id (a dynamic TGI/growth model)")
    ip.add_argument("--tumor-type", default="NSCLC")
    ip.add_argument("--line", default="first")
    ip.add_argument("--drug-effect", type=float, default=1.0)
    ip.add_argument("--schedule", default=None,
                    help="comma-separated scan times in weeks (default: 0,6,12,18,24,36,48)")
    ip.add_argument("--sigma-prop", type=float, default=0.2,
                    help="proportional residual error (CV) of the tumor-size assay")
    ip.add_argument("--sigma-add", type=float, default=0.0, help="additive residual error")
    ip.add_argument("--json", action="store_true", help="emit the result as JSON")
    ip.set_defaults(func=_cmd_identify)

    dp = sub.add_parser(
        "design",
        help="D-optimal trial design: the best sampling schedule a fixed budget allows",
    )
    dp.add_argument("record", help="record id (a dynamic TGI/growth model)")
    dp.add_argument("--tumor-type", default="NSCLC")
    dp.add_argument("--line", default="first")
    dp.add_argument("--drug-effect", type=float, default=1.0)
    dp.add_argument("--n-samples", type=int, default=7, help="measurement budget (# scans)")
    dp.add_argument("--horizon", type=float, default=48.0, help="design horizon in weeks")
    dp.add_argument("--sigma-prop", type=float, default=0.2,
                    help="proportional residual error (CV) of the tumor-size assay")
    dp.add_argument("--sigma-add", type=float, default=0.0, help="additive residual error")
    dp.add_argument("--json", action="store_true", help="emit the result as JSON")
    dp.set_defaults(func=_cmd_design)

    rsp = sub.add_parser(
        "response",
        help="RECIST best response / ORR, and the contested ORR -> OS surrogate",
    )
    rsp.add_argument("record", nargs="?", help="record id (omit with --surrogate)")
    rsp.add_argument("--tumor-type", default="NSCLC")
    rsp.add_argument("--line", default="first")
    rsp.add_argument("--drug-effect", type=float, default=1.0)
    rsp.add_argument("--survival-link", default=None,
                     help="non-default OS link (e.g. survival_link.nsclc_os_growth_rate)")
    rsp.add_argument("--surrogate", action="store_true",
                     help="ORR -> OS discordance across the in-context TGI models")
    rsp.add_argument("--durability", action="store_true",
                     help="breadth (ORR) vs durability (median DoR) table across the models")
    rsp.add_argument("--n", type=int, default=300, help="ensemble depth")
    rsp.add_argument("--seed", type=int, default=0)
    rsp.add_argument("--json", action="store_true", help="emit the result as JSON")
    rsp.set_defaults(func=_cmd_response)

    pp = sub.add_parser(
        "pfs",
        help="progression-free survival two ways: mechanistic RECIST TTP vs the statistical link",
    )
    pp.add_argument("record", nargs="?", help="record id (omit with --routes)")
    pp.add_argument("--tumor-type", default="NSCLC")
    pp.add_argument("--line", default="first")
    pp.add_argument("--drug-effect", type=float, default=1.0)
    pp.add_argument("--landmark", type=float, default=24.0,
                    help="landmark horizon (weeks) for the progression-free rate")
    pp.add_argument("--routes", action="store_true",
                    help="two-route PFS table + route-discordance across the in-context models")
    pp.add_argument("--n", type=int, default=300, help="ensemble depth")
    pp.add_argument("--seed", type=int, default=0)
    pp.add_argument("--json", action="store_true", help="emit the result as JSON")
    pp.set_defaults(func=_cmd_pfs)

    jp = sub.add_parser(
        "joint",
        help="joint (current-value) vs two-stage survival: the non-proportional-hazard axis",
    )
    jp.add_argument("--tumor-type", default="NSCLC")
    jp.add_argument("--line", default="first")
    jp.add_argument("--drug-effect", type=float, default=1.0)
    jp.add_argument("--alpha", type=float, default=1.0,
                    help="association between log tumor size and log hazard (DECLARED, not fitted)")
    jp.set_defaults(func=_cmd_joint)

    bp = sub.add_parser(
        "budget",
        help="model-selection budget: variance split across the structural choices",
    )
    bp.add_argument("--tumor-type", default="NSCLC")
    bp.add_argument("--line", default="first")
    bp.add_argument("--endpoint", default="OS", choices=["OS", "PFS"])
    bp.add_argument("--n", type=int, default=200, help="within-cell ensemble depth")
    bp.add_argument("--seed", type=int, default=0)
    bp.add_argument("--json", action="store_true", help="emit the result as JSON")
    bp.set_defaults(func=_cmd_budget)

    xp = sub.add_parser(
        "interactions",
        help="drug-combination divergence: the interaction model as a model-selection axis",
    )
    xp.add_argument("record", help="record id (a TGI model driven by a drug effect)")
    xp.add_argument("--tumor-type", default="NSCLC")
    xp.add_argument("--line", default="first")
    xp.add_argument("--effect-a", type=float, default=0.6, help="single-agent effect of drug A")
    xp.add_argument("--effect-b", type=float, default=0.6, help="single-agent effect of drug B")
    xp.add_argument("--psi", type=float, default=0.5,
                    help="synergy/antagonism bracket magnitude (DECLARED assumption, not fitted)")
    xp.add_argument("--json", action="store_true", help="emit the result as JSON")
    xp.set_defaults(func=_cmd_interactions)

    lp = sub.add_parser(
        "loewe",
        help="dose-level additivity references (HSA / Bliss / Loewe) as a model-selection axis",
    )
    lp.add_argument("record", help="record id (a TGI model driven by a drug effect)")
    lp.add_argument("--tumor-type", default="NSCLC")
    lp.add_argument("--line", default="first")
    lp.add_argument("--er-a", default="exposure_response.emax_generic", help="ER record for drug A")
    lp.add_argument("--er-b", default="exposure_response.dacomitinib_egfr.emax",
                    help="ER record for drug B")
    lp.add_argument("--dose-a", type=float, default=150.0, help="dose of drug A (exposure units)")
    lp.add_argument("--dose-b", type=float, default=90.0, help="dose of drug B (exposure units)")
    lp.add_argument("--json", action="store_true", help="emit the result as JSON")
    lp.set_defaults(func=_cmd_loewe)

    drp = sub.add_parser(
        "dose-response",
        help="exposure-response model choice as a dose-extrapolation model-selection axis",
    )
    drp.add_argument("record", help="record id (a TGI model driven by a drug effect)")
    drp.add_argument("--tumor-type", default="NSCLC")
    drp.add_argument("--line", default="first")
    drp.add_argument("--c-ref", type=float, default=150.0, help="reference exposure (the studied dose)")
    drp.add_argument("--e-ref", type=float, default=1.0, help="effect the shapes share at C_ref")
    drp.add_argument("--json", action="store_true", help="emit the result as JSON")
    drp.set_defaults(func=_cmd_dose_response)

    esp = sub.add_parser(
        "early-surrogate",
        help="early-surrogate readout timing: how landmark week trades against durable-benefit fidelity",
    )
    esp.add_argument("--tumor-type", default="NSCLC")
    esp.add_argument("--line", default="first")
    esp.add_argument("--reference-link", default=None,
                     help="durable-benefit reference link (default: the context's k_g OS link)")
    esp.add_argument("--json", action="store_true", help="emit the result as JSON")
    esp.set_defaults(func=_cmd_early_surrogate)

    dcp = sub.add_parser(
        "discriminability",
        help="model discriminability: required trial events to tell competing models apart",
    )
    dcp.add_argument("--tumor-type", default="NSCLC")
    dcp.add_argument("--line", default="first")
    dcp.add_argument("--survival-link", default=None, help="survival link (default: week-8)")
    dcp.add_argument("--power", type=float, default=0.8)
    dcp.add_argument("--alpha", type=float, default=0.05)
    dcp.add_argument("--json", action="store_true", help="emit the result as JSON")
    dcp.set_defaults(func=_cmd_discriminability)

    atp = sub.add_parser(
        "atlas",
        help="model-selection atlas: a one-call survey of every axis's headline for a context",
    )
    atp.add_argument("--tumor-type", default="NSCLC")
    atp.add_argument("--line", default="first")
    atp.add_argument("--json", action="store_true", help="emit the result as JSON")
    atp.set_defaults(func=_cmd_atlas)

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
