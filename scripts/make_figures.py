#!/usr/bin/env python3
"""Generate the figures embedded in the README (deterministic, regenerable)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import onkos
from onkos.export.reference import effect
from onkos.export.registry import get_kernel, kernel_values

plt.switch_backend("Agg")

OUT = Path(__file__).resolve().parents[1] / "docs" / "images"
OUT.mkdir(parents=True, exist_ok=True)

PALETTE = ["#2b6cb0", "#c05621", "#2f855a", "#6b46c1", "#b7791f"]


def divergence_figure() -> None:
    ds = onkos.load()
    t = np.linspace(0.0, 104.0, 209)
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    cmp = onkos.compare(ds, purpose="tgi", context=ctx, drug_effect=1.0, t=t)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, tr in enumerate(cmp.included):
        c = PALETTE[i % len(PALETTE)]
        ax1.plot(t, tr.tumor_size, color=c, label=f"{tr.record_id.split('.')[1]} [{tr.tier}]")
        if tr.os_curve is not None:
            ax2.plot(t, tr.os_curve, color=c, label=tr.record_id.split(".")[1])

    ax1.set_title("Tumor-size trajectories (NSCLC, 1L, E=1.0)")
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("tumor size (mm, SLD)")
    ax1.legend(fontsize=8)

    ax2.axhline(0.5, ls=":", color="grey", lw=1)
    ax2.set_title(f"Population OS — divergence={cmp.os_divergence:.2f}")
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("survival fraction")
    ax2.set_ylim(0, 1.02)
    ax2.legend(fontsize=8)

    fig.suptitle(
        f"Virtual-trial divergence (NSCLC, 1L) — model choice moves median OS across "
        f"{cmp.median_os_range[0]:.0f}-{cmp.median_os_range[1]:.0f} wk; "
        f"{len(cmp.excluded)} models greyed out for out-of-context transport",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "divergence.png", dpi=120)
    plt.close(fig)


def tier_figure() -> None:
    ds = onkos.load()
    counts = ds.by_tier()
    tiers = ["A", "B", "C", "D"]
    vals = [counts.get(k, 0) for k in tiers]
    colors = ["#2f855a", "#3182ce", "#dd6b20", "#c53030"]
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(tiers, vals, color=colors)
    ax.set_title("Record confidence-tier distribution")
    ax.set_ylabel("records")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.05, str(v), ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "tiers.png", dpi=120)
    plt.close(fig)


def exposure_response_figure() -> None:
    """ER curves (E vs C) and a PK-driven tumor trajectory (Hypnos composability)."""
    ds = onkos.load()
    c = np.linspace(0, 600, 200)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    for i, rid in enumerate(
        [
            "exposure_response.emax_generic",
            "exposure_response.sigmoid_emax_generic",
            "exposure_response.power_generic",
            "exposure_response.dacomitinib_egfr.emax",
        ]
    ):
        r = ds[rid]
        e = effect(get_kernel(r), c, kernel_values(r))
        ax1.plot(c, e, color=PALETTE[i % len(PALETTE)], label=f"{rid.split('.', 1)[1]} [{r.tier}]")
    ax1.set_title("Exposure-response transforms: C → E")
    ax1.set_xlabel("exposure C (µg/L)")
    ax1.set_ylabel("drug effect E (effect-unit)")
    ax1.legend(fontsize=7)

    # PK-driven tumor dynamics: a constant vs a declining exposure profile.
    t = np.linspace(0.0, 104.0, 209)
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    er = "exposure_response.emax_generic"
    const = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                           exposure=300.0, exposure_response=er, t=t)
    decaying = 300.0 * np.exp(-0.02 * t)
    pk = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                        exposure=decaying, exposure_response=er, t=t)
    ax2.plot(t, const.tumor_size, color=PALETTE[0], label="constant C=300 µg/L")
    ax2.plot(t, pk.tumor_size, color=PALETTE[1], label="declining PK profile")
    ax2b = ax2.twinx()
    ax2b.plot(t, decaying, color=PALETTE[1], ls=":", lw=1, alpha=0.7)
    ax2b.set_ylabel("exposure C(t) (µg/L)", color=PALETTE[1])
    ax2.set_title("PK-driven tumor dynamics (Hypnos composability)")
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("tumor size (mm, SLD)")
    ax2.legend(fontsize=7, loc="upper left")

    fig.suptitle(
        "Exposure-response (Phase B): exposure metric drives the kill term; "
        "a declining PK profile lets resistance-like regrowth emerge",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "exposure_response.png", dpi=120)
    plt.close(fig)


def context_library_figure() -> None:
    """Cross-tumor-context divergence: per-context OS spread + model-choice risk."""
    ds = onkos.load()
    contexts = ["NSCLC", "breast", "CRC", "HCC", "melanoma"]
    t = np.linspace(0.0, 156.0, 313)  # 3 years, so long-OS contexts reach their median

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    divergences = []
    for i, tt in enumerate(contexts):
        cmp = onkos.compare(
            ds, purpose="tgi", context={"tumor_type": tt, "line": "first"}, drug_effect=1.0, t=t
        )
        color = PALETTE[i % len(PALETTE)]
        # Plot the OS envelope (min..max across eligible models) for this context.
        curves = np.vstack([tr.os_curve for tr in cmp.included])
        ax1.fill_between(t, curves.min(axis=0), curves.max(axis=0), color=color, alpha=0.25)
        ax1.plot(t, curves.mean(axis=0), color=color, lw=1.5, label=f"{tt} (n={len(cmp.included)})")
        divergences.append((tt, cmp.os_divergence))

    ax1.axhline(0.5, ls=":", color="grey", lw=1)
    ax1.set_title("Population OS by tumor context (band = model-choice spread)")
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("survival fraction")
    ax1.set_ylim(0, 1.02)
    ax1.legend(fontsize=8)

    labels = [d[0] for d in divergences]
    vals = [d[1] for d in divergences]
    ax2.bar(labels, vals, color=[PALETTE[i % len(PALETTE)] for i in range(len(labels))])
    ax2.set_title("Model-selection risk (max OS divergence) by context")
    ax2.set_ylabel("OS divergence (max pointwise)")
    for i, v in enumerate(vals):
        ax2.text(i, v + 0.004, f"{v:.2f}", ha="center", fontsize=8)

    fig.suptitle(
        "Tumor-context library (Phase C): each context has its own baseline + survival link "
        "+ >=2 eligible TGI models",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "context_library.png", dpi=120)
    plt.close(fig)


def preclinical_figure() -> None:
    """Simeoni preclinical model: exp->linear growth, transit-chain delayed death."""
    ds = onkos.load()
    rid = "preclinical_translation.simeoni_2004.xenograft"
    spec = get_kernel(rid_record := ds[rid])
    t = np.linspace(0.0, 45.0, 451)  # days
    ctx = {"tumor_type": "ovarian_xenograft"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: unperturbed vs treated total tumor weight (semilog).
    untreated = onkos.simulate(ds, rid, context=ctx, drug_effect=0.0, t=t)
    treated = onkos.simulate(ds, rid, context=ctx, drug_effect=120.0, t=t)
    ax1.semilogy(t, untreated.tumor_size, color=PALETTE[2], label="untreated (exp→linear growth)")
    ax1.semilogy(t, treated.tumor_size, color=PALETTE[0], label="treated (E=120 ng/mL)")
    ax1.set_title("Simeoni xenograft TGI — total tumor weight")
    ax1.set_xlabel("days")
    ax1.set_ylabel("tumor weight w (g, log)")
    ax1.legend(fontsize=8)

    # Right: the transit chain x1..x4 under treatment (delayed cell death).
    vals = kernel_values(rid_record)
    vals["w0"] = 0.25
    vals["E"] = 120.0
    from onkos.export.reference import init_vector
    from scipy.integrate import solve_ivp
    sol = solve_ivp(
        lambda tt, yy: spec.rhs(tt, yy, vals), (t[0], t[-1]),
        init_vector(spec, vals), t_eval=t, rtol=1e-8, atol=1e-10, method="LSODA",
    )
    labels = ["x1 proliferating", "x2 transit", "x3 transit", "x4 transit"]
    for i in range(4):
        ax2.plot(t, sol.y[i], color=PALETTE[i % len(PALETTE)], label=labels[i])
    ax2.plot(t, sol.y.sum(axis=0), color="black", lw=1.2, ls="--", label="w = Σx (observed)")
    ax2.set_title("Signal-distribution transit chain (delayed cell death)")
    ax2.set_xlabel("days")
    ax2.set_ylabel("compartment weight (g)")
    ax2.legend(fontsize=7)

    fig.suptitle(
        "Preclinical translation (Phase D): Simeoni 2004 — damaged cells traverse "
        "x1→x2→x3→x4 before dying, so kill is delayed",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "preclinical.png", dpi=120)
    plt.close(fig)


def immuno_oncology_figure() -> None:
    """Hypothesis-tier tumor-immune QSP: control vs escape, and checkpoint rescue."""
    ds = onkos.load()
    t = np.linspace(0.0, 200.0, 801)

    def sim(rid, E):
        return onkos.simulate(ds, rid, context={"tumor_type": "x", "y0": 10.0},
                              drug_effect=E, t=t).tumor_size

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    ax1.plot(t, sim("immuno_oncology.kuznetsov_1994.tumor_immune", 0.0),
             color=PALETTE[2], label="immunogenic — immune control")
    ax1.plot(t, sim("immuno_oncology.poorly_immunogenic.hypothesis", 0.0),
             color=PALETTE[1], label="poorly immunogenic — escape")
    ax1.axhline(500, ls=":", color="grey", lw=1)
    ax1.text(5, 510, "carrying capacity (1/β)", fontsize=7, color="grey")
    ax1.set_title("Tumor-immune dynamics: control vs escape")
    ax1.set_xlabel("time (nondimensional)")
    ax1.set_ylabel("tumor (nondimensional)")
    ax1.legend(fontsize=8)

    poorly = "immuno_oncology.poorly_immunogenic.hypothesis"
    for E in [0.0, 6.0, 12.0]:
        ax2.plot(t, sim(poorly, E), label=f"checkpoint effect E={E:g}")
    ax2.set_title("Checkpoint blockade rescues escape (bistable threshold)")
    ax2.set_xlabel("time (nondimensional)")
    ax2.set_ylabel("tumor (nondimensional)")
    ax2.legend(fontsize=8)

    fig.suptitle(
        "Immuno-oncology (Phase E): HYPOTHESIS-TIER (tier D) — qualitative shapes "
        "only, NOT FOR PREDICTION",
        fontsize=10, color="#c53030",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "immuno_oncology.png", dpi=120)
    plt.close(fig)


def coverage_figure() -> None:
    """Tier x subsystem coverage heatmap + external-validation gauge (Phase F)."""
    from onkos.report import external_validation_coverage

    ds = onkos.load()
    subsystems = sorted({r.subsystem for r in ds})
    tiers = ["A", "B", "C", "D"]
    grid = np.zeros((len(subsystems), len(tiers)), dtype=int)
    for r in ds:
        grid[subsystems.index(r.subsystem), tiers.index(r.tier)] += 1

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [3, 1]})

    im = ax1.imshow(grid, cmap="Blues", aspect="auto")
    ax1.set_xticks(range(len(tiers)), tiers)
    ax1.set_yticks(range(len(subsystems)), subsystems)
    ax1.set_xlabel("confidence tier")
    ax1.set_title("Records by subsystem × tier")
    for i in range(len(subsystems)):
        for j in range(len(tiers)):
            if grid[i, j]:
                ax1.text(j, i, str(grid[i, j]), ha="center", va="center",
                         color="white" if grid[i, j] > grid.max() / 2 else "black", fontsize=9)
    fig.colorbar(im, ax=ax1, shrink=0.8, label="records")

    n_val, n_elig, frac = external_validation_coverage(ds)
    ax2.barh([0], [1.0], color="#e2e8f0")
    ax2.barh([0], [frac], color="#2f855a")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(-1, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("fraction")
    ax2.set_title("External-validation coverage")
    ax2.text(0.5, 0, f"{n_val}/{n_elig}\n({frac * 100:.0f}%)", ha="center", va="center",
             fontsize=12, fontweight="bold")

    fig.suptitle("Dataset health (Phase F): tier/subsystem coverage + external-validation",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "coverage.png", dpi=120)
    plt.close(fig)


def uncertainty_figure() -> None:
    """Monte-Carlo parameter-uncertainty bands from the stored IIV CVs."""
    ds = onkos.load()
    t = np.linspace(0.0, 104.0, 209)
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    ens = onkos.simulate_ensemble(
        ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0, t=t, n=400, seed=0
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    b = ens.tumor_size
    ax1.fill_between(t, b.lo, b.hi, color=PALETTE[0], alpha=0.25, label=f"{ens.ci[0]:g}–{ens.ci[1]:g}%")
    ax1.plot(t, b.median, color=PALETTE[0], lw=1.6, label="median")
    ax1.set_title("Tumor size under parameter uncertainty (IIV)")
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("tumor size (mm, SLD)")
    ax1.legend(fontsize=8)

    o = ens.os_curve
    ax2.fill_between(t, o.lo, o.hi, color=PALETTE[2], alpha=0.25, label=f"{ens.ci[0]:g}–{ens.ci[1]:g}%")
    ax2.plot(t, o.median, color=PALETTE[2], lw=1.6, label="median")
    ax2.axhline(0.5, ls=":", color="grey", lw=1)
    mos = ens.metrics["median_os_weeks"]
    ax2.set_title(f"Population OS — median {mos['median']:.0f} wk  [{mos['lo']:.0f}, {mos['hi']:.0f}]")
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("survival fraction")
    ax2.set_ylim(0, 1.02)
    ax2.legend(fontsize=8)

    w8 = ens.metrics["week8_relative_change"]
    fig.suptitle(
        "Parameter uncertainty (Claret NSCLC) — the ~90% CV resistance/kill terms give a wide "
        f"week-8 change band [{w8['lo'] * 100:.0f}%, {w8['hi'] * 100:.0f}%]; n={ens.n}",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "uncertainty.png", dpi=120)
    plt.close(fig)


def tgi_metrics_figure() -> None:
    """Annotated TGI-metric panel extracted from a simulated trajectory."""
    ds = onkos.load()
    t = np.linspace(0.0, 156.0, 313)
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    tr = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential", context=ctx, drug_effect=1.2, t=t)
    v = tr.tumor_size
    m = tr.metrics
    y0 = v[0]

    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.semilogy(t, v, color=PALETTE[0], lw=1.8, label="tumor size (log)")
    ax.axhline(y0, ls=":", color="grey", lw=1)
    ax.axhline(0.7 * y0, ls="--", color="#2f855a", lw=1, label="RECIST PR (−30% from baseline)")

    nadir, tnad = m["nadir_tumor_size"], m["time_to_nadir_weeks"]
    ax.axhline(1.2 * nadir, ls="--", color="#c53030", lw=1, label="RECIST PD (+20% from nadir)")
    ax.scatter([tnad], [nadir], color="#c53030", zorder=5)
    ax.annotate(
        f"nadir / time-to-growth = {tnad:.0f} wk\ndepth of response = {m['depth_of_response'] * 100:.0f}%",
        (tnad, nadir), textcoords="offset points", xytext=(10, 14), fontsize=8,
    )
    ax.annotate(f"shrink phase\nk_s ≈ {m['tumor_shrinkage_rate_ks']:.3f}/wk",
                (1.0, y0 * 0.9), fontsize=8, color="#2f855a", va="top")
    ax.annotate(f"regrowth phase\nk_g ≈ {m['tumor_growth_rate_kg']:.3f}/wk",
                (t[-1] * 0.62, v[-1] * 0.5), fontsize=8, color="#c05621")
    dor = m["duration_of_response_weeks"]
    title_dor = f"duration of response = {dor:.0f} wk" if np.isfinite(dor) else "no RECIST PR"
    ax.set_title(f"TGI-metric extraction (Stein/Bruno panel) — {title_dor}", fontsize=10)
    ax.set_xlabel("weeks")
    ax.set_ylabel("tumor size (mm, SLD, log)")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "tgi_metrics.png", dpi=120)
    plt.close(fig)


def sensitivity_figure() -> None:
    """Tornado: which parameter's IIV drives the survival prediction."""
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    res = onkos.sensitivity(ds, "resistance.claret_2009.tgi", context=ctx,
                            target="median_os_weeks", n=600, seed=0)
    rows = list(reversed(res.indices))  # largest at top
    labels = [f"{p.symbol}\n(CV {p.iiv_cv_percent:.0f}%)" for p in rows]
    contrib = [p.contribution * 100 for p in rows]
    colors = ["#2f855a" if p.src > 0 else "#c53030" for p in rows]

    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    ax.barh(labels, contrib, color=colors)
    for i, p in enumerate(rows):
        ax.text(p.contribution * 100 + 1, i, f"{p.contribution * 100:.0f}%  (SRC {p.src:+.2f})",
                va="center", fontsize=8)
    ax.set_xlim(0, 105)
    ax.set_xlabel("contribution to median-OS variance (%)")
    ax.set_title(
        f"Parameter sensitivity (Claret NSCLC) — first-order R²={res.r_squared:.2f}; "
        f"verify '{res.dominant.symbol}' first",
        fontsize=10,
    )
    # legend for sign
    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(color="#2f855a", label="↑ param → ↑ OS"),
                 Patch(color="#c53030", label="↑ param → ↓ OS")],
        fontsize=8, loc="lower right",
    )
    fig.tight_layout()
    fig.savefig(OUT / "sensitivity.png", dpi=120)
    plt.close(fig)


def survival_endpoints_figure() -> None:
    """OS and PFS population curves from the same TGI metric (PFS < OS)."""
    ds = onkos.load()
    t = np.linspace(0.0, 156.0, 313)
    ctx = {"tumor_type": "NSCLC", "line": "first"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    for i, rid in enumerate(["resistance.claret_2009.tgi", "tgi_metrics.wang_2009.biexponential"]):
        tr = onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t)
        c = PALETTE[i]
        label = rid.split(".")[1]
        ax1.plot(t, tr.os_curve, color=c, lw=1.7, label=f"{label} OS")
        ax1.plot(t, tr.pfs_curve, color=c, lw=1.4, ls="--", label=f"{label} PFS")
    ax1.axhline(0.5, ls=":", color="grey", lw=1)
    ax1.set_title("OS (solid) vs PFS (dashed) — NSCLC, 1L")
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("survival fraction")
    ax1.set_ylim(0, 1.02)
    ax1.legend(fontsize=7)

    # median OS vs PFS by tumor context
    contexts = ["NSCLC", "breast", "CRC", "HCC", "melanoma"]
    rid_for = {  # an in-context TGI model per tumor type
        "NSCLC": "resistance.claret_2009.tgi",
        "breast": "tgi_metrics.bruno_2020.breast_biexponential",
        "CRC": "resistance.crc_first_line.claret",
        "HCC": "resistance.hcc_first_line.claret",
        "melanoma": "resistance.melanoma_first_line.claret",
    }
    tlong = np.linspace(0.0, 312.0, 625)
    os_meds, pfs_meds = [], []
    for tt in contexts:
        tr = onkos.simulate(ds, rid_for[tt], context={"tumor_type": tt, "line": "first"},
                            drug_effect=1.0, t=tlong)
        os_meds.append(tr.median_os or np.nan)
        pfs_meds.append(tr.median_pfs or np.nan)
    x = np.arange(len(contexts))
    ax2.bar(x - 0.2, os_meds, 0.4, color=PALETTE[0], label="median OS")
    ax2.bar(x + 0.2, pfs_meds, 0.4, color=PALETTE[1], label="median PFS")
    ax2.set_xticks(x, contexts, fontsize=8)
    ax2.set_ylabel("weeks")
    ax2.set_title("Median OS vs PFS by tumor context")
    ax2.legend(fontsize=8)

    fig.suptitle("Survival endpoints (spec §2, §6): parametric OS and PFS links from the week-8 "
                 "TGI metric", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "survival_endpoints.png", dpi=120)
    plt.close(fig)


def survival_model_choice_figure() -> None:
    """Weibull (parametric) vs Cox (nonparametric baseline) OS from the same metric."""
    ds = onkos.load()
    t = np.linspace(0.0, 208.0, 417)
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    w = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0, t=t)
    c = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0, t=t,
                       survival_link="survival_link.nsclc_os_cox")
    cox = ds["survival_link.nsclc_os_cox"].structure["baseline_survival"]

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.plot(t, w.os_curve, color=PALETTE[0], lw=1.8, label="parametric Weibull-PH OS")
    ax.plot(t, c.os_curve, color=PALETTE[1], lw=1.8, label="Cox-PH OS (tabulated baseline)")
    ax.plot(cox["times"], cox["survival"], "o", color=PALETTE[2], ms=4, alpha=0.7,
            label="Cox baseline S₀(t) (x=0)")
    ax.fill_between(t, np.minimum(w.os_curve, c.os_curve), np.maximum(w.os_curve, c.os_curve),
                    color="grey", alpha=0.15)
    ax.axhline(0.5, ls=":", color="grey", lw=1)
    spread = float(np.max(np.abs(w.os_curve - c.os_curve)))
    ax.set_title(f"Survival-model choice (NSCLC OS, same week-8 metric) — "
                 f"max spread {spread:.2f}; median {w.median_os:.0f} vs {c.median_os:.0f} wk",
                 fontsize=10)
    ax.set_xlabel("weeks")
    ax.set_ylabel("survival fraction")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "survival_model_choice.png", dpi=120)
    plt.close(fig)


def line_of_therapy_figure() -> None:
    """First vs second line: shorter survival + the line-aware model set (NSCLC)."""
    ds = onkos.load()
    t = np.linspace(0.0, 208.0, 417)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    for i, line in enumerate(["first", "second"]):
        tr = onkos.simulate(ds, "tgi_metrics.wang_2009.biexponential",
                            context={"tumor_type": "NSCLC", "line": line}, drug_effect=1.0, t=t)
        c = PALETTE[i]
        ax1.plot(t, tr.os_curve, color=c, lw=1.8, label=f"{line} line OS")
        ax1.plot(t, tr.pfs_curve, color=c, lw=1.3, ls="--", label=f"{line} line PFS")
    ax1.axhline(0.5, ls=":", color="grey", lw=1)
    ax1.set_title("Same model, two lines — second-line survival is shorter")
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("survival fraction")
    ax1.set_ylim(0, 1.02)
    ax1.legend(fontsize=7)

    # median OS/PFS by line, with the eligible model count annotated
    lines = ["first", "second"]
    os_meds, pfs_meds, counts = [], [], []
    for line in lines:
        cmp = onkos.compare(ds, purpose="tgi", context={"tumor_type": "NSCLC", "line": line},
                            drug_effect=1.0, t=t)
        curves_os = [tr.median_os for tr in cmp.included if tr.median_os]
        curves_pfs = [tr.median_pfs for tr in cmp.included if tr.median_pfs]
        os_meds.append(np.mean(curves_os) if curves_os else 0)
        pfs_meds.append(np.mean(curves_pfs) if curves_pfs else 0)
        counts.append(len(cmp.included))
    x = np.arange(len(lines))
    ax2.bar(x - 0.2, os_meds, 0.4, color=PALETTE[0], label="mean median OS")
    ax2.bar(x + 0.2, pfs_meds, 0.4, color=PALETTE[1], label="mean median PFS")
    for i, cnt in enumerate(counts):
        ax2.text(i, max(os_meds) * 1.02, f"{cnt} models", ha="center", fontsize=8)
    ax2.set_xticks(x, [f"{x_} line" for x_ in lines])
    ax2.set_ylabel("weeks")
    ax2.set_title("NSCLC median OS / PFS by line of therapy")
    ax2.legend(fontsize=8)

    fig.suptitle("Line of therapy (Phase C): line-aware survival matching — a 2L context never "
                 "borrows a 1L model", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "line_of_therapy.png", dpi=120)
    plt.close(fig)


def composability_chain_figure() -> None:
    """The full PK -> exposure -> tumor-dynamics -> survival chain (Hypnos compose)."""
    from onkos import pk
    from onkos.export.reference import effect
    from onkos.export.registry import get_kernel, kernel_values

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    er_id = "exposure_response.emax_generic"
    er = ds[er_id]
    er_spec, er_vals = get_kernel(er), kernel_values(er)
    pkargs = dict(ka=0.5, ke=0.05, v=5.0, f=0.8)
    doses = [600, 1200, 2400]
    t = np.linspace(0.0, 104.0, 209)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    (a_pk, a_er), (a_tum, a_os) = axes

    # 1) PK concentration profile (multiple-dose) for the middle dose + C_avg
    tp = np.linspace(0, 168, 600)  # one week, hourly-ish
    conc = pk.concentration_profile(1200, 24, 7, **pkargs, t=tp)
    cavg_mid = pk.steady_state_metrics(dose=1200, tau=24, **pkargs)["c_avg"]
    a_pk.plot(tp, conc, color=PALETTE[0])
    a_pk.axhline(cavg_mid, ls="--", color=PALETTE[1], label=f"C_avg = {cavg_mid:.0f} µg/L")
    a_pk.set_title("1. PK: oral regimen → concentration C(t)")
    a_pk.set_xlabel("hours")
    a_pk.set_ylabel("C (µg/L)")
    a_pk.legend(fontsize=8)

    # 2) ER transform C -> E, with each dose's C_avg operating point
    c_grid = np.linspace(0, 400, 200)
    a_er.plot(c_grid, effect(er_spec, c_grid, er_vals), color=PALETTE[2])
    for i, dose in enumerate(doses):
        cavg = pk.steady_state_metrics(dose=dose, tau=24, **pkargs)["c_avg"]
        e = float(effect(er_spec, cavg, er_vals))
        a_er.plot([cavg], [e], "o", color=PALETTE[i % 5], ms=7)
        a_er.annotate(f"{dose} mg", (cavg, e), textcoords="offset points", xytext=(4, -10), fontsize=7)
    a_er.set_title("2. Exposure-response: C_avg → effect E")
    a_er.set_xlabel("C_avg (µg/L)")
    a_er.set_ylabel("drug effect E")

    # 3) + 4) tumor dynamics and OS per dose
    for i, dose in enumerate(doses):
        cavg = pk.steady_state_metrics(dose=dose, tau=24, **pkargs)["c_avg"]
        tr = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx,
                            exposure=cavg, exposure_response=er_id, t=t)
        c = PALETTE[i % 5]
        a_tum.plot(t, tr.tumor_size, color=c, label=f"{dose} mg")
        a_os.plot(t, tr.os_curve, color=c, label=f"{dose} mg (mOS {tr.median_os:.0f})"
                  if tr.median_os else f"{dose} mg")
    a_tum.set_title("3. Tumor dynamics driven by E")
    a_tum.set_xlabel("weeks")
    a_tum.set_ylabel("tumor size (mm, SLD)")
    a_tum.legend(fontsize=8)
    a_os.axhline(0.5, ls=":", color="grey", lw=1)
    a_os.set_title("4. Population OS")
    a_os.set_xlabel("weeks")
    a_os.set_ylabel("survival fraction")
    a_os.set_ylim(0, 1.02)
    a_os.legend(fontsize=8)

    fig.suptitle("Hypnos composability: dose → exposure → tumor dynamics → survival, "
                 "one open tier-annotated chain (illustrative PK)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / "composability_chain.png", dpi=120)
    plt.close(fig)


def model_average_figure() -> None:
    """Model averaging (Axis 3): the averaged OS curve with its between-model band,
    and the within/between variance split ranking contexts by model-choice risk."""
    ds = onkos.load()
    t = np.linspace(0.0, 156.0, 313)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: NSCLC 1L — per-model OS curves, the model-averaged S̄(t), and its
    # between-model ±1σ band (the disagreement that travels with the average).
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    cmp = onkos.compare(ds, purpose="tgi", context=ctx, drug_effect=1.0, t=t)
    ma = cmp.model_average(target="median_os_weeks", endpoint="OS", weights="equal", n=300)
    for i, tr in enumerate(cmp.included):
        ax1.plot(t, tr.os_curve, color=PALETTE[i % len(PALETTE)], lw=1.0, alpha=0.7,
                 label=f"{tr.record_id.split('.')[1]} [{tr.tier}]")
    band = ma.between_band
    ax1.fill_between(t, np.clip(ma.curve - band, 0, 1), np.clip(ma.curve + band, 0, 1),
                     color="black", alpha=0.12, label="between-model ±1σ")
    ax1.plot(t, ma.curve, color="black", lw=2.2, label=f"model average [{ma.tier}]")
    ax1.axhline(0.5, ls=":", color="grey", lw=1)
    ax1.set_title(
        f"Model-averaged OS (NSCLC 1L) — median {ma.point:.0f} wk, "
        f"model-selection {ma.model_selection_fraction:.0%}",
        fontsize=9,
    )
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("survival fraction")
    ax1.set_ylim(0, 1.02)
    ax1.legend(fontsize=7)

    # Right: per-context within vs between variance (stacked) — the between share
    # IS the irreducible model-choice risk that ranks curation value.
    contexts = [("NSCLC", "first"), ("NSCLC", "second"), ("breast", "first"),
                ("CRC", "first"), ("HCC", "first"), ("melanoma", "first")]
    labels, within, between, fracs = [], [], [], []
    for tt, ln in contexts:
        c = onkos.compare(ds, purpose="tgi", context={"tumor_type": tt, "line": ln},
                          drug_effect=1.0, t=t)
        if len(c.included) < 2:
            continue
        m = c.model_average(target="median_os_weeks", endpoint="OS", weights="equal", n=200)
        labels.append(f"{tt}\n{ln[:2]}")
        within.append(m.within_var)
        between.append(m.between_var)
        fracs.append(m.model_selection_fraction)
    order = np.argsort(fracs)[::-1]
    labels = [labels[i] for i in order]
    within = np.array([within[i] for i in order])
    between = np.array([between[i] for i in order])
    fracs = [fracs[i] for i in order]
    x = np.arange(len(labels))
    ax2.bar(x, within, color="#3182ce", label="within (parameter noise)")
    ax2.bar(x, between, bottom=within, color="#c53030", label="between (model-selection)")
    for i, f in enumerate(fracs):
        ax2.text(i, within[i] + between[i], f"{f:.0%}", ha="center", va="bottom", fontsize=8)
    ax2.set_xticks(x, labels, fontsize=8)
    ax2.set_ylabel("variance of median OS (week²)")
    ax2.set_title("Uncertainty split by context (label = model-selection %)", fontsize=9)
    ax2.legend(fontsize=8)

    fig.suptitle(
        "Model-selection uncertainty (Axis 3): the divergence becomes a variance "
        "decomposition + an honestly-weighted central forecast — NOT a clinical prediction",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "model_average.png", dpi=120)
    plt.close(fig)


def identifiability_figure() -> None:
    """Practical identifiability (design axis): predicted RSE vs stored IIV CV, and
    how follow-up length identifies (or fails to identify) each parameter."""
    from onkos.identify import identifiability

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    rid = "resistance.claret_2009.tgi"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: per-parameter predicted RSE (Cramér-Rao) next to the stored IIV CV under
    # a realistic RECIST cadence. The kill rate kD is well identified; the growth and
    # resistance terms are flat — so λ's ~90% CV is partly a flat-likelihood artifact.
    res = identifiability(ds, rid, context=ctx)
    order = ["kL", "kD", "lambda"]  # growth, kill, resistance (kernel order)
    by = {p.symbol: p for p in res.params}
    rse = [by[s].rse_percent for s in order]
    cv = [by[s].iiv_cv_percent or 0 for s in order]
    x = np.arange(len(order))
    ax1.bar(x - 0.2, rse, 0.4, color="#c53030", label="predicted RSE (design)")
    ax1.bar(x + 0.2, cv, 0.4, color="#3182ce", label="stored IIV CV")
    ax1.axhline(res.rse_ceiling_percent, ls=":", color="grey", lw=1)
    ax1.text(2.3, res.rse_ceiling_percent * 1.05, "RSE ceiling 50%", fontsize=7, color="grey",
             ha="right")
    ax1.set_yscale("log")
    ax1.set_xticks(x, ["kL\ngrowth", "kD\nkill", "λ\nresistance"], fontsize=8)
    ax1.set_ylabel("percent (log)")
    ax1.set_title("Predicted RSE vs stored IIV CV (NSCLC, RECIST cadence)", fontsize=9)
    ax1.legend(fontsize=8)
    ax1.annotate("flat likelihood:\nCV is partly an artifact", (2, max(rse[2], cv[2])),
                 textcoords="offset points", xytext=(-4, 18), fontsize=7, color="#c53030",
                 ha="center")

    # Right: predicted RSE for each parameter as the trial follow-up lengthens
    # (q6w to the horizon). kD is identifiable from any short trial; λ only drops
    # below the ceiling once follow-up runs past resistance-driven regrowth; kL
    # (growth, masked by treatment) stays the hardest to pin down.
    horizons = [24, 36, 48, 72, 104, 156]
    curves = {s: [] for s in order}
    for h in horizons:
        sched = list(np.arange(0.0, h + 1e-9, 6.0))
        r = identifiability(ds, rid, context=ctx, schedule=sched)
        b = {p.symbol: p for p in r.params}
        for s in order:
            curves[s].append(b[s].rse_percent)
    labels = {"kL": "kL growth", "kD": "kD kill", "lambda": "λ resistance"}
    colors = {"kL": "#b7791f", "kD": "#2f855a", "lambda": "#c53030"}
    for s in order:
        ax2.plot(horizons, curves[s], "o-", color=colors[s], label=labels[s], ms=4)
    ax2.axhline(50, ls=":", color="grey", lw=1)
    ax2.text(26, 56, "identifiable below here", fontsize=7, color="grey")
    ax2.set_yscale("log")
    ax2.set_xlabel("trial follow-up horizon (weeks, q6w sampling)")
    ax2.set_ylabel("predicted RSE % (log)")
    ax2.set_title("Identifiability vs follow-up length", fontsize=9)
    ax2.legend(fontsize=8)

    fig.suptitle(
        "Practical identifiability (design axis): the Fisher information says which "
        "parameters a realistic trial can estimate — NOT a patient-level quantity",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "identifiability.png", dpi=120)
    plt.close(fig)


def combination_interaction_figure() -> None:
    """Drug-combination interaction as a model-selection axis: how the combined effect
    and the predicted OS depend on the (unmeasured) interaction assumption."""
    from onkos.interaction import combine_effects, compare_interactions

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    rid = "resistance.claret_2009.tgi"
    ea = eb = 0.6

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: combined effect E_AB vs the interaction parameter ψ, with the HSA and
    # additive nulls as horizontal references. ψ is a declared assumption, not fitted.
    psis = np.linspace(-1.0, 2.0, 121)
    greco = [combine_effects(ea, eb, model="greco", psi=p) for p in psis]
    hsa = combine_effects(ea, eb, model="hsa")
    add = combine_effects(ea, eb, model="additive")
    ax1.plot(psis, greco, color="#6b46c1", lw=1.8, label="greco E_AB(ψ)")
    ax1.axhline(add, ls="--", color="#2b6cb0", lw=1.2, label=f"additive / Bliss null = {add:g}")
    ax1.axhline(hsa, ls=":", color="#2f855a", lw=1.2, label=f"HSA (highest single agent) = {hsa:g}")
    ax1.axvline(0, color="grey", lw=0.8)
    ax1.fill_betweenx([0, max(greco)], -1.0, 0, color="#c53030", alpha=0.06)
    ax1.fill_betweenx([0, max(greco)], 0, 2.0, color="#2f855a", alpha=0.06)
    ax1.text(-0.92, max(greco) * 0.95, "antagonism", fontsize=7, color="#c53030")
    ax1.text(1.2, max(greco) * 0.95, "synergy", fontsize=7, color="#2f855a")
    ax1.set_title(f"Combined effect of E_A={ea:g} + E_B={eb:g} vs interaction ψ", fontsize=9)
    ax1.set_xlabel("interaction parameter ψ (DECLARED assumption)")
    ax1.set_ylabel("combined drug effect E_AB")
    ax1.legend(fontsize=7, loc="upper left")

    # Right: population OS under each interaction model for the SAME single-agent
    # activity — the survival divergence driven purely by the interaction assumption.
    t = np.linspace(0.0, 156.0, 313)
    cmp = compare_interactions(ds, rid, context=ctx, effect_a=ea, effect_b=eb, psi=0.5, t=t)
    colors = {"hsa": "#2f855a", "additive": "#2b6cb0", "greco+0.5": "#6b46c1",
              "greco-0.5": "#c05621"}
    for label, tr in cmp.trajectories.items():
        mos = f" (mOS {tr.median_os:.0f})" if tr.median_os else ""
        ax2.plot(t, tr.os_curve, color=colors.get(label, "grey"), lw=1.7,
                 label=f"{label}{mos}")
    ax2.axhline(0.5, ls=":", color="grey", lw=1)
    ax2.set_ylim(0, 1.02)
    ax2.set_title(f"Population OS by interaction model — divergence {cmp.os_divergence:.2f}",
                  fontsize=9)
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("survival fraction")
    ax2.legend(fontsize=7)

    rng = cmp.median_os_range
    fig.suptitle(
        "Combination therapy: the interaction model is itself a model-selection axis — same "
        f"single-agent activity, median OS {rng[0]:.0f}-{rng[1]:.0f} wk by assumption",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "combination_interaction.png", dpi=120)
    plt.close(fig)


def two_population_resistance_figure() -> None:
    """Mechanistic resistance: the sensitive/resistant clone split that produces the
    nadir-then-regrowth, and the phenomenological-vs-mechanistic resistance divergence."""
    from onkos.export.reference import init_vector
    from onkos.export.registry import get_kernel, kernel_values
    from scipy.integrate import solve_ivp

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    rid = "resistance.nsclc_first_line.two_population"
    t = np.linspace(0.0, 156.0, 313)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: the two clones and the observed tumor. The drug crushes the sensitive
    # clone (nadir), then the untouched resistant clone outgrows — the mechanistic
    # origin of resistance-driven regrowth, with R0 a biologically interpretable seed.
    spec = get_kernel(ds[rid])
    vals = kernel_values(ds[rid])
    vals["V0"] = onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t).tumor_size[0] - vals["R0"]
    vals["E"] = 1.0
    sol = solve_ivp(lambda tt, yy: spec.rhs(tt, yy, vals), (t[0], t[-1]),
                    init_vector(spec, vals), t_eval=t, rtol=1e-8, atol=1e-10, method="LSODA")
    s, r = sol.y
    ax1.semilogy(t, s, color="#2b6cb0", lw=1.6, label="sensitive clone S (killed)")
    ax1.semilogy(t, r, color="#c53030", lw=1.6, label="resistant clone R (outgrows)")
    ax1.semilogy(t, s + r, color="black", lw=2.0, ls="--", label="observed tumor V = S + R")
    nadir_i = int(np.argmin(s + r))
    ax1.scatter([t[nadir_i]], [(s + r)[nadir_i]], color="black", zorder=5)
    ax1.annotate("nadir, then\nresistant-driven regrowth", (t[nadir_i], (s + r)[nadir_i]),
                 textcoords="offset points", xytext=(14, 6), fontsize=8)
    ax1.set_title("Two-population resistance — clone decomposition", fontsize=9)
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("tumor size (mm, SLD, log)")
    ax1.legend(fontsize=7, loc="lower right")

    # Right: same context, two resistance MECHANISMS — phenomenological decay-of-effect
    # (Claret λ) vs mechanistic resistant-subclone. Tuned to share the early kill, they
    # agree at week 8 (and hence on the week-8-driven OS) yet diverge in the tumor TAIL:
    # the resistance-model risk is real but nearly invisible to a short-trial surrogate.
    mech = onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t)
    claret = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0, t=t)
    ax2.semilogy(t, claret.tumor_size, color="#c05621", lw=1.8,
                 label=f"phenomenological (Claret λ)  mOS {claret.median_os:.0f}")
    ax2.semilogy(t, mech.tumor_size, color="#2f855a", lw=1.8,
                 label=f"mechanistic (resistant subclone)  mOS {mech.median_os:.0f}")
    ax2.fill_between(t, np.minimum(claret.tumor_size, mech.tumor_size),
                     np.maximum(claret.tumor_size, mech.tumor_size), color="grey", alpha=0.15)
    ax2.axvline(8, ls=":", color="grey", lw=1)
    ax2.text(9, ax2.get_ylim()[1] * 0.5, "week 8\n(agree → same OS)", fontsize=7, color="grey")
    ratio = mech.tumor_size[-1] / claret.tumor_size[-1]
    ax2.set_title(f"Resistance model as a divergence axis — {ratio:.1f}x tumor split by 3 yr",
                  fontsize=9)
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("tumor size (mm, SLD, log)")
    ax2.legend(fontsize=8, loc="lower right")

    fig.suptitle(
        "Mechanistic resistance (Goldie-Coldman): the resistance MODEL is a model-selection "
        "axis — same early kill & OS, divergent tumor tail",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "two_population_resistance.png", dpi=120)
    plt.close(fig)


def survival_metric_choice_figure() -> None:
    """Which on-treatment metric drives the OS link is a model-selection axis: the
    early week-8 surrogate vs the tail-sensitive growth-rate constant k_g re-rank —
    and invert — which model looks better."""
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0.0, 260.0, 521)
    models = [
        ("resistance.claret_2009.tgi", "Claret (phenom. resistance)", "#c05621"),
        ("resistance.nsclc_first_line.two_population", "two-population (mechanistic)", "#2f855a"),
        ("drug_effect.norton_simon.nsclc", "Norton-Simon (complete responder)", "#2b6cb0"),
        ("tgi_metrics.wang_2009.biexponential", "Wang biexponential", "#6b46c1"),
    ]
    kg_link = "survival_link.nsclc_os_growth_rate"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)

    for rid, label, color in models:
        w = onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t)            # week-8 link
        g = onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t, survival_link=kg_link)
        wm = f"{w.median_os:.0f}" if w.median_os else "n/r"
        gm = f"{g.median_os:.0f}" if g.median_os else "n/r"
        ax1.plot(t, w.os_curve, color=color, lw=1.8, label=f"{label}  (mOS {wm})")
        ax2.plot(t, g.os_curve, color=color, lw=1.8, label=f"{label}  (mOS {gm})")

    for ax, title in (
        (ax1, "OS from week-8 change (default surrogate)"),
        (ax2, "OS from growth-rate constant k_g (tail-sensitive)"),
    ):
        ax.axhline(0.5, ls=":", color="grey", lw=1)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("weeks")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=7, loc="upper right")
    ax1.set_ylabel("survival fraction")

    fig.suptitle(
        "Survival-metric choice is a model-selection axis: week-8 says two-population > Claret "
        "and ranks the complete responder last; k_g inverts both",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "survival_metric_choice.png", dpi=120)
    plt.close(fig)


def model_selection_budget_figure() -> None:
    """The model-selection budget: total forecast variance split across the structural
    choices (TGI model, survival link, their interaction) vs parameter noise — the
    capstone of the model-selection arc."""
    from onkos.budget import model_selection_budget

    ds = onkos.load()
    comp_colors = {
        "parameter": "#3182ce",
        "tgi_model": "#b7791f",
        "survival_link": "#2f855a",
        "interaction": "#c53030",
    }
    comp_labels = {
        "parameter": "parameter noise (IIV)",
        "tgi_model": "TGI-model choice",
        "survival_link": "survival-link choice",
        "interaction": "model × link interaction",
    }
    order = ["parameter", "tgi_model", "survival_link", "interaction"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: NSCLC 1L — the full 4x3 grid, one stacked bar of the four components.
    b = model_selection_budget(ds, context={"tumor_type": "NSCLC", "line": "first"},
                               endpoint="OS", n=300)
    left = 0.0
    for k in order:
        frac = b.fractions[k]
        ax1.barh([0], [frac], left=left, color=comp_colors[k], label=comp_labels[k])
        if frac > 0.04:
            ax1.text(left + frac / 2, 0, f"{frac * 100:.0f}%", ha="center", va="center",
                     color="white", fontsize=9, fontweight="bold")
        left += frac
    ax1.axvline(b.fractions["parameter"], color="black", lw=1.2, ls=":")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-0.5, 0.5)
    ax1.set_yticks([])
    ax1.set_xlabel("share of total forecast variance")
    ax1.set_title(
        f"NSCLC 1L OS budget — {b.structural_fraction * 100:.0f}% structural "
        f"(irreducible), {b.fractions['parameter'] * 100:.0f}% parameter",
        fontsize=9,
    )
    ax1.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)

    # Right: per-context parameter-vs-structural split, ranked by structural share.
    contexts = [("NSCLC", "first"), ("CRC", "first"), ("HCC", "first"), ("melanoma", "first"),
                ("breast", "first"), ("NSCLC", "second")]
    rows = []
    for tt, ln in contexts:
        bb = model_selection_budget(ds, context={"tumor_type": tt, "line": ln},
                                    endpoint="OS", n=120)
        rows.append((f"{tt}\n{ln[:2]}", bb.fractions["parameter"], bb.structural_fraction))
    rows.sort(key=lambda r: r[2], reverse=True)
    labels = [r[0] for r in rows]
    param = np.array([r[1] for r in rows])
    struct = np.array([r[2] for r in rows])
    x = np.arange(len(labels))
    ax2.bar(x, param, color="#3182ce", label="parameter (reducible by more data)")
    ax2.bar(x, struct, bottom=param, color="#c53030", label="structural (irreducible)")
    ax2.axhline(0.5, ls=":", color="grey", lw=1)
    for i, s in enumerate(struct):
        ax2.text(i, param[i] + s + 0.01, f"{s * 100:.0f}%", ha="center", fontsize=8)
    ax2.set_xticks(x, labels, fontsize=8)
    ax2.set_ylim(0, 1.12)
    ax2.set_ylabel("share of forecast variance")
    ax2.set_title("Structural vs parameter share by context", fontsize=9)
    ax2.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1)

    fig.suptitle(
        "The model-selection budget (capstone): which structural assumption — not the "
        "parameters — drives the forecast, and where standardization buys the most",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "model_selection_budget.png", dpi=120)
    plt.close(fig)


def response_orr_figure() -> None:
    """RECIST best-response distribution per model, and the ORR -> OS surrogate that is
    concordant under the week-8 link but discordant under the tail-sensitive k_g link."""
    from onkos.response import objective_response_rate, response_vs_survival

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0.0, 156.0, 313)
    models = [
        ("drug_effect.norton_simon.nsclc", "Norton-Simon"),
        ("resistance.claret_2009.tgi", "Claret"),
        ("resistance.nsclc_first_line.two_population", "two-population"),
        ("tgi_metrics.wang_2009.biexponential", "Wang biexp"),
    ]
    cat_colors = {"CR": "#2f855a", "PR": "#68d391", "SD": "#f6ad55", "PD": "#c53030"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: stacked RECIST best-response distribution (the population readout) per model.
    labels, dists, orrs = [], [], []
    for rid, name in models:
        rr = objective_response_rate(ds, rid, context=ctx, t=t, n=400)
        labels.append(f"{name}\nORR {rr.orr * 100:.0f}%")
        dists.append(rr.distribution)
        orrs.append(rr.orr)
    x = np.arange(len(labels))
    bottoms = np.zeros(len(labels))
    for cat in ("CR", "PR", "SD", "PD"):
        vals = np.array([d[cat] for d in dists])
        ax1.bar(x, vals, bottom=bottoms, color=cat_colors[cat], label=cat)
        bottoms += vals
    ax1.axhline(0.5, ls=":", color="grey", lw=1)
    ax1.set_xticks(x, labels, fontsize=8)
    ax1.set_ylim(0, 1.0)
    ax1.set_ylabel("fraction of trial (population)")
    ax1.set_title("RECIST best-response distribution (NSCLC 1L)", fontsize=9)
    ax1.legend(fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.16))

    # Right: ORR vs median OS under two survival links — the conditional surrogate.
    rs_w = response_vs_survival(ds, context=ctx, t=t, n=400)
    rs_k = response_vs_survival(ds, context=ctx, survival_link="survival_link.nsclc_os_growth_rate",
                               t=t, n=400)
    by_w = {r["record_id"]: r for r in rs_w.rows}
    by_k = {r["record_id"]: r for r in rs_k.rows}
    for rid, name in models:
        o = by_w[rid]["orr"]
        ax2.plot([o], [by_w[rid]["median_os_weeks"]], "o", color="#2b6cb0", ms=8)
        ax2.plot([o], [by_k[rid]["median_os_weeks"]], "s", color="#c53030", ms=8)
        ax2.annotate(name, (o, by_w[rid]["median_os_weeks"]), textcoords="offset points",
                     xytext=(6, 4), fontsize=7, color="#2b6cb0")
    # Trend lines (sorted by ORR) to show concordance vs inversion.
    order = sorted(models, key=lambda m: by_w[m[0]]["orr"])
    ox = [by_w[m[0]]["orr"] for m in order]
    ax2.plot(ox, [by_w[m[0]]["median_os_weeks"] for m in order], color="#2b6cb0", lw=1,
             alpha=0.6, label=f"week-8 link (concordant, {rs_w.discordant_pairs}/{rs_w.total_pairs})")
    ax2.plot(ox, [by_k[m[0]]["median_os_weeks"] for m in order], color="#c53030", lw=1,
             alpha=0.6, label=f"k_g link (discordant, {rs_k.discordant_pairs}/{rs_k.total_pairs})")
    ax2.set_xlabel("ORR (objective response rate)")
    ax2.set_ylabel("median OS (weeks)")
    ax2.set_title("ORR → OS surrogate: faithful under week-8, inverted under k_g", fontsize=9)
    ax2.legend(fontsize=7, loc="center right")

    fig.suptitle(
        "RECIST response & the ORR→OS surrogate: the phase-2 endpoint can rank models the "
        "opposite way survival does — conditional on the survival mechanism",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "response_orr.png", dpi=120)
    plt.close(fig)


def duration_of_response_figure() -> None:
    """Depth is not durability: ORR (breadth) and median DoR (durability) dissociate, and
    under the tail-sensitive k_g link durability tracks survival where breadth inverts it."""
    from onkos.response import response_vs_survival

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0.0, 208.0, 417)
    names = {
        "drug_effect.norton_simon.nsclc": "Norton-Simon",
        "resistance.claret_2009.tgi": "Claret",
        "resistance.nsclc_first_line.two_population": "two-population",
        "resistance.nsclc_first_line.acquired": "acquired",
        "tgi_metrics.wang_2009.biexponential": "Wang biexp",
    }

    def label_of(rid):
        return names.get(rid, rid.split(".")[1])
    # k_g link: the tail-sensitive endpoint where ORR mis-ranks OS (v0.27).
    rs = response_vs_survival(ds, context=ctx, survival_link="survival_link.nsclc_os_growth_rate",
                              t=t, n=500)
    rows = [r for r in rs.rows if r["median_dor_weeks"] is not None]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: ORR (breadth) vs median DoR (durability), coloured by tail-driven OS.
    orr = np.array([r["orr"] for r in rows])
    dor = np.array([r["median_dor_weeks"] for r in rows])
    osv = np.array([r["median_os_weeks"] for r in rows])
    sc = ax1.scatter(orr, dor, c=osv, cmap="viridis", s=140, edgecolor="black", zorder=3)
    for r in rows:
        ax1.annotate(label_of(r["record_id"]), (r["orr"], r["median_dor_weeks"]),
                     textcoords="offset points", xytext=(8, 4), fontsize=8)
    ax1.axhspan(0, 35, xmin=0.55, color="#c53030", alpha=0.07)
    ax1.text(0.99, 24, "broad but brief\n(durability failure)", fontsize=7, color="#c53030",
             ha="right")
    ax1.set_xlabel("ORR — breadth (how many respond)")
    ax1.set_ylabel("median DoR — durability (weeks)")
    ax1.set_title("Breadth and durability dissociate", fontsize=9)
    fig.colorbar(sc, ax=ax1, shrink=0.85, label="median OS under k_g (wk)")

    # Right: per model sorted by tail-driven OS — the highest-ORR (broadest) model is the
    # worst survivor, and its responses are among the briefest; the best survivor's
    # responses are durable. Breadth (ORR) does not order OS; durability (DoR) does better.
    order = sorted(rows, key=lambda r: r["median_os_weeks"])
    labels = [f"{label_of(r['record_id'])}\nOS {r['median_os_weeks']:.0f}" for r in order]
    x = np.arange(len(order))
    ax2.bar(x - 0.2, [r["orr"] for r in order], 0.4, color="#c53030", label="ORR (breadth)")
    ax2b = ax2.twinx()
    ax2b.bar(x + 0.2, [r["median_dor_weeks"] for r in order], 0.4, color="#2b6cb0",
             label="median DoR (durability)")
    ax2.set_xticks(x, labels, fontsize=7)
    ax2.set_ylabel("ORR", color="#c53030")
    ax2.set_ylim(0, 1.1)
    ax2b.set_ylabel("median DoR (wk)", color="#2b6cb0")
    ax2b.set_ylim(0, 75)
    ax2.set_title("Sorted by survival (k_g): broadest ≠ longest-lived", fontsize=9)
    lines = ax2.containers + ax2b.containers
    ax2.legend(lines, [c.get_label() for c in lines], fontsize=7, loc="upper center")

    fig.suptitle(
        "Duration of response — depth is not durability: the highest response rate can be "
        "the least durable, which is the mechanism of the ORR→OS surrogate failure",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "duration_of_response.png", dpi=120)
    plt.close(fig)


def cross_context_generalization_figure() -> None:
    """The headline findings are not NSCLC artifacts: the ORR->OS surrogate inversion and
    the budget's survival-link axis reproduce across five solid-tumor contexts."""
    from onkos.budget import model_selection_budget
    from onkos.response import response_vs_survival

    ds = onkos.load()
    contexts = ["NSCLC", "breast", "CRC", "HCC", "melanoma"]
    t = np.linspace(0.0, 312.0, 625)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: ORR->OS discordant fraction, week-8 (faithful) vs k_g (inverted), per context.
    w_disc, k_disc = [], []
    for tt in contexts:
        ctx = {"tumor_type": tt, "line": "first"}
        kg = f"survival_link.{tt.lower()}_os_growth_rate"
        w_disc.append(response_vs_survival(ds, context=ctx, t=t, n=200).discordant_fraction)
        k_disc.append(response_vs_survival(ds, context=ctx, survival_link=kg, t=t, n=200)
                      .discordant_fraction)
    x = np.arange(len(contexts))
    ax1.bar(x - 0.2, w_disc, 0.4, color="#2b6cb0", label="week-8 link (shrinkage surrogate)")
    ax1.bar(x + 0.2, k_disc, 0.4, color="#c53030", label="k_g link (tail-sensitive)")
    ax1.set_xticks(x, contexts, fontsize=8)
    ax1.set_ylim(0, 1.05)
    ax1.set_ylabel("ORR→OS discordant fraction")
    ax1.set_title("ORR→OS surrogate inverts under k_g — in every context", fontsize=9)
    ax1.legend(fontsize=7, loc="upper left")

    # Right: model-selection budget — the survival-link axis, empty before v0.29, now real.
    vlink, vstruct = [], []
    for tt in contexts:
        b = model_selection_budget(ds, context={"tumor_type": tt, "line": "first"},
                                   endpoint="OS", n=120)
        vlink.append(b.fractions["survival_link"])
        vstruct.append(b.structural_fraction)
    ax2.bar(x, vstruct, 0.6, color="#e2e8f0", label="all structural choices")
    ax2.bar(x, vlink, 0.6, color="#2f855a", label="survival-link axis (new in v0.29)")
    ax2.axhline(0.5, ls=":", color="grey", lw=1)
    for i, v in enumerate(vlink):
        ax2.text(i, vstruct[i] + 0.02, f"{v * 100:.0f}%", ha="center", fontsize=8, color="#2f855a")
    ax2.set_xticks(x, contexts, fontsize=8)
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("share of forecast variance")
    ax2.set_title("Budget survival-link axis now populated across tumors", fontsize=9)
    ax2.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        "Cross-context generalization (v0.29): the resistance-mechanism, surrogate, and "
        "budget findings reproduce across five solid-tumor contexts — not NSCLC artifacts",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "cross_context_generalization.png", dpi=120)
    plt.close(fig)


def kill_mechanism_figure() -> None:
    """Norton-Simon (kill ∝ growth) vs log-kill (Claret, kill ∝ size) mechanisms."""
    from onkos.export.registry import get_kernel, kernel_values

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0.0, 156.0, 313)
    ns = onkos.simulate(ds, "drug_effect.norton_simon.nsclc", context=ctx, drug_effect=1.0, t=t)
    cl = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0, t=t)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    ax1.plot(t, ns.tumor_size, color=PALETTE[2], lw=1.8, label="Norton-Simon (kill ∝ growth)")
    ax1.plot(t, cl.tumor_size, color=PALETTE[1], lw=1.8, label="Claret log-kill + resistance")
    ax1.set_title("Same context, different kill mechanism")
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("tumor size (mm, SLD)")
    ax1.legend(fontsize=8)
    ax1.annotate("eradication\n(no resistance)", (110, 5), fontsize=8, color=PALETTE[2])
    ax1.annotate("resistance-driven\nregrowth", (110, cl.tumor_size[220] + 20), fontsize=8,
                 color=PALETTE[1])

    # The Norton-Simon signature: fractional kill rises as the tumor shrinks.
    spec = get_kernel(ds["drug_effect.norton_simon.nsclc"])
    v = kernel_values(ds["drug_effect.norton_simon.nsclc"])
    v["E"] = 1.0
    sizes = np.linspace(5, 195, 100)
    frac = np.array([-spec.rhs(0.0, [s], v)[0] / s for s in sizes])
    ax2.plot(sizes, frac, color=PALETTE[2])
    ax2.axhline(0, ls=":", color="grey", lw=1)
    ax2.set_title("Norton-Simon hypothesis: smaller tumor → more sensitive")
    ax2.set_xlabel("tumor size (mm)")
    ax2.set_ylabel("fractional kill rate (1/week)")

    fig.suptitle("Drug-effect subsystem (spec §3): the assumed KILL MECHANISM is itself a "
                 "model-selection axis", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "kill_mechanism.png", dpi=120)
    plt.close(fig)


def pfs_routes_figure() -> None:
    """PFS has two routes — the statistical week-8-keyed hazard link and the mechanistic
    RECIST time-to-progression — and they invert the model ranking for resistance dynamics,
    in every solid-tumor context. The week-8 link is blind to the regrowth tail."""
    from onkos.response import pfs_route_divergence, progression_free_survival

    ds = onkos.load()
    nsclc = {"tumor_type": "NSCLC", "line": "first"}
    names = {
        "resistance.claret_2009.tgi": "Claret",
        "drug_effect.norton_simon.nsclc": "Norton-Simon",
        "resistance.nsclc_first_line.two_population": "two-population",
        "tgi_metrics.wang_2009.biexponential": "Wang biexp",
    }
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # Left: NSCLC two-route bars per model, sorted by mechanistic TTP. The two-population
    # model is shortest mechanically yet longest statistically — the route inversion.
    div = pfs_route_divergence(ds, context=nsclc, n=300)
    rows = sorted(div.rows, key=lambda r: -(r["median_ttp_weeks"] or 0))
    labels = [names.get(r["record_id"], r["record_id"].split(".")[1]) for r in rows]
    x = np.arange(len(rows))
    mech = [r["median_ttp_weeks"] for r in rows]
    stat = [r["median_pfs_link_weeks"] for r in rows]
    ax1.bar(x - 0.2, mech, 0.4, color="#2f855a", label="mechanistic (RECIST TTP)")
    ax1.bar(x + 0.2, stat, 0.4, color="#c05621", label="statistical (week-8 link)")
    tp = next(i for i, r in enumerate(rows)
              if r["record_id"] == "resistance.nsclc_first_line.two_population")
    ax1.annotate("route inverts\nthis model's rank", (tp, stat[tp]),
                 textcoords="offset points", xytext=(0, 16), fontsize=7, ha="center",
                 color="#c53030", arrowprops=dict(arrowstyle="->", color="#c53030", lw=1))
    ax1.set_xticks(x, labels, fontsize=8)
    ax1.set_ylabel("median PFS (weeks)")
    ax1.set_title(f"NSCLC 1L: the two routes disagree "
                  f"({div.discordant_pairs}/{div.total_pairs} pairs inverted)", fontsize=9)
    ax1.legend(fontsize=7, loc="upper right")

    # Right: across all contexts, the Claret-vs-two-population route inversion. Mechanistic
    # ranks Claret >> two-pop; statistical ranks them level/reversed — everywhere.
    contexts = [("NSCLC", "first"), ("breast", "first"), ("CRC", "first"),
                ("HCC", "first"), ("melanoma", "first")]
    clab = [c[0] for c in contexts]
    cm, tm, cs, ts = [], [], [], []
    for tt, ln in contexts:
        ctx = {"tumor_type": tt, "line": ln}
        cl = progression_free_survival(ds, "resistance.claret_2009.tgi", context=ctx, n=300)
        tp_id = next(r.id for r in ds if r.id.endswith("two_population")
                     and r.derivation_context and r.derivation_context.tumor_type == tt)
        twp = progression_free_survival(ds, tp_id, context=ctx, n=300)
        cm.append(cl.median_ttp_weeks)
        tm.append(twp.median_ttp_weeks)
        cs.append(cl.median_pfs_link_weeks)
        ts.append(twp.median_pfs_link_weeks)
    xc = np.arange(len(contexts))
    ax2.plot(xc, cm, "o-", color="#2f855a", label="Claret — mechanistic")
    ax2.plot(xc, tm, "o--", color="#2f855a", alpha=0.5, label="two-pop — mechanistic")
    ax2.plot(xc, cs, "s-", color="#c05621", label="Claret — statistical")
    ax2.plot(xc, ts, "s--", color="#c05621", alpha=0.5, label="two-pop — statistical")
    ax2.set_xticks(xc, clab, fontsize=8)
    ax2.set_ylabel("median PFS (weeks)")
    ax2.set_title("Mechanistic ranks Claret ≫ two-pop; statistical levels them — every context",
                  fontsize=8.5)
    ax2.legend(fontsize=6.5, ncol=2, loc="upper right")

    fig.suptitle(
        "PFS endpoint — the route is a model-selection axis: the week-8 hazard link is blind "
        "to the resistant-clone regrowth the mechanism sees, so the two routes invert",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "pfs_routes.png", dpi=120)
    plt.close(fig)


def optimal_design_figure() -> None:
    """D-optimal trial design: the best schedule a fixed budget allows concentrates samples at
    the kill phase and regrowth onset, rescues the borderline resistance term, but cannot
    rescue the structurally flat growth rate — separating circumstantial from structural
    unidentifiability."""
    from onkos.design import optimal_schedule

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    od = optimal_schedule(ds, "resistance.claret_2009.tgi", context=ctx, n_samples=7, horizon=48.0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # Left: where the two designs sample, over the tumor trajectory. The optimal design
    # clusters at the kill phase and the regrowth onset.
    t = np.linspace(0, 48, 193)
    traj = onkos.simulate(ds, "resistance.claret_2009.tgi", context=ctx, drug_effect=1.0, t=t)
    ax1.plot(t, traj.tumor_size, color="#718096", lw=1.6, zorder=1)
    ax1.scatter(od.uniform.schedule,
                [np.interp(s, t, traj.tumor_size) for s in od.uniform.schedule],
                marker="v", s=90, color="#c05621", label="uniform", zorder=3)
    ax1.scatter(od.optimal.schedule,
                [np.interp(s, t, traj.tumor_size) for s in od.optimal.schedule],
                marker="o", s=90, color="#2f855a", edgecolor="black", label="D-optimal", zorder=4)
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("tumor size (mm, SLD)")
    ax1.set_title(f"Where the budget is spent (N={od.n_samples}, D-eff {od.d_efficiency:.2f}x)",
                  fontsize=9)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.annotate("kill phase", (10, np.interp(10, t, traj.tumor_size)),
                 textcoords="offset points", xytext=(4, -22), fontsize=7, color="#2f855a")
    ax1.annotate("regrowth onset", (33, np.interp(33, t, traj.tumor_size)),
                 textcoords="offset points", xytext=(-6, 12), fontsize=7, color="#2f855a")

    # Right: per-parameter RSE, uniform vs D-optimal, with the 50% identifiability line.
    syms = list(od.uniform.rse_percent)
    x = np.arange(len(syms))
    u = [min(od.uniform.rse_percent[s], 260) for s in syms]
    o = [min(od.optimal.rse_percent[s], 260) for s in syms]
    ax2.bar(x - 0.2, u, 0.4, color="#c05621", label="uniform")
    ax2.bar(x + 0.2, o, 0.4, color="#2f855a", label="D-optimal")
    ax2.axhline(od.rse_ceiling_percent, ls="--", color="#c53030", lw=1)
    ax2.text(len(syms) - 0.5, od.rse_ceiling_percent + 6, "identifiability line (50%)",
             fontsize=7, color="#c53030", ha="right")
    ax2.set_xticks(x, syms, fontsize=9)
    ax2.set_ylabel("predicted RSE (%)  [capped at 260]")
    ax2.set_title("λ rescued across the line; kL stays structurally flat", fontsize=9)
    ax2.legend(fontsize=8, loc="upper center")

    fig.suptitle(
        "D-optimal trial design — the best schedule a fixed budget allows rescues the "
        "circumstantially flat λ but cannot rescue the structurally flat kL (v0.22 capstone)",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "optimal_design.png", dpi=120)
    plt.close(fig)


def acquired_resistance_figure() -> None:
    """The resistance ORIGIN as a model-selection axis: acquired (drug-induced switching, R
    generated from zero) vs pre-existing (a baseline resistant subclone). Matched on kg/kd/kgr
    they agree at week 8 and on the week-8 OS surrogate, but the acquired model has a shallower
    nadir and earlier progression — a tail divergence the surrogate cannot see."""
    from onkos.export.reference import init_vector
    from onkos.export.registry import get_kernel, kernel_values
    from onkos.response import time_to_progression
    from scipy.integrate import solve_ivp

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    acq_id = "resistance.nsclc_first_line.acquired"
    pre_id = "resistance.nsclc_first_line.two_population"
    t = np.linspace(0.0, 156.0, 313)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Left: the acquired clone decomposition — the resistant pool is GENERATED from zero by
    # drug-induced switching (contrast two_population, where R starts at a pre-existing R0).
    spec = get_kernel(ds[acq_id])
    vals = kernel_values(ds[acq_id])
    vals["V0"] = 100.0
    vals["E"] = 1.0
    sol = solve_ivp(lambda tt, yy: spec.rhs(tt, yy, vals), (t[0], t[-1]),
                    init_vector(spec, vals), t_eval=t, rtol=1e-8, atol=1e-10, method="LSODA")
    s, r = sol.y
    ax1.semilogy(t, np.maximum(s, 1e-3), color="#2b6cb0", lw=1.6, label="sensitive S (killed + switching out)")
    ax1.semilogy(t, np.maximum(r, 1e-3), color="#c53030", lw=1.6, label="resistant R (generated from 0)")
    ax1.semilogy(t, s + r, color="black", lw=2.0, ls="--", label="observed tumor V = S + R")
    ax1.set_ylim(1e-2, 2e2)
    ax1.set_title("Acquired resistance — R is generated by treatment", fontsize=9)
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("tumor size (mm, SLD, log)")
    ax1.legend(fontsize=7, loc="lower right")
    ax1.annotate("no resistance\nat baseline (R0=0)", (2, 1e-1), fontsize=7, color="#c53030")

    # Right: the two origins' observed tumors — agree at week 8 / on OS, diverge in the tail.
    acq = onkos.simulate(ds, acq_id, context=ctx, drug_effect=1.0, t=t)
    pre = onkos.simulate(ds, pre_id, context=ctx, drug_effect=1.0, t=t)
    ttp_a = time_to_progression(t, acq.tumor_size)
    ttp_p = time_to_progression(t, pre.tumor_size)
    ax2.semilogy(t, pre.tumor_size, color="#c05621", lw=1.8,
                 label=f"pre-existing subclone  nadir {pre.tumor_size.min():.1f}, TTP {ttp_p:.0f}w")
    ax2.semilogy(t, acq.tumor_size, color="#2f855a", lw=1.8,
                 label=f"acquired switching  nadir {acq.tumor_size.min():.1f}, TTP {ttp_a:.0f}w")
    ax2.axvline(8, ls=":", color="grey", lw=1)
    ax2.text(9, 60, "week 8\n(agree → same OS)", fontsize=7, color="grey")
    ax2.set_title(f"Same week-8 OS ({acq.median_os:.0f} vs {pre.median_os:.0f} wk), different tail",
                  fontsize=9)
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("tumor size (mm, SLD, log)")
    ax2.legend(fontsize=7.5, loc="lower right")

    fig.suptitle(
        "Acquired vs pre-existing resistance: the resistance ORIGIN is a model-selection axis — "
        "matched kinetics, same week-8 OS, but a shallower nadir & earlier progression",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "acquired_resistance.png", dpi=120)
    plt.close(fig)


def burden_auc_figure() -> None:
    """The integrated tumor burden as a THIRD TGI→OS bridge metric. Left: the log relative
    tumor-size curve whose time-average the metric integrates (depth lowers it, the tail raises
    it). Right: median OS across the three bridge metrics — week-8 (depth-only), k_g (tail-only),
    burden (both) — re-rank the model set THREE distinct ways. The complete responder inverts from
    mid-pack to first; the minimal responder (Wang) is ranked 2nd by tail-only k_g but last by the
    depth-aware burden metric."""
    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0.0, 260.0, 521)
    models = [
        ("resistance.claret_2009.tgi", "Claret (phenom.)", "#c05621"),
        ("resistance.nsclc_first_line.two_population", "two-population", "#2f855a"),
        ("drug_effect.norton_simon.nsclc", "Norton-Simon (complete resp.)", "#2b6cb0"),
        ("tgi_metrics.wang_2009.biexponential", "Wang biexp (minimal resp.)", "#6b46c1"),
        ("resistance.nsclc_first_line.acquired", "acquired resistance", "#b7791f"),
    ]
    links = [
        ("week-8\n(depth-only)", "survival_link.nsclc_os_week8"),
        ("k_g\n(tail-only)", "survival_link.nsclc_os_growth_rate"),
        ("burden\n(depth+tail)", "survival_link.nsclc_os_burden_auc"),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Left: the log relative tumor-size curve the burden metric integrates.
    for rid, label, color in models:
        tr = onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t)
        rel = np.maximum(tr.tumor_size / tr.tumor_size[0], 1e-3)
        b = tr.metrics["log_burden_auc"]
        ax1.plot(t, np.log(rel), color=color, lw=1.7, label=f"{label}  (burden {b:+.1f})")
    ax1.axhline(0.0, ls=":", color="grey", lw=1)
    ax1.text(5, 0.25, "baseline (burden 0)", fontsize=7, color="grey")
    ax1.set_title("log relative tumor size — the burden metric is its time-average", fontsize=9)
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("log(v / y0)")
    ax1.set_ylim(-7.2, 5.5)
    ax1.legend(fontsize=7, loc="upper left")

    # Right: median OS slopegraph across the three bridge metrics.
    os_by_link = {
        rid: [onkos.simulate(ds, rid, context=ctx, drug_effect=1.0, t=t, survival_link=link).median_os
              for _, link in links]
        for rid, _, _ in models
    }
    xs = np.arange(len(links))
    for rid, label, color in models:
        ys = [o if o else np.nan for o in os_by_link[rid]]
        ax2.plot(xs, ys, color=color, lw=1.8, marker="o", ms=5, label=label)
        ax2.annotate(f"{ys[-1]:.0f}", (xs[-1], ys[-1]), xytext=(6, 0),
                     textcoords="offset points", fontsize=7, color=color, va="center")
    ax2.set_xticks(xs)
    ax2.set_xticklabels([name for name, _ in links], fontsize=8)
    ax2.set_xlim(-0.3, len(links) - 0.3)
    ax2.set_ylabel("median OS (weeks)")
    ax2.set_title("each metric re-ranks the models a different way", fontsize=9)
    ax2.legend(fontsize=7, loc="upper left")

    fig.suptitle(
        "The integrated tumor burden is a third TGI→OS bridge metric: it sees both depth and tail, "
        "so it ranks the models differently from week-8 AND from k_g",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "burden_auc.png", dpi=120)
    plt.close(fig)


def joint_survival_figure() -> None:
    """Joint (current-value) vs two-stage survival. Left: the hazard ratio is TIME-VARYING —
    suppressed while the tumor is small, then rising steeply as a resistant clone regrows — a
    non-proportional hazard the two-stage links (a constant HR) cannot represent. Right: the
    resulting OS curves; the joint link bends down in the tail for the regrowing models, and
    it inverts the week-8 ranking of the light- vs heavy-tail resistance models."""
    from onkos.joint import joint_survival

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    t = np.linspace(0.0, 260.0, 521)
    models = [
        ("resistance.claret_2009.tgi", "Claret (phenom.)", "#c05621"),
        ("resistance.nsclc_first_line.two_population", "two-population", "#2f855a"),
        ("drug_effect.norton_simon.nsclc", "Norton-Simon (complete resp.)", "#2b6cb0"),
        ("resistance.nsclc_first_line.acquired", "acquired resistance", "#b7791f"),
    ]
    js = {rid: joint_survival(ds, rid, context=ctx, drug_effect=1.0, t=t, alpha=1.0)
          for rid, _, _ in models}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Left: the time-varying hazard ratio (the non-proportionality).
    for rid, label, color in models:
        ax1.semilogy(t, np.maximum(js[rid].hazard_ratio, 1e-3), color=color, lw=1.7,
                     label=f"{label}  (PH-viol {js[rid].ph_violation:.0f}x)")
    ax1.axhline(1.0, ls=":", color="grey", lw=1)
    ax1.text(5, 1.25, "HR = 1 (baseline size)\n— a two-stage link is FLAT here", fontsize=7, color="grey")
    ax1.set_title("the current-value hazard ratio is time-varying (non-proportional)", fontsize=9)
    ax1.set_xlabel("weeks")
    ax1.set_ylabel("hazard ratio HR(t) = exp(α·log(v/y0))")
    ax1.set_ylim(1e-2, 3e2)
    ax1.legend(fontsize=7, loc="upper left")

    # Right: joint (solid) vs two-stage week-8 (dashed) OS for the light- vs heavy-tail pair.
    for rid, label, color in models[:2]:
        j = js[rid]
        jm = f"{j.median_os:.0f}" if j.median_os else "n/r"
        tm = f"{j.two_stage_median_os:.0f}" if j.two_stage_median_os else "n/r"
        ax2.plot(t, j.os_curve, color=color, lw=1.9, label=f"{label} — joint (mOS {jm})")
        ax2.plot(t, j.two_stage_curve, color=color, lw=1.3, ls="--",
                 label=f"{label} — two-stage week-8 (mOS {tm})")
    ax2.axhline(0.5, ls=":", color="grey", lw=1)
    ax2.set_title("week-8 ranks two-pop > Claret; the joint link inverts it", fontsize=9)
    ax2.set_xlabel("weeks")
    ax2.set_ylabel("survival fraction")
    ax2.set_ylim(0, 1.02)
    ax2.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        "Joint (current-value) vs two-stage survival: the hazard tracks the tumor in real time, so "
        "a regrowing clone makes the hazard ratio rise — a non-proportional hazard PH links can't encode",
        fontsize=9.0,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "joint_survival.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    divergence_figure()
    tier_figure()
    exposure_response_figure()
    context_library_figure()
    preclinical_figure()
    immuno_oncology_figure()
    coverage_figure()
    uncertainty_figure()
    tgi_metrics_figure()
    sensitivity_figure()
    survival_endpoints_figure()
    survival_model_choice_figure()
    line_of_therapy_figure()
    composability_chain_figure()
    kill_mechanism_figure()
    model_average_figure()
    identifiability_figure()
    combination_interaction_figure()
    two_population_resistance_figure()
    survival_metric_choice_figure()
    model_selection_budget_figure()
    response_orr_figure()
    duration_of_response_figure()
    cross_context_generalization_figure()
    pfs_routes_figure()
    optimal_design_figure()
    acquired_resistance_figure()
    burden_auc_figure()
    joint_survival_figure()
    print(f"Wrote figures to {OUT}")
