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
    print(f"Wrote figures to {OUT}")
