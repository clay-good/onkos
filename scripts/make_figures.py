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

    excluded = ", ".join(rid.split(".")[1] for rid, _ in cmp.excluded)
    fig.suptitle(
        f"Virtual-trial divergence — model choice moves median OS across "
        f"{cmp.median_os_range[0]:.0f}-{cmp.median_os_range[1]:.0f} wk  "
        f"(greyed out for out-of-context transport: {excluded})",
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


if __name__ == "__main__":
    divergence_figure()
    tier_figure()
    exposure_response_figure()
    print(f"Wrote figures to {OUT}")
