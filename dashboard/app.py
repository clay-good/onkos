"""Onkos Streamlit dashboard — the virtual-trial divergence view, plus the
parameter-uncertainty and sensitivity lenses, over the package API.

Run: ``streamlit run dashboard/app.py``  (needs the ``dashboard`` extra)

Population/trial-level forward simulation only. NOT a prognostic tool, NOT a
treatment recommender. The divergence view answers "how do published models
disagree at the trial level?" — never "your tumor / your survival." All data
comes from the tested package API (`onkos.compare`, `simulate_ensemble`,
`sensitivity`); this file is a thin presentation layer.
"""

from __future__ import annotations

import numpy as np
import onkos
import streamlit as st
from onkos._const import CLINICAL_USE

st.set_page_config(page_title="Onkos — virtual-trial divergence", layout="wide")


@st.cache_data
def _load():
    return onkos.load()


ds = _load()
st.title("Onkos — virtual-trial divergence view")
st.error(f"NOT FOR CLINICAL USE — clinicalUse = {CLINICAL_USE}")

tumor_types = sorted(
    {r.derivation_context.tumor_type for r in ds
     if r.derivation_context and r.derivation_context.tumor_type
     and r.purpose == "tgi" and r.subsystem not in ("preclinical_translation", "immuno_oncology")}
)

with st.sidebar:
    st.header("Context")
    tumor_type = st.selectbox("Tumor type", tumor_types,
                              index=tumor_types.index("NSCLC") if "NSCLC" in tumor_types else 0)
    line = st.selectbox("Line of therapy", ["first", "second"], index=0)
    drug_effect = st.slider("Drug-effect size (E)", 0.0, 2.0, 1.0, 0.1)
    horizon = st.slider("Horizon (weeks)", 26, 312, 156, 13)

t = np.linspace(0.0, float(horizon), 2 * horizon + 1)
ctx = {"tumor_type": tumor_type, "line": line}
cmp = onkos.compare(ds, purpose="tgi", context=ctx, drug_effect=drug_effect, t=t)

tab_div, tab_model, tab_browse = st.tabs(
    ["Divergence view", "Analyze a model", "Browse dataset"]
)

with tab_div:
    if not cmp.included:
        st.warning("No eligible models for this context.")
    else:
        c1, c2, c3 = st.columns(3)
        c1.subheader("Tumor size")
        c1.line_chart({tr.record_id.split(".", 1)[1]: tr.tumor_size for tr in cmp.included})
        c2.subheader("Population OS")
        c2.line_chart({tr.record_id.split(".", 1)[1]: tr.os_curve for tr in cmp.included
                       if tr.os_curve is not None})
        c3.subheader("Population PFS")
        c3.line_chart({tr.record_id.split(".", 1)[1]: tr.pfs_curve for tr in cmp.included
                       if tr.pfs_curve is not None})

        m = st.columns(4)
        m[0].metric("Eligible models", len(cmp.included))
        m[1].metric("OS divergence", f"{cmp.os_divergence:.3f}")
        m[2].metric("PFS divergence", f"{cmp.pfs_divergence:.3f}")
        rng = cmp.median_os_range
        m[3].metric("Median-OS spread (wk)", f"{(rng[1] - rng[0]):.1f}" if rng else "n/a")

        st.subheader("Included models")
        st.table([
            {
                "id": tr.record_id, "tier": tr.tier,
                "week8 Δ": round(tr.metrics["week8_relative_change"], 3),
                "depth": round(tr.metrics["depth_of_response"], 3),
                "k_g": round(tr.metrics["tumor_growth_rate_kg"], 4)
                if np.isfinite(tr.metrics["tumor_growth_rate_kg"]) else None,
                "median OS": round(tr.median_os, 1) if tr.median_os else None,
                "median PFS": round(tr.median_pfs, 1) if tr.median_pfs else None,
            }
            for tr in cmp.included
        ])
        if cmp.excluded:
            st.subheader("Excluded (out-of-context transport — greyed out)")
            st.table([{"id": rid, "reason": reason} for rid, reason in cmp.excluded])

with tab_model:
    options = [tr.record_id for tr in cmp.included]
    if not options:
        st.info("No eligible models to analyze for this context.")
    else:
        rid = st.selectbox("Model", options)
        col = st.columns(2)
        n = col[0].slider("Monte-Carlo samples", 50, 500, 200, 50)
        target = col[1].selectbox("Sensitivity target",
                                  ["median_os_weeks", "median_pfs_weeks",
                                   "week8_relative_change", "depth_of_response"])

        ens = onkos.simulate_ensemble(ds, rid, context=ctx, drug_effect=drug_effect, t=t, n=n)
        b, o = ens.tumor_size, ens.os_curve
        u1, u2 = st.columns(2)
        u1.subheader("Tumor size — uncertainty band")
        u1.line_chart({"lo": b.lo, "median": b.median, "hi": b.hi})
        if o is not None:
            u2.subheader("Population OS — uncertainty band")
            u2.line_chart({"lo": o.lo, "median": o.median, "hi": o.hi})

        try:
            res = onkos.sensitivity(ds, rid, context=ctx, target=target, n=max(n, 200))
            st.subheader(f"Sensitivity — what drives {target} (R²={res.r_squared:.2f})")
            st.bar_chart({p.symbol: p.contribution for p in res.indices})
            if res.dominant:
                st.caption(f"Verify **{res.dominant.symbol}** first "
                           f"({res.dominant.contribution * 100:.0f}% of variance).")
        except ValueError:
            st.caption("This model has no inter-individual variability to analyze.")

with tab_browse:
    st.subheader(f"{len(ds)} records")
    st.table([
        {
            "id": r.id, "subsystem": r.subsystem, "purpose": r.purpose, "tier": r.tier,
            "review": r.review_status,
            "tumor_type": r.derivation_context.tumor_type if r.derivation_context else None,
            "citation": r.primary_citation.key if r.primary_citation else None,
        }
        for r in ds
    ])
