"""Onkos Streamlit dashboard — browse + virtual-trial divergence view.

Run: ``streamlit run dashboard/app.py``

Population/trial-level forward simulation only. NOT a prognostic tool, NOT a
treatment recommender. The divergence view answers "how do published models
disagree at the trial level?" — never "your tumor / your survival."
"""

from __future__ import annotations

import numpy as np
import streamlit as st

import onkos
from onkos._const import CLINICAL_USE

st.set_page_config(page_title="Onkos — virtual-trial divergence", layout="wide")


@st.cache_data
def _load():
    ds = onkos.load()
    return ds


ds = _load()

st.title("Onkos — virtual-trial divergence view")
st.error(f"NOT FOR CLINICAL USE — clinicalUse = {CLINICAL_USE}")

tumor_types = sorted({r.derivation_context.tumor_type for r in ds
                      if r.derivation_context and r.derivation_context.tumor_type})

with st.sidebar:
    st.header("Context")
    tumor_type = st.selectbox("Tumor type", tumor_types, index=tumor_types.index("NSCLC")
                              if "NSCLC" in tumor_types else 0)
    line = st.selectbox("Line of therapy", ["first", "second"], index=0)
    drug_effect = st.slider("Drug-effect size (E)", 0.0, 2.0, 1.0, 0.1)
    horizon = st.slider("Horizon (weeks)", 26, 156, 104, 13)

t = np.linspace(0.0, float(horizon), 2 * horizon + 1)
ctx = {"tumor_type": tumor_type, "line": line}
cmp = onkos.compare(ds, purpose="tgi", context=ctx, drug_effect=drug_effect, t=t)

tab_div, tab_browse = st.tabs(["Divergence view", "Browse dataset"])

with tab_div:
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Tumor-size trajectories")
        ts = {tr.record_id: tr.tumor_size for tr in cmp.included}
        if ts:
            st.line_chart({k: v for k, v in ts.items()}, x=None)
    with c2:
        st.subheader("Population OS curves")
        os_ = {tr.record_id: tr.os_curve for tr in cmp.included if tr.os_curve is not None}
        if os_:
            st.line_chart(os_)

    m1, m2, m3 = st.columns(3)
    m1.metric("Eligible models", len(cmp.included))
    m2.metric("OS divergence (max pointwise)", f"{cmp.os_divergence:.3f}")
    rng = cmp.median_os_range
    m3.metric("Median-OS spread (wk)", f"{(rng[1] - rng[0]):.1f}" if rng else "n/a")

    st.subheader("Included models")
    st.table([
        {
            "id": tr.record_id,
            "tier": tr.tier,
            "week8_change": round(tr.metrics["week8_relative_change"], 3),
            "depth_of_response": round(tr.metrics["depth_of_response"], 3),
            "median_OS_wk": round(tr.median_os, 1) if tr.median_os else None,
        }
        for tr in cmp.included
    ])

    if cmp.excluded:
        st.subheader("Excluded (out-of-context transport — greyed out)")
        st.table([{"id": rid, "reason": reason} for rid, reason in cmp.excluded])

with tab_browse:
    st.subheader(f"{len(ds)} records")
    st.table([
        {
            "id": r.id,
            "subsystem": r.subsystem,
            "purpose": r.purpose,
            "tier": r.tier,
            "review": r.review_status,
            "tumor_type": r.derivation_context.tumor_type if r.derivation_context else None,
            "citation": r.primary_citation.key if r.primary_citation else None,
        }
        for r in ds
    ])
