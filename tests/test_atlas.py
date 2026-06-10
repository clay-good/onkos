"""Model-selection atlas — the synthesis index + per-context survey (v0.39.0).

The atlas is the navigational layer over the model-selection axes: a declarative registry
(AXES) and a one-call per-context survey that reports each applicable axis's NATIVE headline.
It is deliberately a survey, not a decomposition — different units, not orthogonal — so it
flags comparable=False and points to the budget for the rigorous partition.

These landmarks check the registry is well-formed, the survey runs and agrees with the
underlying modules, and the honesty guardrails (comparable flag, budget pointer, tier) hold.
"""

import numpy as np
import onkos
from onkos.atlas import AXES, Atlas, Axis, model_selection_atlas

NSCLC = {"tumor_type": "NSCLC", "line": "first"}


# ---- registry landmarks ----------------------------------------------------

def test_registry_is_well_formed():
    """Every axis has all fields, unique keys, a valid scope, and a real module path."""
    assert len(AXES) >= 7
    keys = [a.key for a in AXES]
    assert len(keys) == len(set(keys))  # unique
    for a in AXES:
        assert isinstance(a, Axis)
        assert a.key and a.label and a.varies and a.finding
        assert a.module.startswith("onkos.")
        assert a.scope in ("single-agent", "combination")
        assert a.version.startswith("v0.")


def test_registry_covers_the_shipped_axes():
    """The registry names the axes shipped across the project (catches an un-registered
    new axis or a renamed module)."""
    keys = {a.key for a in AXES}
    expected = {"tgi_model", "survival_link", "survival_structure", "exposure_response",
                "readout_timing", "model_discriminability", "additivity_reference"}
    assert expected <= keys
    # each named module is importable
    import importlib
    for a in AXES:
        importlib.import_module(a.module)


# ---- survey landmarks ------------------------------------------------------

def test_atlas_surveys_the_single_agent_axes():
    a = model_selection_atlas(onkos.load(), context=NSCLC)
    assert isinstance(a, Atlas)
    surveyed = {e.key for e in a.entries}
    # the single-agent axes are run; the combination axis (needs a regimen) is not
    assert {"tgi_model", "survival_link", "survival_structure", "exposure_response",
            "readout_timing", "model_discriminability"} <= surveyed
    assert "additivity_reference" not in surveyed


def test_atlas_headlines_agree_with_the_underlying_modules():
    """The survey is a thin, faithful orchestration — its headlines equal the modules' when
    given the same horizon the atlas uses (t = linspace(0, 312, 625))."""
    ds = onkos.load()
    t = np.linspace(0.0, 312.0, 625)
    a = model_selection_atlas(ds, context=NSCLC)

    from onkos.discriminability import model_discriminability
    md = model_discriminability(ds, context=NSCLC, t=t)
    assert a.get("model_discriminability").headline == md.n_indistinguishable

    from onkos.early_surrogate import surrogate_timing_fidelity
    st = surrogate_timing_fidelity(ds, context=NSCLC, t=t)
    assert a.get("readout_timing").headline == st.earliest_discordance

    from onkos.compare import compare
    rng = compare(ds, purpose="tgi", context=NSCLC, t=t).median_os_range
    assert np.isclose(a.get("tgi_model").headline, rng[1] - rng[0], atol=0.11)


def test_os_swing_axes_are_grouped_and_positive():
    """The weeks-unit axes are the loosely-comparable OS-swing group; each is a real spread."""
    a = model_selection_atlas(onkos.load(), context=NSCLC)
    swing = a.os_swing_axes
    assert {e.key for e in swing} == {"tgi_model", "survival_link", "survival_structure",
                                      "exposure_response"}
    assert all(e.headline is None or e.headline >= 0 for e in swing)
    # the survival-structure / metric choices are large swings here
    assert a.get("survival_structure").headline > 20.0


def test_atlas_is_a_survey_not_a_decomposition():
    """The honesty guardrails: comparable=False and an explicit pointer to the budget for the
    rigorous orthogonal partition — so the atlas is never read as a variance decomposition."""
    a = model_selection_atlas(onkos.load(), context=NSCLC)
    assert a.comparable is False
    assert "budget" in a.note.lower()


def test_atlas_inherits_tier_and_carries_clinical_use():
    a = model_selection_atlas(onkos.load(), context=NSCLC)
    assert a.tier in ("A", "B", "C", "D")
    import json

    d = json.loads(a.to_json())
    assert d["onkos:clinicalUse"].startswith("PROHIBITED")
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert d["comparable"] is False


def test_atlas_reproduces_across_contexts():
    """The survey runs for the other solid-tumor contexts too (they have fewer models, so
    fewer/smaller axes, but the survey is well-formed)."""
    ds = onkos.load()
    for tt in ("breast", "CRC"):
        a = model_selection_atlas(ds, context={"tumor_type": tt, "line": "first"})
        assert len(a.entries) >= 4
        assert a.get("model_discriminability") is not None
