"""Dashboard smoke test: the Streamlit app renders every curated context without raising.

Runs headlessly through Streamlit's ``AppTest``; skipped when the optional ``dashboard``
extra is not installed (CI installs it in the exports job).
"""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "dashboard" / "app.py")


def _select(app, label):
    return next(s for s in app.selectbox if s.label == label)


def _render(tumor_type=None, line="first", drug_effect=1.0):
    app = AppTest.from_file(APP, default_timeout=120).run()
    if tumor_type is not None:
        _select(app, "Tumor type").set_value(tumor_type)
    _select(app, "Line of therapy").set_value(line)
    next(s for s in app.slider if s.label.startswith("Drug-effect")).set_value(drug_effect)
    return app.run()


def _tumor_types():
    return list(_select(AppTest.from_file(APP, default_timeout=120).run(), "Tumor type").options)


def test_dashboard_renders_every_tumor_type():
    tumor_types = _tumor_types()
    assert {"NSCLC", "breast", "CRC", "HCC", "melanoma"} <= set(tumor_types)
    for tumor_type in tumor_types:
        app = _render(tumor_type)
        assert not app.exception, (tumor_type, [e.value for e in app.exception])
        assert any("NOT FOR CLINICAL USE" in e.value for e in app.error)  # the safety banner


@pytest.mark.parametrize("drug_effect", [0.0, 2.0])
def test_dashboard_handles_the_effect_slider_extremes(drug_effect):
    assert not _render("NSCLC", drug_effect=drug_effect).exception


def test_dashboard_reports_a_context_with_no_models_instead_of_failing():
    app = _render("breast", line="second")  # no second-line breast models are curated
    assert not app.exception
    assert any("No eligible models" in w.value for w in app.warning)
