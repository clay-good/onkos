"""Phase F — dataset health report and reproducibility gate."""

from pathlib import Path

import onkos
import pytest
from onkos._data import dataset_dir
from onkos.report import build_report, external_validation_coverage, stats


def test_external_validation_fully_covered():
    ds = onkos.load()
    n_val, n_elig, frac = external_validation_coverage(ds)
    assert n_elig >= 15
    assert n_val == n_elig  # every eligible clinical model has an external metric
    assert frac == 1.0


def test_stats_shape():
    s = stats(onkos.load())
    assert s["by_tier"]["D"] >= 2  # the hypothesis-tier IO records
    assert set(s["by_review_status"]) <= {"unverified", "verified", "contested"}
    assert "immuno_oncology" in s["tier_by_subsystem"]
    # IO is entirely tier D
    assert set(s["tier_by_subsystem"]["immuno_oncology"]) == {"D"}


def test_report_contains_key_sections():
    md = build_report(onkos.load())
    for heading in (
        "# Onkos dataset health report",
        "## Confidence tiers",
        "## Tier × subsystem",
        "## External-validation backlog",
        "NOT FOR PREDICTION",
    ):
        assert heading in md


def _committed_report() -> Path:
    return dataset_dir().parent / "docs" / "dataset-health.md"


@pytest.mark.skipif(not _committed_report().exists(), reason="source checkout only")
def test_committed_report_is_in_sync():
    """The committed report must match a freshly generated one (reproducibility)."""
    md = build_report(onkos.load())
    committed = _committed_report().read_text()
    assert md == committed, "docs/dataset-health.md is stale — run `onkos report --output docs/dataset-health.md`"
