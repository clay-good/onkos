"""COMBINE archives must be reproducible independently of the build clock."""

import zipfile

import onkos
import pytest
from onkos.export.combine import build_omex


@pytest.mark.parametrize("record_id", [
    "resistance.claret_2009.tgi",
    "survival_link.nsclc_os_week8",
])
def test_archive_is_independent_of_build_time(tmp_path, monkeypatch, record_id):
    record = onkos.load()[record_id]
    monkeypatch.setattr(zipfile.time, "localtime", lambda *args: (2024, 1, 1, 0, 0, 0, 0, 1, -1))
    first = build_omex(record, tmp_path / "first.omex")
    monkeypatch.setattr(zipfile.time, "localtime", lambda *args: (2026, 9, 28, 12, 0, 0, 0, 1, -1))
    second = build_omex(record, tmp_path / "second.omex")
    assert first.read_bytes() == second.read_bytes()
    with zipfile.ZipFile(first) as archive:
        assert archive.testzip() is None
        assert "manifest.xml" in archive.namelist()
        for entry in archive.infolist():
            assert entry.compress_type == zipfile.ZIP_DEFLATED
            assert archive.read(entry)
