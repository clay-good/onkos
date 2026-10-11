"""CLI smoke tests."""

from onkos.cli import main


def test_version(capsys):
    assert main(["version"]) == 0
    assert "onkos" in capsys.readouterr().out


def test_validate(capsys):
    assert main(["validate"]) == 0
    assert "OK" in capsys.readouterr().out


def test_info(capsys):
    assert main(["info"]) == 0
    out = capsys.readouterr().out
    assert "By tier" in out and "By subsystem" in out


def test_report(capsys):
    assert main(["report"]) == 0
    out = capsys.readouterr().out
    assert "dataset health report" in out
    assert "External-validation coverage" in out


def test_report_to_file(tmp_path):
    out = tmp_path / "health.md"
    assert main(["report", "--output", str(out)]) == 0
    assert "Confidence tiers" in out.read_text()


def test_uncertainty(capsys):
    assert main(["uncertainty", "resistance.claret_2009.tgi", "--n", "40"]) == 0
    out = capsys.readouterr().out
    assert "bands" in out and "median_os_weeks" in out


def test_sensitivity(capsys):
    assert main(["sensitivity", "resistance.claret_2009.tgi", "--n", "120"]) == 0
    out = capsys.readouterr().out
    assert "contribution" in out and "verify" in out


def test_simulate_compare(capsys):
    assert main(["simulate", "--compare", "--tumor-type", "NSCLC", "--line", "first"]) == 0
    out = capsys.readouterr().out
    assert "OS  divergence" in out and "PFS divergence" in out
    assert "EXCLUDED" in out


def test_simulate_compare_json(capsys):
    import json

    assert main(["simulate", "--compare", "--tumor-type", "NSCLC", "--json"]) == 0
    d = json.loads(capsys.readouterr().out)
    assert d["NOT_FOR_CLINICAL_USE"] is True
    assert "os_divergence" in d and "pfs_divergence" in d and d["included"]


def test_export_sbml(tmp_path, capsys):
    out = tmp_path / "sbml"
    assert main(["export", "--format", "sbml", "--output", str(out)]) == 0
    files = list(out.glob("*.xml"))
    assert files
    assert "clinicalUse" in files[0].read_text()


def test_export_omex(tmp_path):
    out = tmp_path / "omex"
    assert main(["export", "--format", "omex", "--output", str(out)]) == 0
    assert list(out.glob("*.omex"))


def test_simulate_with_exposure(capsys):
    rc = main([
        "simulate", "resistance.claret_2009.tgi",
        "--tumor-type", "NSCLC", "--line", "first",
        "--exposure", "200",
        "--exposure-response", "exposure_response.dacomitinib_egfr.emax",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "exposure=200.0 via exposure_response.dacomitinib_egfr.emax" in out


def test_export_pharmml_includes_er(tmp_path):
    out = tmp_path / "pharmml"
    assert main(["export", "--format", "pharmml", "--output", str(out)]) == 0
    assert (out / "exposure_response.emax_generic.pharmml").exists()


def test_export_pharmml_so(tmp_path):
    out = tmp_path / "so"
    assert main(["export", "--format", "so", "--output", str(out)]) == 0
    f = out / "resistance.claret_2009.tgi.so.xml"
    assert f.exists()
    text = f.read_text()
    assert "PopulationEstimates" in text and "clinicalUse" in text


def test_export_jsonld(tmp_path):
    out = tmp_path / "jsonld"
    assert main(["export", "--format", "jsonld", "--output", str(out)]) == 0
    f = out / "resistance.claret_2009.tgi.jsonld"
    assert f.exists()
    import json

    doc = json.loads(f.read_text())
    assert "@context" in doc and doc["@id"].endswith("resistance.claret_2009.tgi")


def test_list_shows_record_ids_and_filters(capsys):
    assert main(["list"]) == 0
    everything = capsys.readouterr().out
    assert "resistance.claret_2009.tgi" in everything and "growth_laws.gompertz" in everything

    assert main(["list", "--tumor-type", "NSCLC", "--purpose", "tgi"]) == 0
    filtered = capsys.readouterr().out
    assert "resistance.claret_2009.tgi" in filtered and "NSCLC/second" in filtered
    assert "breast" not in filtered and "growth_laws" not in filtered


def test_list_json_matches_the_dataset(capsys):
    import json

    import onkos

    assert main(["list", "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert [r["id"] for r in rows] == onkos.load().ids()
    assert {"purpose", "kernel", "tier", "review_status", "tumor_type"} <= set(rows[0])


def test_list_with_no_match_exits_nonzero(capsys):
    assert main(["list", "--tier", "A"]) == 1  # nothing is externally validated yet
    assert "No records match" in capsys.readouterr().err


def test_review_queue_links_each_record_to_its_source(capsys):
    import json

    assert main(["review-queue"]) == 0
    assert "source: https://doi.org/" in capsys.readouterr().out
    assert main(["review-queue", "--json"]) == 0
    queue = json.loads(capsys.readouterr().out)
    assert queue and all(r["source_links"] for r in queue)
