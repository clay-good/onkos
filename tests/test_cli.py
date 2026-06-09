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


def test_simulate_compare(capsys):
    assert main(["simulate", "--compare", "--tumor-type", "NSCLC", "--line", "first"]) == 0
    out = capsys.readouterr().out
    assert "OS divergence" in out
    assert "EXCLUDED" in out


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
