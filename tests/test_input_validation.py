"""Input-validation contract: bad input fails loudly and legibly.

A typo'd record id, an out-of-domain number, or a context with no curated model must
never produce a raw traceback, a table of ``nan``, or — worst — a confident conclusion
drawn from zero models. The library raises a descriptive ``ValueError``/``KeyError``;
the CLI turns that into a one-line ``onkos: error:`` message and exit status 2.
"""

import numpy as np
import onkos
import pytest
from onkos.cli import main
from onkos.design import optimal_schedule
from onkos.discriminability import required_events
from onkos.identify import identifiability

CLARET = "resistance.claret_2009.tgi"
NSCLC = {"tumor_type": "NSCLC", "line": "first"}


@pytest.fixture(scope="module")
def ds():
    return onkos.load()


# --- library -----------------------------------------------------------------


def test_unknown_record_id_suggests_close_matches(ds):
    with pytest.raises(KeyError, match="unknown record id") as exc:
        ds["resistance.claret_2009.tgj"]
    assert CLARET in exc.value.args[0]


def test_simulate_rejects_non_dynamic_records(ds):
    with pytest.raises(ValueError, match="'survival' kernel"):
        onkos.simulate(ds, "survival_link.nsclc_os_growth_rate", context=NSCLC)


@pytest.mark.parametrize("effect", [float("nan"), float("inf")])
def test_simulate_rejects_non_finite_effect(ds, effect):
    with pytest.raises(ValueError, match="finite"):
        onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=effect)


@pytest.mark.parametrize("n", [0, -3])
def test_ensemble_needs_at_least_one_sample(ds, n):
    with pytest.raises(ValueError, match="at least one sample"):
        onkos.simulate_ensemble(ds, CLARET, context=NSCLC, n=n)


def test_sensitivity_rejects_unknown_target(ds):
    with pytest.raises(ValueError, match="unknown target 'bogus'"):
        onkos.sensitivity(ds, CLARET, context=NSCLC, target="bogus", n=5)


@pytest.mark.parametrize("sigma", [(0.0, 0.0), (-1.0, 0.0), (0.2, -1.0), (float("nan"), 0.0)])
def test_design_tools_need_a_real_residual_error_model(ds, sigma):
    prop, add = sigma
    with pytest.raises(ValueError, match="residual error"):
        identifiability(ds, CLARET, context=NSCLC, sigma_prop=prop, sigma_add=add)
    with pytest.raises(ValueError, match="residual error"):
        optimal_schedule(ds, CLARET, context=NSCLC, sigma_prop=prop, sigma_add=add)


def test_identifiability_rejects_negative_scan_times(ds):
    with pytest.raises(ValueError, match="non-negative"):
        identifiability(ds, CLARET, context=NSCLC, schedule=[-4.0, 0.0, 8.0])


def test_identifiability_does_not_depend_on_scan_order(ds):
    a = identifiability(ds, CLARET, context=NSCLC, schedule=[0, 6, 12, 24, 48])
    b = identifiability(ds, CLARET, context=NSCLC, schedule=[48, 0, 24, 6, 12])
    assert [p.rse_percent for p in a.params] == [p.rse_percent for p in b.params]


@pytest.mark.parametrize("horizon", [0.0, -5.0, float("inf")])
def test_optimal_schedule_needs_a_positive_horizon(ds, horizon):
    with pytest.raises(ValueError, match="horizon"):
        optimal_schedule(ds, CLARET, context=NSCLC, horizon=horizon)


@pytest.mark.parametrize("kw", [{"power": 0.0}, {"power": 1.0}, {"alpha": 0.0}, {"alpha": 2.0}])
def test_required_events_needs_power_and_alpha_in_unit_interval(kw):
    with pytest.raises(ValueError, match=r"\(0, 1\)"):
        required_events(0.5, **kw)


def test_valid_boundary_inputs_still_work(ds):
    """The guards reject only the meaningless: an untreated arm (E = 0), a single-sample
    ensemble, and an additive-only error model all remain legal."""
    assert np.isfinite(onkos.simulate(ds, CLARET, context=NSCLC, drug_effect=0.0).tumor_size).all()
    assert onkos.simulate_ensemble(ds, CLARET, context=NSCLC, n=1).n == 1
    assert identifiability(ds, CLARET, context=NSCLC, sigma_prop=0.0, sigma_add=2.0).params


def test_pure_cores_reject_malformed_inputs():
    """The landmark-tested cores are public too: mismatched or out-of-domain arguments
    must not broadcast or clamp their way to a plausible-looking number."""
    from onkos.combine import decompose
    from onkos.early_surrogate import landmark_response
    from onkos.joint import current_value_survival

    with pytest.raises(ValueError, match="same length"):
        decompose([60.0, 90.0], [4.0], [0.5, 0.5])  # used to broadcast the lone variance

    t = np.linspace(0.0, 52.0, 53)
    with pytest.raises(ValueError, match="non-negative and aligned"):
        current_value_survival(t, -np.ones_like(t), shape=1.3, scale=60.0)
    with pytest.raises(ValueError, match="non-negative and aligned"):
        current_value_survival(t, np.ones(3), shape=1.3, scale=60.0)
    with pytest.raises(ValueError, match="shape > 0"):
        current_value_survival(t, np.ones_like(t), shape=0.0, scale=60.0)

    short, v = np.linspace(0.0, 4.0, 5), np.array([100.0, 90.0, 80.0, 70.0, 60.0])
    assert landmark_response(short, v, 4.0) == pytest.approx(-0.4)
    assert np.isnan(landmark_response(short, v, 8.0))  # off the grid: unobserved, not -0.4


# --- CLI ---------------------------------------------------------------------

BAD_INVOCATIONS = [
    ["simulate", "no.such.record"],
    ["simulate", "survival_link.nsclc_os_growth_rate"],
    ["simulate", CLARET, "--drug-effect", "nan"],
    ["simulate", CLARET, "--drug-effect=-1"],
    ["simulate", CLARET, "--exposure", "100", "--exposure-response", "no.such"],
    ["compare", "--average", "--n", "0"],
    ["compare", "--average", "--target", "bogus"],
    ["uncertainty", CLARET, "--n", "0"],
    ["uncertainty", "no.such"],
    ["sensitivity", CLARET, "--target", "bogus", "--n", "5"],
    ["sensitivity", "growth_laws.gompertz", "--n", "5"],
    ["identify", CLARET, "--schedule", "a,b"],
    ["identify", CLARET, "--schedule=-4,0,8"],
    ["identify", CLARET, "--sigma-prop", "0", "--sigma-add", "0"],
    ["identify", "survival_link.nsclc_os_growth_rate"],
    ["design", CLARET, "--n-samples", "1"],
    ["design", CLARET, "--n-samples", "100000"],
    ["design", CLARET, "--horizon", "0"],
    ["response", CLARET, "--n", "0"],
    ["response", CLARET, "--survival-link", "bogus", "--n", "5"],
    ["pfs", CLARET, "--n", "0"],
    ["pfs", CLARET, "--landmark=-1", "--n", "5"],
    ["joint", "--alpha", "nan"],
    ["budget", "--n", "0"],
    ["interactions", CLARET, "--effect-a", "nan"],
    ["interactions", "no.such"],
    ["loewe", CLARET, "--dose-a=-1"],
    ["loewe", CLARET, "--er-a", "bogus"],
    ["dose-response", CLARET, "--c-ref", "0"],
    ["dose-response", CLARET, "--e-ref", "0"],
    ["early-surrogate", "--reference-link", "bogus"],
    ["discriminability", "--power", "0"],
    ["discriminability", "--alpha", "2"],
    ["discriminability", "--survival-link", "bogus"],
]

# Context-wide analyses over a context with no curated model.
NO_MODEL_INVOCATIONS = [
    ["atlas"], ["joint"], ["discriminability"], ["early-surrogate"], ["budget"],
    ["pfs", "--routes"], ["response", "--surrogate"], ["response", "--durability"],
]


def _run(argv, capsys):
    """Exit status + stderr, whether the error surfaced in argparse or in the command."""
    try:
        code = main(argv)
    except SystemExit as e:
        code = e.code
    return code, capsys.readouterr().err


@pytest.mark.parametrize("argv", BAD_INVOCATIONS, ids=" ".join)
def test_cli_bad_input_is_a_clean_usage_error(argv, capsys):
    code, err = _run(argv, capsys)
    assert code == 2
    assert "error:" in err and "Traceback" not in err


@pytest.mark.parametrize("argv", NO_MODEL_INVOCATIONS, ids=" ".join)
def test_cli_context_without_models_fails_instead_of_concluding(argv, capsys):
    code, err = _run([*argv, "--tumor-type", "nsclc"], capsys)  # wrong case: no such context
    assert code == 2
    assert "no eligible" in err and "NSCLC/first" in err  # names the curated contexts


def test_cli_unknown_record_names_the_near_miss(capsys):
    code, err = _run(["simulate", "resistance.claret_2009.tgj"], capsys)
    assert code == 2 and CLARET in err


def test_cli_unwritable_export_path_is_a_clean_error(tmp_path, capsys):
    blocker = tmp_path / "file"
    blocker.write_text("")
    code, err = _run(["export", "--format", "csv", "--output", str(blocker / "x")], capsys)
    assert code == 2 and "error:" in err
