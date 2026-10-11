"""Every ``--json`` output is strictly valid JSON.

Non-finite results are routine here — a regrowth rate that never occurred (``nan``), an
inestimable parameter (``inf`` RSE) — and Python would serialize them as the non-JSON
tokens ``NaN`` / ``Infinity`` that ``jq`` and ``JSON.parse`` reject. They must be ``null``.
"""

import json

import pytest
from onkos._serialize import to_strict_json
from onkos.cli import main

CLARET = "resistance.claret_2009.tgi"
TWO_POP = "resistance.nsclc_first_line.two_population"

JSON_INVOCATIONS = [
    ["simulate", "--compare", "--json"],  # Norton-Simon eradicates: k_g is nan
    ["simulate", "--compare", "--json", "--include-curves"],
    ["compare", "--json"],
    ["compare", "--average", "--json", "--n", "20"],
    ["identify", CLARET, "--json"],
    ["identify", TWO_POP, "--json"],  # rank-deficient design: inf RSE / collinearity
    ["identify", CLARET, "--json", "--drug-effect", "0"],
    ["design", TWO_POP, "--json"],
    ["response", CLARET, "--json", "--n", "20"],
    ["response", "--surrogate", "--json", "--n", "20"],
    ["pfs", CLARET, "--json", "--n", "20", "--drug-effect", "5"],  # never progresses
    ["pfs", "--routes", "--json", "--n", "20"],
    ["budget", "--json", "--n", "10"],
    ["interactions", CLARET, "--json"],
    ["loewe", CLARET, "--json"],
    ["dose-response", CLARET, "--json"],
    ["early-surrogate", "--json"],
    ["discriminability", "--json"],
    ["atlas", "--json"],
    ["review-queue", "--json"],
    ["simulate", CLARET, "--json"],
    ["simulate", "drug_effect.norton_simon.nsclc", "--json", "--include-curves"],  # nan k_g
    ["uncertainty", CLARET, "--json", "--n", "20"],
    ["uncertainty", "drug_effect.norton_simon.nsclc", "--json", "--n", "20", "--include-curves"],
    ["sensitivity", CLARET, "--json", "--n", "30"],
    ["joint", "--json"],
    ["info", "--json"],
    ["audit", "--json"],
    ["validate", "--json"],
    ["list", "--json"],
]


def _reject(token):
    raise AssertionError(f"non-JSON token {token} in output")


@pytest.mark.parametrize("argv", JSON_INVOCATIONS, ids=" ".join)
def test_cli_json_is_strict(argv, capsys):
    assert main(argv) == 0
    json.loads(capsys.readouterr().out, parse_constant=_reject)


def test_non_finite_values_serialize_as_null():
    payload = {"a": float("nan"), "b": [1.5, float("inf"), (-float("inf"),)], "c": "nan", "d": 2}
    assert json.loads(to_strict_json(payload)) == {
        "a": None, "b": [1.5, None, [None]], "c": "nan", "d": 2,
    }


def test_every_command_with_a_result_offers_json():
    """Machine-readable output is the rule, not a per-command courtesy: every subcommand
    takes ``--json`` except the three whose output is already a file format or a bare
    string (``report`` writes Markdown, ``export`` writes artifacts, ``version``)."""
    from onkos.cli import build_parser

    subparsers = next(
        a for a in build_parser()._actions if a.__class__.__name__ == "_SubParsersAction"
    )
    without = {
        name for name, sp in subparsers.choices.items()
        if "--json" not in sp._option_string_actions
    }
    assert without == {"version", "report", "export"}


def test_result_objects_round_trip_through_json():
    import onkos

    ds = onkos.load()
    ctx = {"tumor_type": "NSCLC", "line": "first"}
    tr = onkos.simulate(ds, CLARET, context=ctx)
    d = json.loads(tr.to_json(include_curves=True))
    assert d["NOT_FOR_CLINICAL_USE"] is True and d["tier"] == tr.tier
    assert d["median_os_weeks"] == tr.median_os and len(d["tumor_size"]) == tr.t.size
    assert set(d["survival"]) == {"OS", "PFS"}

    ens = onkos.simulate_ensemble(ds, CLARET, context=ctx, n=10)
    e = json.loads(ens.to_json())
    assert e["n"] == 10 and e["metrics"]["median_os_weeks"] == ens.median_os
    assert "tumor_size" not in e  # curves only on request

    res = onkos.sensitivity(ds, CLARET, context=ctx, n=40)
    s = json.loads(res.to_json())
    assert s["dominant"] == res.dominant.symbol
    assert sum(p["contribution"] for p in s["parameters"]) == pytest.approx(1.0)
