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
