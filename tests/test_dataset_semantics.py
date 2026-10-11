"""Semantic integrity: ``onkos validate`` rejects records that are well-formed JSON but
inconsistent or not computable.

The JSON Schema proves shape, not meaning. Each case below corrupts exactly one record of
an otherwise-valid copy of the dataset in a way the schema accepts, and the validator
must name the problem — instead of the record crashing at simulation time or, worse,
silently changing results (a typo'd ``link_metric`` used to quietly fall back to the
baseline hazard).
"""

import copy
import json
import shutil

import pytest
from onkos._data import dataset_dir
from onkos.validate import validate_dataset

CLARET = "resistance.claret_2009.tgi"
WEIBULL = "survival_link.nsclc_os_week8"
COX = "survival_link.nsclc_os_cox"
BASELINE = "tumor_type_baselines.nsclc_first_line"
EMAX = "exposure_response.emax_generic"


def _param(record, symbol):
    return next(p for p in record["parameters"] if p["symbol"] == symbol)


def _drop(symbol):
    def mutate(r):
        r["parameters"] = [p for p in r["parameters"] if p["symbol"] != symbol]
    return mutate


def _rename(symbol, new):
    def mutate(r):
        _param(r, symbol)["symbol"] = new
    return mutate


def _duplicate(symbol, as_symbol=None):
    def mutate(r):
        clone = copy.deepcopy(_param(r, symbol))
        clone["symbol"] = as_symbol or symbol
        r["parameters"].append(clone)
    return mutate


def _set(symbol, **fields):
    def mutate(r):
        p = _param(r, symbol)
        for key, value in fields.items():
            if key in ("central", "low", "high", "units"):
                p["value"][key] = value
            else:
                p[key] = value
    return mutate


def _structure(**fields):
    def mutate(r):
        r["structure"].update(fields)
    return mutate


def _table(key, index, value):
    def mutate(r):
        r["structure"]["baseline_survival"][key][index] = value
    return mutate


def _truncate_table(r):
    r["structure"]["baseline_survival"]["survival"].pop()


def _no_table(r):
    del r["structure"]["baseline_survival"]


def _no_kernel(r):
    del r["kernel"]


def _derived_elsewhere(r):
    r["derivation_context"]["tumor_type"] = "breast"


def _clone_as(new_id):
    def mutate(r):
        r["id"] = new_id
    return mutate


CASES = [
    (CLARET, _drop("lambda"), "requires parameter(s) ['lambda']"),
    (CLARET, _rename("lambda", "lam_typo"), "requires parameter(s) ['lambda']"),
    (CLARET, _duplicate("kL"), "'kL' is declared 2 times"),
    (CLARET, _duplicate("kL", "mystery"), "['mystery'] are not used by kernel"),
    (CLARET, _no_kernel, "no bound kernel"),
    (CLARET, _set("kL", central=-0.02), "is negative"),
    (CLARET, _set("lambda", central=0.0), "does not compute at its central values"),
    (CLARET, _set("kD", iiv_cv_percent=-40), "iiv_cv_percent must be positive"),
    (CLARET, _set("kL", low=0.5, high=0.9), "outside [low, high]"),
    (CLARET, _set("kL", units=" "), "empty units"),
    (CLARET, _derived_elsewhere, "not in its own validated_tumor_types"),
    (WEIBULL, _set("weibull_shape", central=0.0), "not a valid survival function"),
    (WEIBULL, _structure(link_metric="no_such_metric"), "is not a TGI metric"),
    (COX, _structure(default=True), "2 default OS survival links for NSCLC/first"),
    (COX, _truncate_table, "equal length"),
    (COX, _table("survival", 3, 0.99), "never increase"),
    (COX, _table("survival", 0, 1.4), "must start at 1"),
    (COX, _table("times", 3, 5), "strictly increase"),
    (COX, _no_table, "needs structure.baseline_survival"),
    (BASELINE, _drop("baseline_tumor_size"), "no 'baseline_tumor_size' parameter"),
    (BASELINE, _set("baseline_tumor_size", central=0.0, low=None), "must be positive"),
    (BASELINE, _clone_as(BASELINE + "_copy"), "2 context_baseline records for NSCLC/first"),
    (EMAX, _set("EC50", central=0.0), "does not compute at its central values"),
]


@pytest.fixture(scope="module")
def pristine(tmp_path_factory):
    path = tmp_path_factory.mktemp("pristine") / "dataset"
    shutil.copytree(dataset_dir(), path)
    return path


@pytest.mark.parametrize(
    ("record_id", "mutate", "expected"), CASES, ids=[expected for *_, expected in CASES]
)
def test_validate_rejects_semantically_broken_record(
    pristine, tmp_path, record_id, mutate, expected
):
    base = tmp_path / "dataset"
    shutil.copytree(pristine, base)
    source = base / "records" / f"{record_id}.json"
    record = json.loads(source.read_text())
    mutate(record)
    (base / "records" / f"{record['id']}.json").write_text(json.dumps(record))

    errors = validate_dataset(str(base))

    assert any(expected in error for error in errors), errors


def test_pristine_copy_is_valid(pristine):
    """The control: the mutations above, not the copying, are what fail."""
    assert validate_dataset(str(pristine)) == []
