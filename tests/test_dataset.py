"""Dataset integrity: schema validity, referential integrity, invariants."""

import onkos
from onkos.validate import validate_dataset


def test_dataset_validates():
    assert validate_dataset() == []


def test_load_nonempty():
    ds = onkos.load()
    assert len(ds) >= 7
    assert "resistance.claret_2009.tgi" in ds


def test_record_tier_is_worst_parameter_tier():
    ds = onkos.load()
    for r in ds:
        if not r.parameters:
            continue
        worst = max(p.tier for p in r.parameters)  # 'A' < 'B' < 'C' < 'D' lexicographically
        assert r.tier >= worst, f"{r.id}: record tier {r.tier} better than worst param {worst}"


def test_every_record_has_resolvable_primary_citation():
    ds = onkos.load()
    for r in ds:
        assert r.primary_citation is not None, r.id
        assert r.primary_citation.doi, r.id


def test_high_uncertainty_terms_carry_iiv():
    ds = onkos.load()
    claret = ds["resistance.claret_2009.tgi"]
    assert claret["lambda"].iiv_cv_percent is not None
    assert claret["lambda"].iiv_cv_percent > 50  # resistance term is poorly identified


def test_parameter_access_by_symbol():
    ds = onkos.load()
    claret = ds["resistance.claret_2009.tgi"]
    assert claret["kL"].central > 0
    assert "lambda" in claret
