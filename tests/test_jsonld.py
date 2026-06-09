"""JSON-LD linked-data export — validated by expanding to RDF triples with rdflib."""

import json

import onkos
import pytest
from onkos.export.jsonld import dataset_jsonld, load_context, record_node, to_jsonld

ONKOS_NS = "https://onkos.dev/ns#"
BQBIOL = "http://biomodels.net/biology-qualifiers/"
REC = "https://onkos.dev/record/resistance.claret_2009.tgi"


def test_context_is_loaded_from_dataset():
    ctx = load_context()
    assert ctx["onkos"] == ONKOS_NS
    assert ctx["clinicalUse"] == "onkos:clinicalUse"
    # isDescribedBy is typed as an IRI so its values become RDF resources
    assert ctx["isDescribedBy"]["@type"] == "@id"


def test_record_node_structure():
    node = record_node(onkos.load()["resistance.claret_2009.tgi"])
    assert node["@id"] == REC
    assert node["@type"] == "Model"
    assert node["confidenceTier"] == "C"
    assert node["hasParameter"] and node["hasParameter"][0]["symbol"] == "kL"


def test_document_is_valid_json():
    doc = json.loads(to_jsonld(onkos.load()["resistance.claret_2009.tgi"]))
    assert "@context" in doc and doc["@id"] == REC


def test_hypothesis_tier_carries_prediction_status():
    node = record_node(onkos.load()["immuno_oncology.kuznetsov_1994.tumor_immune"])
    assert node["predictionStatus"] == "DO NOT USE FOR PREDICTION (hypothesis-tier)"


# --- RDF-level validation (rdflib) -----------------------------------------
rdflib = pytest.importorskip("rdflib")


def _graph(data):
    g = rdflib.Graph()
    g.parse(data=data, format="json-ld")
    return g


def test_record_expands_to_expected_triples():
    g = _graph(to_jsonld(onkos.load()["resistance.claret_2009.tgi"]))
    rec = rdflib.URIRef(REC)
    onk = rdflib.Namespace(ONKOS_NS)
    bq = rdflib.Namespace(BQBIOL)

    clinical = list(g.objects(rec, onk.clinicalUse))
    assert clinical and "PROHIBITED" in str(clinical[0])
    assert str(next(g.objects(rec, onk.confidenceTier))) == "C"

    # DOI/PMID become real RDF resources (IRIs), not strings
    described = {str(o) for o in g.objects(rec, bq.isDescribedBy)}
    assert "https://identifiers.org/doi/10.1200/JCO.2008.21.0807" in described
    assert all(isinstance(o, rdflib.URIRef) for o in g.objects(rec, bq.isDescribedBy))


def test_dataset_graph_covers_every_record():
    ds = onkos.load()
    g = _graph(dataset_jsonld(ds))
    onk = rdflib.Namespace(ONKOS_NS)
    record_ids = {str(s) for s in g.subjects(onk.recordId, None)}
    assert len(record_ids) == len(ds)


def test_vt_json_is_valid_jsonld():
    from onkos.export import to_virtual_trial_json

    g = _graph(to_virtual_trial_json(onkos.load()["resistance.claret_2009.tgi"]))
    onk = rdflib.Namespace(ONKOS_NS)
    rec = rdflib.URIRef(REC)
    assert "PROHIBITED" in str(next(g.objects(rec, onk.clinicalUse)))
