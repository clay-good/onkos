"""MIRIAM-style annotations carried by every export.

Each export carries: a dataset-version pin; the propagated confidence tier and
any transportability/failure-mode warnings; ``bqbiol:isDescribedBy`` DOI/PMID
links that survive even if a tool strips the custom ``onkos:`` predicates; and a
universal, machine-readable ``onkos:clinicalUse`` prohibition.
"""

from __future__ import annotations

from .._const import CLINICAL_USE
from ..models import Record


def clinical_use_rdf() -> str:
    return f'onkos:clinicalUse "{CLINICAL_USE}"'


def identifier_uris(record: Record) -> list[str]:
    uris: list[str] = []
    cit = record.primary_citation
    if cit and cit.doi:
        uris.append(f"https://identifiers.org/doi/{cit.doi}")
    if cit and cit.pmid:
        uris.append(f"https://identifiers.org/pubmed:{cit.pmid}")
    return uris


def annotations_block(
    record: Record, *, tier: str | None = None, dataset_version: str = "0.1.0", warnings=None
) -> str:
    """A compact, format-agnostic annotation block (RDF-flavored text).

    Format builders embed this verbatim (as comments or as an RDF island) so the
    same provenance and prohibition travel with every artifact.
    """
    tier = tier or record.tier
    lines = [
        f"onkos:datasetVersion {dataset_version}",
        f"onkos:confidenceTier {tier}",
        clinical_use_rdf(),
    ]
    for uri in identifier_uris(record):
        lines.append(f"bqbiol:isDescribedBy <{uri}>")
    dc = record.derivation_context
    if dc:
        lines.append(
            "onkos:derivationContext "
            f"drug={dc.drug}; class={dc.drug_class}; tumor={dc.tumor_type}; line={dc.line_of_therapy}"
        )
    tp = record.transportability
    if tp:
        lines.append(
            "onkos:transportability "
            f"tumor_types={list(tp.validated_tumor_types)}; lines={list(tp.validated_lines)}"
        )
    for w in warnings or []:
        lines.append(f"onkos:warning {w}")
    return "\n".join(lines)


def sbml_rdf_xml(record: Record, *, tier: str | None = None, dataset_version: str = "0.1.0") -> str:
    """RDF/XML island embeddable in an SBML <annotation> element."""
    tier = tier or record.tier
    desc = []
    for uri in identifier_uris(record):
        desc.append(f'      <rdf:li rdf:resource="{uri}"/>')
    bag = "\n".join(desc)
    return (
        '    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n'
        '             xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"\n'
        '             xmlns:onkos="https://onkos.dev/ns#">\n'
        f'      <rdf:Description rdf:about="#{record.id}">\n'
        f'        <onkos:clinicalUse>{CLINICAL_USE}</onkos:clinicalUse>\n'
        f"        <onkos:confidenceTier>{tier}</onkos:confidenceTier>\n"
        f"        <onkos:datasetVersion>{dataset_version}</onkos:datasetVersion>\n"
        "        <bqbiol:isDescribedBy>\n"
        "          <rdf:Bag>\n"
        f"{bag}\n"
        "          </rdf:Bag>\n"
        "        </bqbiol:isDescribedBy>\n"
        "      </rdf:Description>\n"
        "    </rdf:RDF>"
    )
