"""COMBINE .omex archive — bundles SBML + PharmML + virtual-trial JSON +
provenance + the universal prohibition manifest into one citable artifact."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

from .._const import CLINICAL_USE
from ..models import Record
from .jsonld import to_jsonld
from .pharmml import to_pharmml
from .pharmml_so import to_pharmml_so
from .registry import get_kernel
from .sbml import to_sbml
from .virtual_trial_json import to_virtual_trial_json

_MANIFEST = """<?xml version="1.0" encoding="UTF-8"?>
<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">
  <content location="." format="http://identifiers.org/combine.specifications/omex"/>
{entries}
</omexManifest>
"""

_FORMATS = {
    ".xml": "http://identifiers.org/combine.specifications/sbml",
    ".pharmml": "http://purl.org/NET/mediatypes/application/pharmml+xml",
    ".so.xml": "http://purl.org/NET/mediatypes/application/pharmml-so+xml",
    ".json": "http://purl.org/NET/mediatypes/application/json",
    ".jsonld": "http://purl.org/NET/mediatypes/application/ld+json",
}


def _suffix(name: str) -> str:
    return ".so.xml" if name.endswith(".so.xml") else Path(name).suffix


def build_omex(record: Record, out_path: str, *, tier: str | None = None) -> Path:
    spec = get_kernel(record)
    pid = record.id.replace(".", "_")
    files = {f"{pid}.json": to_virtual_trial_json(record, tier=tier)}
    files[f"{pid}.jsonld"] = to_jsonld(record, tier=tier)  # linked-data provenance
    files[f"{pid}.pharmml"] = to_pharmml(record, tier=tier)
    files[f"{pid}.so.xml"] = to_pharmml_so(record, tier=tier)  # PharmML Standard Output
    if spec.kind == "ode":
        files[f"{pid}.xml"] = to_sbml(record, tier=tier)

    files["provenance.json"] = json.dumps(
        {
            "onkos:clinicalUse": CLINICAL_USE,
            "record": record.id,
            "primary_citation": record.primary_citation.key if record.primary_citation else None,
            "tier": tier or record.tier,
        },
        indent=2,
        ensure_ascii=False,
    )

    entries = "\n".join(
        f'  <content location="./{name}" format="{_FORMATS.get(_suffix(name), "")}"/>'
        for name in files
    )
    files["manifest.xml"] = _MANIFEST.format(entries=entries)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, content in files.items():
            zf.writestr(name, content)
    return out
