"""PharmML export (structural model + parameters) — the durable interop anchor.

A compact PharmML 0.9-flavored document: the standardized pharmacometric markup
plus the Onkos provenance/prohibition annotations. Generated, never hand-edited.
"""

from __future__ import annotations

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def to_pharmml(record: Record, *, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    vals = kernel_values(record)
    if "E" in " ".join(spec.rhs_infix.values()):
        vals["E"] = float(drug_effect)

    pop_params = "\n".join(
        f'      <PopulationParameter symbId="{k}">\n'
        f"        <ct:Real>{v}</ct:Real>\n"
        f"      </PopulationParameter>"
        for k, v in vals.items()
    )
    state = spec.states[0]
    deriv = (
        f'      <ct:DerivativeVariable symbId="{state}"/>'
        if spec.kind == "ode"
        else f'      <ct:Variable symbId="{state}"/>'
    )
    ann = "\n".join(f"    <!-- {ln} -->" for ln in annotations_block(record, tier=tier).splitlines())

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!-- Onkos PharmML — GENERATED, do not hand-edit. -->
<PharmML xmlns="http://www.pharmml.org/pharmml/0.9/PharmML"
         xmlns:ct="http://www.pharmml.org/pharmml/0.9/CommonTypes"
         xmlns:onkos="https://onkos.dev/ns#">
{ann}
  <ModelDefinition>
    <StructuralModel blkId="sm_{_pid(record.id)}">
      <ct:name>{record.name}</ct:name>
{deriv}
    </StructuralModel>
    <ParameterModel blkId="pm_{_pid(record.id)}">
{pop_params}
    </ParameterModel>
  </ModelDefinition>
</PharmML>
"""


def _pid(s: str) -> str:
    return s.replace(".", "_")
