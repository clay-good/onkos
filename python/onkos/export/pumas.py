"""Pumas (Julia) export — open-source simulation/estimation target."""

from __future__ import annotations

import re

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def to_pumas(record: Record, *, y0: float = 100.0, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError("Pumas export supports ODE kernels only")

    vals = kernel_values(record)
    all_infix = " ".join(spec.rhs_infix.values())
    if "E" in all_infix:
        vals["E"] = float(drug_effect)
    if "y0" in all_infix:
        vals["y0"] = float(y0)

    # Julia uses exp/log natively; rename infix function `ln` -> `log`.
    tv = "\n".join(f"        tv{k} = {v}" for k, v in vals.items())
    pre = "\n".join(f"        {k} = tv{k}" for k in vals)
    init = "\n".join(
        f"        {s} = {y0 if i == 0 else 0.0}" for i, s in enumerate(spec.states)
    )
    dyn = "\n".join(
        f"        {s}' = {spec.rhs_infix[s].replace('ln(', 'log(')}" for s in spec.states
    )
    ann = "\n".join(f"# {ln}" for ln in annotations_block(record, tier=tier).splitlines())

    return f"""# Onkos Pumas model — GENERATED, do not hand-edit.
# {record.id}: {record.name}
{ann}
using Pumas

onkos_model = @model begin
    @param begin
{tv}
    end
    @pre begin
{pre}
    end
    @init begin
{init}
    end
    @dynamics begin
{dyn}
    end
end
"""


def parse_pumas_params(text: str) -> dict:
    """Re-read the @param ``tv<name> = value`` typical values (for round-trip)."""
    return {m.group(1): float(m.group(2)) for m in re.finditer(r"tv(\w+)\s*=\s*([-\d.eE+]+)", text)}
