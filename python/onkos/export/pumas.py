"""Pumas (Julia) export — open-source simulation/estimation target."""

from __future__ import annotations

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def to_pumas(record: Record, *, y0: float = 100.0, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError("Pumas export supports ODE kernels only")

    vals = kernel_values(record)
    infix = spec.rhs_infix[spec.states[0]]
    if "E" in infix:
        vals["E"] = float(drug_effect)
    if "y0" in infix:
        vals["y0"] = float(y0)

    state = spec.states[0]
    # Julia uses exp/log natively; rename infix function `ln` -> `log`.
    ode = infix.replace("ln(", "log(")
    tv = "\n".join(f"        tv{k} = {v}" for k, v in vals.items())
    pre = "\n".join(f"        {k} = tv{k}" for k in vals)
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
        {state} = {y0}
    end
    @dynamics begin
        {state}' = {ode}
    end
end
"""
