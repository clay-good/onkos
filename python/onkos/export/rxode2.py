"""nlmixr2 / rxode2 (R) export — open-source simulation/estimation target."""

from __future__ import annotations

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def _to_r_ode(infix: str, state: str) -> str:
    expr = infix.replace(state, state)  # rxode2 uses the state name directly
    return f"d/dt({state}) = {expr}"


def to_rxode2(record: Record, *, y0: float = 100.0, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError("rxode2 export supports ODE kernels only")

    vals = kernel_values(record)
    infix = spec.rhs_infix[spec.states[0]]
    if "E" in infix:
        vals["E"] = float(drug_effect)
    if "y0" in infix:
        vals["y0"] = float(y0)

    params = ", ".join(f"{k} = {v}" for k, v in vals.items())
    ode = _to_r_ode(infix, spec.states[0])
    ann = "\n".join(f"# {ln}" for ln in annotations_block(record, tier=tier).splitlines())
    state = spec.states[0]

    return f"""# Onkos rxode2/nlmixr2 model — GENERATED, do not hand-edit.
# {record.id}: {record.name}
{ann}
library(rxode2)

onkos_params <- c({params})

onkos_model <- rxode2({{
  {ode}
}})

onkos_inits <- c({state} = {y0})
# ev <- et(seq(0, 104, by = 0.5))            # weeks
# sim <- rxSolve(onkos_model, onkos_params, ev, inits = onkos_inits)
"""
