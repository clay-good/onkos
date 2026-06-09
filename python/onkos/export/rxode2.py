"""nlmixr2 / rxode2 (R) export — open-source simulation/estimation target."""

from __future__ import annotations

import re

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def to_rxode2(record: Record, *, y0: float = 100.0, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError("rxode2 export supports ODE kernels only")

    vals = kernel_values(record)
    all_infix = " ".join(spec.rhs_infix.values())
    if "E" in all_infix:
        vals["E"] = float(drug_effect)
    if "y0" in all_infix:
        vals["y0"] = float(y0)

    params = ", ".join(f"{k} = {v}" for k, v in vals.items())
    odes = "\n  ".join(f"d/dt({s}) = {spec.rhs_infix[s]}" for s in spec.states)
    inits = ", ".join(f"{s} = {y0 if i == 0 else 0.0}" for i, s in enumerate(spec.states))
    obs = f"\n  tumor_size = {spec.observable}" if spec.observable else ""
    ann = "\n".join(f"# {ln}" for ln in annotations_block(record, tier=tier).splitlines())

    return f"""# Onkos rxode2/nlmixr2 model — GENERATED, do not hand-edit.
# {record.id}: {record.name}
{ann}
library(rxode2)

onkos_params <- c({params})

onkos_model <- rxode2({{
  {odes}{obs}
}})

onkos_inits <- c({inits})
# ev <- et(seq(0, 104, by = 0.5))            # weeks
# sim <- rxSolve(onkos_model, onkos_params, ev, inits = onkos_inits)
"""


def parse_rxode2_params(text: str) -> dict:
    """Re-read the parameter vector ``c(name = value, ...)`` (for round-trip)."""
    m = re.search(r"onkos_params\s*<-\s*c\(([^)]*)\)", text)
    if not m:
        return {}
    out = {}
    for pair in m.group(1).split(","):
        if "=" in pair:
            k, v = pair.split("=")
            out[k.strip()] = float(v)
    return out
