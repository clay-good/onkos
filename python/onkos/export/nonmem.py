"""NONMEM control-stream export (ADVAN13 general-ODE form).

TGI and survival-link models are overwhelmingly NONMEM-fit; this is the format
pharmacometricians read first. The generated stream is a faithful, runnable
skeleton with the dataset's parameter values as ``$THETA`` initial estimates and
the kernel dynamics in ``$DES``.
"""

from __future__ import annotations

import re

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def _des_line(infix: str, state: str, names: list[str]) -> str:
    expr = infix.replace(state, "A(1)")
    # uppercase known identifiers (longest first to avoid partial overlaps)
    for n in sorted(names, key=len, reverse=True):
        expr = re.sub(rf"\b{re.escape(n)}\b", n.upper(), expr)
    expr = expr.replace("exp(", "EXP(").replace("ln(", "LOG(")
    return f"  DADT(1) = {expr}"


def to_nonmem(record: Record, *, y0: float = 100.0, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError(f"NONMEM export supports ODE kernels only; '{record.kernel}' is {spec.kind}")

    vals: dict[str, float] = kernel_values(record)
    infix = spec.rhs_infix[spec.states[0]]
    if "E" in infix:
        vals["E"] = float(drug_effect)
    if "y0" in infix:
        vals["y0"] = float(y0)

    names = list(vals)
    pk_lines = "\n".join(f"  {n.upper()} = THETA({i + 1})" for i, n in enumerate(names))
    theta_lines = "\n".join(
        f"  (0, {vals[n]})  ; {n.upper()}" for n in names
    )
    des = _des_line(infix, spec.states[0], names)
    ann = "\n".join(f"; {ln}" for ln in annotations_block(record, tier=tier).splitlines())

    return f""";; Onkos NONMEM control stream — GENERATED, do not hand-edit.
;; Model: {record.id} ({record.name})
{ann}
$PROBLEM {record.id}
$INPUT ID TIME DV MDV
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9
$MODEL
  COMP=(TUMOR)
$PK
{pk_lines}
  A_0(1) = {y0}
$DES
{des}
$ERROR
  IPRED = A(1)
  Y = IPRED*(1 + EPS(1))
$THETA
{theta_lines}
$OMEGA 0 FIX
$SIGMA 0.04
$ESTIMATION METHOD=1 INTER MAXEVAL=9999
"""


def parse_nonmem_thetas(text: str) -> list[float]:
    """Extract $THETA initial estimates (for round-trip verification)."""
    out: list[float] = []
    in_block = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("$THETA"):
            in_block = True
            continue
        if in_block:
            if stripped.startswith("$"):
                break
            m = re.search(r"\(0,\s*([-\d.eE+]+)\)", stripped)
            if m:
                out.append(float(m.group(1)))
    return out
