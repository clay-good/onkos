"""SBML L3v2 export with bidirectional infix<->MathML conversion.

The MathML round-trip (``tests/test_roundtrip.py``) re-parses the generated
rate law back to an evaluable expression and checks it numerically against the
reference kernel's hand-written ``rhs`` — so a serialization bug cannot ship.
"""

from __future__ import annotations

import ast
import xml.etree.ElementTree as ET

from ..models import Record
from .annotate import sbml_rdf_xml
from .registry import get_kernel, kernel_values

MATHML_NS = "http://www.w3.org/1998/Math/MathML"
TIME_URL = "http://www.sbml.org/sbml/symbols/time"

_BINOP = {ast.Add: "plus", ast.Sub: "minus", ast.Mult: "times", ast.Div: "divide", ast.Pow: "power"}
_FUNC = {"exp": "exp", "ln": "ln", "log": "ln"}


def _node_to_mathml(node) -> str:
    if isinstance(node, ast.Expression):
        return _node_to_mathml(node.body)
    if isinstance(node, ast.BinOp):
        op = _BINOP[type(node.op)]
        return f"<apply><{op}/>{_node_to_mathml(node.left)}{_node_to_mathml(node.right)}</apply>"
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return f"<apply><minus/>{_node_to_mathml(node.operand)}</apply>"
    if isinstance(node, ast.Call):
        fn = _FUNC[node.func.id]
        return f"<apply><{fn}/>{_node_to_mathml(node.args[0])}</apply>"
    if isinstance(node, ast.Name):
        if node.id == "t":
            return f'<csymbol encoding="text" definitionURL="{TIME_URL}">t</csymbol>'
        return f"<ci>{node.id}</ci>"
    if isinstance(node, ast.Constant):
        return f"<cn>{node.value}</cn>"
    raise ValueError(f"unsupported expression node: {ast.dump(node)}")


def infix_to_mathml(expr: str) -> str:
    body = _node_to_mathml(ast.parse(expr, mode="eval"))
    return f'<math xmlns="{MATHML_NS}">{body}</math>'


def _local(tag: str) -> str:
    return tag.split("}")[-1]


def _elem_to_infix(el) -> str:
    tag = _local(el.tag)
    if tag == "math":
        return _elem_to_infix(list(el)[0])
    if tag == "cn":
        return el.text.strip()
    if tag == "ci":
        return el.text.strip()
    if tag == "csymbol":
        url = el.get("definitionURL", "")
        return "t" if "time" in url else el.text.strip()
    if tag == "apply":
        children = list(el)
        op = _local(children[0].tag)
        args = [_elem_to_infix(c) for c in children[1:]]
        if op == "plus":
            return "(" + " + ".join(args) + ")"
        if op == "minus":
            return f"(-{args[0]})" if len(args) == 1 else f"({args[0]} - {args[1]})"
        if op == "times":
            return "(" + " * ".join(args) + ")"
        if op == "divide":
            return f"({args[0]} / {args[1]})"
        if op == "power":
            return f"({args[0]} ** {args[1]})"
        if op in ("exp", "ln"):
            fn = "exp" if op == "exp" else "ln"
            return f"{fn}({args[0]})"
        raise ValueError(f"unsupported MathML operator: {op}")
    raise ValueError(f"unsupported MathML element: {tag}")


def mathml_to_infix(mathml: str) -> str:
    return _elem_to_infix(ET.fromstring(mathml))


def _sbml_parameters(record: Record, y0: float, drug_effect: float):
    spec = get_kernel(record)
    vals = kernel_values(record)
    params = dict(vals)
    infix = " ".join(spec.rhs_infix.values())
    if "E" in infix:
        params["E"] = float(drug_effect)
    if "y0" in infix:
        params["y0"] = float(y0)
    return params


def to_sbml(record: Record, *, y0: float = 100.0, drug_effect: float = 1.0, tier=None) -> str:
    """Generate an SBML L3v2 document for a (possibly multi-state) ODE record."""
    spec = get_kernel(record)
    if spec.kind != "ode":
        raise ValueError(f"SBML export supports ODE kernels only; '{record.kernel}' is {spec.kind}")

    params = _sbml_parameters(record, y0, drug_effect)
    p_xml = "\n".join(
        f'      <parameter id="{k}" value="{v}" constant="true"/>' for k, v in params.items()
    )
    # The seed input fills the first state's initial amount; other states start at 0.
    species_xml = "\n".join(
        f'      <species id="{s}" compartment="body" initialAmount="{y0 if i == 0 else 0.0}" '
        'hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>'
        for i, s in enumerate(spec.states)
    )
    rules_xml = "\n".join(
        f'      <rateRule variable="{s}">\n        {infix_to_mathml(spec.rhs_infix[s])}\n'
        "      </rateRule>"
        for s in spec.states
    )

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="{_xml_id(record.id)}" name="{_xml_attr(record.name)}">
    <annotation>
{sbml_rdf_xml(record, tier=tier)}
    </annotation>
    <listOfCompartments>
      <compartment id="body" spatialDimensions="3" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
{species_xml}
    </listOfSpecies>
    <listOfParameters>
{p_xml}
    </listOfParameters>
    <listOfRules>
{rules_xml}
    </listOfRules>
  </model>
</sbml>
"""


def _xml_id(s: str) -> str:
    return s.replace(".", "_")


def _xml_attr(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace('"', "&quot;")
