"""PharmML Standard Output (SO) export — the results companion to PharmML.

The COMBINE convention pairs a PharmML *model* with a PharmML *Standard Output*
(SO) carrying the estimation results (spec §7: ".omex bundles SBML + PharmML +
SO + provenance"). For Onkos's curated published parameters the SO carries:

- ``PopulationEstimates`` (MLE) — the central parameter values;
- the inter-individual variability as **random-effect variances**, computed from
  the reported CV as the lognormal variance ``omega = ln(1 + CV^2)``. This is
  honest: IIV is between-subject variability, *not* the precision (RSE) of the
  population estimate — which the dataset does not curate, and which is therefore
  deliberately omitted rather than faked;
- the external-validation performance as ``ModelDiagnostic`` entries;
- the universal Onkos annotations (clinicalUse, tier, DOI, predictionStatus).
"""

from __future__ import annotations

import math
import re

from ..models import Record
from .annotate import annotations_block
from .registry import get_kernel, kernel_values


def _iiv_variance(record: Record) -> dict:
    """kernel-internal name -> lognormal random-effect variance (omega)."""
    spec = get_kernel(record)
    sym_to_kernel = dict(zip(spec.record_symbols, spec.params))
    out = {}
    for p in record.parameters:
        if p.iiv_cv_percent and p.symbol in sym_to_kernel:
            cv = p.iiv_cv_percent / 100.0
            out[sym_to_kernel[p.symbol]] = math.log(1.0 + cv * cv)
    return out


def to_pharmml_so(record: Record, *, drug_effect: float = 1.0, tier=None) -> str:
    spec = get_kernel(record)
    vals = kernel_values(record)
    if "E" in " ".join(spec.rhs_infix.values()):
        vals["E"] = float(drug_effect)
    omegas = _iiv_variance(record)

    mle_rows = "\n".join(
        "          <ct:Row>"
        f"<ct:String>{k}</ct:String><ct:Real>{v}</ct:Real>"
        "</ct:Row>"
        for k, v in vals.items()
    )
    iiv_rows = "\n".join(
        f'      <onkos:randomEffectVariance parameter="{k}" omega="{v}"/>'
        for k, v in omegas.items()
    )
    diag_rows = "\n".join(
        f'      <onkos:externalValidation metric="{pp.metric}" value="{pp.value}" '
        f'population="{pp.population or ""}"/>'
        for pp in record.predictive_performance
    )
    ann = "\n".join(f"  <!-- {ln} -->" for ln in annotations_block(record, tier=tier).splitlines())

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!-- Onkos PharmML Standard Output (SO) — GENERATED, do not hand-edit. -->
<SO xmlns="http://www.pharmml.org/so/0.3/StandardOutput"
    xmlns:ct="http://www.pharmml.org/pharmml/0.9/CommonTypes"
    xmlns:onkos="https://onkos.dev/ns#">
{ann}
  <SOBlock blkId="SO_{_pid(record.id)}">
    <Estimation>
      <PopulationEstimates>
        <MLE>
          <ct:DataSet>
            <ct:Definition>
              <ct:Column columnId="parameter" valueType="string" columnNum="1"/>
              <ct:Column columnId="estimate" valueType="real" columnNum="2"/>
            </ct:Definition>
            <ct:Table>
{mle_rows}
            </ct:Table>
          </ct:DataSet>
        </MLE>
      </PopulationEstimates>
      <onkos:interIndividualVariability note="lognormal random-effect variance omega = ln(1+CV^2); IIV is not estimate precision (RSE)">
{iiv_rows or "      <!-- no IIV reported -->"}
      </onkos:interIndividualVariability>
    </Estimation>
    <ModelDiagnostic>
{diag_rows or "      <!-- no external validation recorded -->"}
    </ModelDiagnostic>
  </SOBlock>
</SO>
"""


def _pid(s: str) -> str:
    return s.replace(".", "_")


def parse_so_estimates(text: str) -> dict:
    """Re-read the MLE parameter estimates from an SO document (for round-trip)."""
    out = {}
    for m in re.finditer(
        r"<ct:String>([^<]+)</ct:String><ct:Real>([-\d.eE+]+)</ct:Real>", text
    ):
        out[m.group(1)] = float(m.group(2))
    return out
