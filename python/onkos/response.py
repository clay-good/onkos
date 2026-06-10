"""RECIST best response & objective response rate (ORR) — the phase-2 endpoint, and
the ORR -> OS surrogate made computable.

ORR (the fraction of a trial achieving an objective RECIST response) is the dominant
go/no-go endpoint in early oncology, yet its validity as a surrogate for overall
survival is famously contested: drugs with a high response rate routinely fail to
extend survival. Onkos has OS and PFS endpoints; this adds the response endpoint and,
crucially, lets a model's ORR and OS be read off the *same* simulated trial so their
(dis)agreement is a measured quantity rather than an assumption — the response-endpoint
analog of the v0.25 survival-metric work.

Best overall response is classified per RECIST 1.1 on the tumor-size (sum-of-longest-
diameters) trajectory, from the patient's *observed baseline* ``v[0]`` (as RECIST
measures it):

* **CR** — (near-)complete disappearance (>= 95% shrinkage; an SLD-continuous proxy);
* **PR** — >= 30% shrinkage from baseline at nadir;
* **PD** — no PR and >= 20% regrowth from nadir (progression);
* **SD** — neither (stable disease).

ORR = P(CR or PR) and DCR (disease-control rate) = P(CR, PR, or SD) are *population /
trial-level* rates, computed over the stored inter-individual variability (the same
ensemble :mod:`onkos.uncertainty` uses) — NOT an individual response probability. The
module never recommends or ranks therapies; the surrogate discordance it reports is a
statement about *models*, not a treatment choice.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from ._const import CLINICAL_USE
from .load import Dataset
from .metrics import _PD_GROWTH, _PR_SHRINK
from .uncertainty import ensemble_samples

__all__ = [
    "RECIST_CATEGORIES",
    "best_response",
    "ResponseRates",
    "objective_response_rate",
    "ResponseSurvival",
    "response_vs_survival",
]

RECIST_CATEGORIES = ("CR", "PR", "SD", "PD")
_CR_SHRINK = 0.95  # (near-)complete disappearance on a continuous SLD trajectory


def best_response(t, v) -> str:
    """RECIST 1.1-style best overall response from a tumor-size trajectory ``v``,
    measured from the observed baseline ``v[0]`` (CR > PR > PD > SD)."""
    v = np.asarray(v, dtype=float)
    base = float(v[0])
    if base <= 0:
        return "SD"
    i = int(np.argmin(v))
    nadir = float(v[i])
    depth = (base - nadir) / base  # fractional shrinkage from baseline at nadir
    if depth >= _CR_SHRINK:
        return "CR"
    if depth >= _PR_SHRINK:
        return "PR"
    # No objective response: progressive disease if it regrew >= 20% above the nadir.
    post_nadir_max = float(v[i:].max())
    if nadir > 0 and post_nadir_max >= (1.0 + _PD_GROWTH) * nadir:
        return "PD"
    return "SD"


# --------------------------------------------------------------------------- #
# Population response rates over the IIV ensemble.                            #
# --------------------------------------------------------------------------- #


@dataclass
class ResponseRates:
    record_id: str
    context: dict
    n: int
    tier: str
    orr: float  # P(CR or PR) — the objective response rate
    dcr: float  # P(CR, PR, or SD) — the disease-control rate
    distribution: dict  # CR/PR/SD/PD -> fraction (sums to 1)
    median_os_weeks: float | None  # read off the SAME trial, for the surrogate question
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "record_id": self.record_id,
            "context": self.context,
            "n": self.n,
            "tier": self.tier,
            "orr": self.orr,
            "dcr": self.dcr,
            "distribution": self.distribution,
            "median_os_weeks": self.median_os_weeks,
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def objective_response_rate(
    ds: Dataset,
    record_id: str,
    *,
    context: dict | None = None,
    drug_effect: float | None = 1.0,
    exposure=None,
    exposure_response: str | None = None,
    survival_link: str | None = None,
    t: np.ndarray | None = None,
    n: int = 300,
    seed: int = 0,
) -> ResponseRates:
    """Population RECIST best-response distribution, ORR, and DCR for ``record_id`` over
    its stored IIV (the :mod:`onkos.uncertainty` ensemble), plus the median OS read off
    the same trial so the ORR -> OS surrogate can be examined. Trial-level only."""
    if t is None:
        t = np.linspace(0.0, 156.0, 313)
    s = ensemble_samples(
        ds, record_id, context=context, drug_effect=drug_effect, exposure=exposure,
        exposure_response=exposure_response, survival_link=survival_link, t=t, n=n, seed=seed,
    )
    cats = [best_response(s.t, s.tumor[i]) for i in range(s.n)]
    counts = Counter(cats)
    dist = {c: counts.get(c, 0) / s.n for c in RECIST_CATEGORIES}
    orr = dist["CR"] + dist["PR"]
    dcr = orr + dist["SD"]

    os_samples = s.median.get("OS")
    median_os = None
    if os_samples is not None:
        finite = os_samples[np.isfinite(os_samples)]
        median_os = float(np.median(finite)) if finite.size else None

    return ResponseRates(
        record_id=record_id,
        context=context or {},
        n=s.n,
        tier=s.tier,
        orr=orr,
        dcr=dcr,
        distribution=dist,
        median_os_weeks=median_os,
        warnings=list(s.warnings),
    )


# --------------------------------------------------------------------------- #
# The ORR -> OS surrogate: is the response rate a faithful ranking of survival? #
# --------------------------------------------------------------------------- #


@dataclass
class ResponseSurvival:
    context: dict
    endpoint: str
    rows: list  # per-model {record_id, tier, orr, dcr, median_os_weeks}
    discordant_pairs: int  # model pairs where ORR order and OS order disagree
    total_pairs: int
    warnings: list = field(default_factory=list)

    @property
    def discordant_fraction(self) -> float:
        """Fraction of model pairs whose ORR ranking contradicts their OS ranking — a
        nonzero value is direct evidence ORR does not faithfully track OS here."""
        return self.discordant_pairs / self.total_pairs if self.total_pairs else 0.0

    @property
    def orr_predicts_os(self) -> bool:
        return self.discordant_pairs == 0

    def to_dict(self) -> dict:
        return {
            "onkos:clinicalUse": CLINICAL_USE,
            "NOT_FOR_CLINICAL_USE": True,
            "context": self.context,
            "endpoint": self.endpoint,
            "rows": self.rows,
            "discordant_pairs": self.discordant_pairs,
            "total_pairs": self.total_pairs,
            "discordant_fraction": self.discordant_fraction,
            "orr_predicts_os": self.orr_predicts_os,
            "warnings": list(self.warnings),
        }

    def to_json(self, *, indent: int = 2) -> str:
        import json

        return json.dumps(self.to_dict(), indent=indent)


def response_vs_survival(
    ds: Dataset,
    *,
    context: dict,
    drug_effect: float = 1.0,
    survival_link: str | None = None,
    t: np.ndarray | None = None,
    n: int = 300,
    seed: int = 0,
) -> ResponseSurvival:
    """For every in-context TGI model, compute ORR and median OS from the same trial and
    count the model pairs whose ORR ranking *contradicts* their OS ranking — making the
    contested ORR -> OS surrogate a measured quantity. A discordant pair is a model with
    a higher response rate yet shorter survival than another (the failure mode that sends
    high-ORR drugs into negative phase-3 survival trials).

    Statement about *models under this context*, never a therapy ranking."""
    from .compare import compare

    if t is None:
        t = np.linspace(0.0, 156.0, 313)
    cmp = compare(ds, purpose="tgi", context=context, drug_effect=drug_effect, t=t)
    rows = []
    for tr in cmp.included:
        rr = objective_response_rate(
            ds, tr.record_id, context=context, drug_effect=drug_effect,
            survival_link=survival_link, t=t, n=n, seed=seed,
        )
        rows.append({
            "record_id": tr.record_id,
            "tier": rr.tier,
            "orr": rr.orr,
            "dcr": rr.dcr,
            "median_os_weeks": rr.median_os_weeks,
        })

    # Count discordant pairs among models with both an ORR and a finite median OS.
    scored = [r for r in rows if r["median_os_weeks"] is not None]
    discordant = 0
    total = 0
    for i in range(len(scored)):
        for j in range(i + 1, len(scored)):
            a, b = scored[i], scored[j]
            d_orr = a["orr"] - b["orr"]
            d_os = a["median_os_weeks"] - b["median_os_weeks"]
            if d_orr == 0 or d_os == 0:
                continue  # a tie is neither concordant nor discordant
            total += 1
            if (d_orr > 0) != (d_os > 0):  # higher ORR but shorter OS (or vice versa)
                discordant += 1

    return ResponseSurvival(
        context=context,
        endpoint="OS",
        rows=rows,
        discordant_pairs=discordant,
        total_pairs=total,
    )
