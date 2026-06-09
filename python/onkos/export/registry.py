"""Bind dataset records to reference kernels and expose tier propagation.

This is the seam between the curated dataset and the computational layer: given
a :class:`~onkos.models.Record`, return its :class:`KernelSpec` and the
kernel-internal parameter values. Tier/transportability propagation is
re-exported from :mod:`onkos.tiers` so format builders import one place.
"""

from __future__ import annotations

from ..models import Record
from ..tiers import propagate, transport_check, worst_tier  # noqa: F401  (re-export)
from .reference import KERNELS, KernelSpec


def get_kernel(record: Record) -> KernelSpec:
    if not record.kernel:
        raise ValueError(f"record '{record.id}' has no bound kernel")
    if record.kernel not in KERNELS:
        raise ValueError(f"record '{record.id}' references unknown kernel '{record.kernel}'")
    return KERNELS[record.kernel]


def kernel_values(record: Record) -> dict[str, float]:
    """Kernel-internal parameter values for ``record`` (record symbols mapped)."""
    spec = get_kernel(record)
    return spec.map_values(record.param_values())
