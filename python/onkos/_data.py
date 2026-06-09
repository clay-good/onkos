"""Locate the dataset directory (source of truth).

Resolution order:
1. ``$ONKOS_DATASET`` environment variable, if set.
2. A ``_dataset`` directory bundled inside the installed package
   (produced by ``scripts/sync_dataset_into_package.py``).
3. A ``dataset`` directory found by walking up from this file
   (the source checkout / editable install).
"""

from __future__ import annotations

import os
from pathlib import Path


def _has_records(path: Path) -> bool:
    return (path / "records").is_dir() and (path / "schema").is_dir()


def dataset_dir() -> Path:
    env = os.environ.get("ONKOS_DATASET")
    if env:
        p = Path(env).expanduser()
        if _has_records(p):
            return p
        raise FileNotFoundError(f"ONKOS_DATASET={env} does not contain records/ and schema/")

    bundled = Path(__file__).resolve().parent / "_dataset"
    if _has_records(bundled):
        return bundled

    for parent in Path(__file__).resolve().parents:
        candidate = parent / "dataset"
        if _has_records(candidate):
            return candidate

    raise FileNotFoundError(
        "Could not locate the Onkos dataset. Set ONKOS_DATASET or run from a source checkout."
    )
