#!/usr/bin/env python3
"""Copy the source-of-truth ``dataset/`` into the package as ``_dataset/`` so a
built wheel ships the data. Run before ``python -m build``.

The dataset at the repo root remains the single source of truth; this produces a
disposable, git-ignored copy inside ``python/onkos/`` for packaging only.
"""

from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "dataset"
DST = ROOT / "python" / "onkos" / "_dataset"


def main() -> int:
    if not (SRC / "records").is_dir():
        raise SystemExit(f"no dataset at {SRC}")
    if DST.exists():
        shutil.rmtree(DST)
    for sub in ("schema", "records", "citations"):
        s = SRC / sub
        if s.is_dir():
            shutil.copytree(s, DST / sub)
    n = len(list((DST / "records").glob("*.json")))
    print(f"Synced {n} records into {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
