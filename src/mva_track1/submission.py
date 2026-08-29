"""Build the Track 1 submission CSV.

The ``epcr`` column is the estimated probability of a causal relationship, on
``(0, 1]``, and only the *relative* ranking of our own rows is scored. Two
paths are supported and the distinction matters for honesty:

* :func:`heuristic_epcr` is a transparent, fully automated score. It is what the
  pipeline emits on a fresh sample where nobody has looked at the output yet.
* ``submitted_epcr`` in the config records the values we actually submitted,
  which were set by expert review. Running with ``--use-submitted-epcr``
  reproduces the submitted file byte-for-byte; running without it shows what the
  unattended pipeline would have said.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pandas as pd

COLUMNS = [
    "proband_id",
    "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2",
    "epcr", "finding_type", "notes",
]


def heuristic_epcr(row: dict[str, Any]) -> float:
    """Automated EPCR in (0, 1]. Deliberately simple and auditable."""
    score = 0.05
    if row.get("clinvar_sig") and "athogenic" in str(row["clinvar_sig"]):
        score += 0.30 + 0.05 * (row.get("clinvar_stars") or 0)
    if str(row.get("impact")) == "HIGH":
        score += 0.25
    elif str(row.get("impact")) == "MODERATE":
        score += 0.10
    af = row.get("af_joint")
    if af is None or pd.isna(af) or af < 1e-5:
        score += 0.15
    elif af < 1e-3:
        score += 0.08
    if row.get("phenotype_match"):
        score += 0.15
    if row.get("phase") == "in_trans":
        score += 0.10
    return round(min(score, 0.99), 3)


def write_submission(rows: list[dict[str, Any]], path: str | Path) -> Path:
    """Write rows sorted by descending EPCR. Blank second-allele columns for
    single-variant rows, as the template specifies."""
    path = Path(path)
    rows = sorted(rows, key=lambda r: -float(r["epcr"]))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})
    return path


def rows_from_config(submitted: list[dict[str, Any]], proband_id: str) -> list[dict[str, Any]]:
    out = []
    for r in submitted:
        row = {"proband_id": proband_id}
        row.update(r)
        out.append(row)
    return out
