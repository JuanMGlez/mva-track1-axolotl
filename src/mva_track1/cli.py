"""Command-line entry point: ``python -m mva_track1 --vcf ... --out results/``."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config
from .pipeline import run


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mva_track1",
        description="MVA Hackathon 2026 Track 1 pipeline (team Axolotl): "
                    "HPO-anchored targeted panel with ACMG/AMP evidence weighting.",
    )
    p.add_argument("--vcf", required=True, type=Path, help="bgzipped single-sample VCF (GRCh38)")
    p.add_argument("--config", type=Path, default=Path("config/panels.yaml"))
    p.add_argument("--out", type=Path, default=Path("results"))
    p.add_argument("--cache-dir", type=Path, default=Path(".cache/annotations"),
                   help="on-disk annotation cache; pass an empty string to disable")
    p.add_argument("--heuristic-epcr", action="store_true",
                   help="emit the unattended pipeline's own EPCR scores instead of "
                        "the expert-reviewed values recorded in the config")
    p.add_argument("--no-figure", action="store_true")
    a = p.parse_args(argv)

    paths = run(
        vcf=a.vcf,
        cfg=load_config(a.config),
        outdir=a.out,
        cache_dir=(a.cache_dir if str(a.cache_dir) else None),
        use_submitted_epcr=not a.heuristic_epcr,
        make_figure=not a.no_figure,
    )
    print("\n".join(f"{k:12s} {v}" for k, v in paths.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
