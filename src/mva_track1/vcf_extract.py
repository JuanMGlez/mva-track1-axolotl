"""Interval extraction and genotype/QC derivation from a single-sample VCF.

Uses ``pysam`` for indexed random access when it is installed, and falls back to
a single streaming pass over the bgzipped file otherwise (no external binaries
required). Contig names are normalised so a ``chr``-prefixed and an unprefixed
callset behave identically.
"""

from __future__ import annotations

import gzip
import importlib.util
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd


@dataclass
class Interval:
    chrom: str
    start: int
    end: int
    gene: str


def _norm(chrom: str) -> str:
    return chrom.replace("chr", "")


def _iter_streaming(vcf: Path, ivs: list[Interval]) -> Iterator[tuple[Interval, list[str]]]:
    by_chrom: dict[str, list[Interval]] = {}
    for iv in ivs:
        by_chrom.setdefault(_norm(iv.chrom), []).append(iv)
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            for iv in by_chrom.get(_norm(f[0]), ()):
                pos = int(f[1])
                if iv.start <= pos <= iv.end:
                    yield iv, f
                    break


def _iter_tabix(vcf: Path, ivs: list[Interval]) -> Iterator[tuple[Interval, list[str]]]:
    import pysam

    tbx = pysam.TabixFile(str(vcf))
    contigs = {_norm(c): c for c in tbx.contigs}
    for iv in ivs:
        contig = contigs.get(_norm(iv.chrom))
        if contig is None:
            continue
        for row in tbx.fetch(contig, iv.start - 1, iv.end):
            yield iv, row.split("\t")


def extract(vcf: str | Path, intervals: Iterable[Interval], prefer_tabix: bool = True) -> pd.DataFrame:
    """All VCF records overlapping ``intervals``, one row per record."""
    vcf, ivs = Path(vcf), list(intervals)
    # importlib.util.find_spec, not try/except around the call: _iter_tabix is a
    # generator, so its ImportError would not surface until first iteration.
    have_pysam = importlib.util.find_spec("pysam") is not None
    use_tabix = prefer_tabix and have_pysam and Path(str(vcf) + ".tbi").exists()
    it = _iter_tabix(vcf, ivs) if use_tabix else _iter_streaming(vcf, ivs)

    rows = []
    for iv, f in it:
        rows.append(
            dict(
                gene=iv.gene,
                chrom=_norm(f[0]),
                pos=int(f[1]),
                vcf_id=f[2],
                ref=f[3],
                alt=f[4],
                qual=f[5],
                filt=f[6],
                fmt=f[8],
                smp=f[9],
            )
        )
    df = pd.DataFrame(rows, columns=["gene", "chrom", "pos", "vcf_id", "ref", "alt", "qual", "filt", "fmt", "smp"])
    return df.drop_duplicates(subset=["gene", "chrom", "pos", "ref", "alt"]).reset_index(drop=True)


def parse_sample(fmt: str, smp: str) -> dict[str, str]:
    return dict(zip(fmt.split(":"), smp.split(":")))


def zygosity(gt: str) -> str:
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return "NO_CALL"
    x, y = sorted(a)
    if x == y == "0":
        return "HOM_REF"
    if x == y:
        return "HOM_ALT"
    if x == "0":
        return "HET"
    return "HET_MULTI"          # two different non-reference alleles at one site


def allele_balance(ad: str, gt: str) -> float | None:
    """Fraction of reads supporting the called alternate allele(s)."""
    try:
        depths = [int(x) for x in ad.split(",")]
    except (ValueError, AttributeError):
        return None
    total = sum(depths)
    if total == 0:
        return None
    idx = {int(i) for i in gt.replace("|", "/").split("/") if i.isdigit() and i != "0"}
    return sum(depths[i] for i in idx if i < len(depths)) / total


def ab_plausible(zyg: str, ab: float | None) -> bool:
    """Is the observed allele balance consistent with the called zygosity?"""
    if ab is None:
        return True
    if zyg == "HET":
        return 0.20 <= ab <= 0.80
    if zyg == "HOM_ALT":
        return ab >= 0.85
    return True


def add_genotype_qc(df: pd.DataFrame, min_dp: int, min_gq: int) -> pd.DataFrame:
    """Append GT / zygosity / DP / GQ / allele balance / phase / QC columns."""
    out = []
    for r in df.itertuples(index=False):
        d = parse_sample(r.fmt, r.smp)
        gt = d.get("GT", "./.")
        zyg = zygosity(gt)
        try:
            dp = int(d.get("DP", "0"))
        except ValueError:
            dp = 0
        try:
            gq = int(float(d.get("GQ", "0")))
        except ValueError:
            gq = 0
        ab = allele_balance(d.get("AD", ""), gt)
        out.append(
            dict(
                GT=gt,
                zygosity=zyg,
                DP=dp,
                GQ=gq,
                AB=ab,
                PGT=d.get("PGT"),
                PID=d.get("PID"),
                ab_ok=ab_plausible(zyg, ab),
                qc_pass=(r.filt == "PASS") and dp >= min_dp and gq >= min_gq,
            )
        )
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(out)], axis=1)
