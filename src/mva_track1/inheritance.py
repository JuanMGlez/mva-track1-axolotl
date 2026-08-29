"""Inheritance-model logic.

The two questions the challenge asks are answered per gene, not per variant:

(a) is there a rare, damaging **homozygous** genotype, and
(b) are there two distinct rare damaging **heterozygous** variants in the same
    gene (a compound-heterozygous candidate)?

Phase is reported honestly. A pair is only ``in_trans`` when the VCF's own
physical-phasing fields (``PGT``/``PID``) place the two alleles on opposite
haplotypes; otherwise the verdict is ``unphased`` and parental segregation is
named as the resolving test.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from itertools import combinations

import pandas as pd

HIGH_OR_MODERATE = {"HIGH", "MODERATE"}
DAMAGING_TERMS = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
    "start_lost", "stop_lost", "transcript_ablation", "missense_variant",
    "inframe_deletion", "inframe_insertion", "protein_altering_variant",
}


def flag_candidates(df: pd.DataFrame, rarity_af_max: float) -> pd.DataFrame:
    """Add ``is_rare`` / ``is_damaging`` / ``candidate`` columns.

    A missing gnomAD frequency counts as rare: absence from 1.6 M alleles is
    evidence of rarity, not of missing data, for a callset of this quality.
    """
    out = df.copy()
    out["is_rare"] = out["af_joint"].isna() | (out["af_joint"] < rarity_af_max)
    out["is_damaging"] = out["impact"].isin(HIGH_OR_MODERATE) | out["most_severe"].isin(DAMAGING_TERMS)
    out["candidate"] = out["is_rare"] & out["is_damaging"] & out["qc_pass"]
    return out


def phase_status(a: pd.Series, b: pd.Series) -> tuple[str, str]:
    """``(status, explanation)`` for a candidate pair."""
    pa, pb = a.get("PGT"), b.get("PGT")
    ida, idb = a.get("PID"), b.get("PID")
    if pa and pb and ida and idb and ida == idb:
        if pa != pb:
            return "in_trans", f"physically phased in trans within phase set {ida} (PGT {pa} vs {pb})"
        return "in_cis", f"physically phased in cis within phase set {ida} (both PGT {pa})"
    dist = abs(int(a["pos"]) - int(b["pos"]))
    return (
        "unphased",
        f"no shared phase set; alleles are {dist:,} bp apart, beyond short-read phasing range. "
        "Parental Sanger segregation is the resolving test (and would activate PM3).",
    )


@dataclass
class GeneVerdict:
    gene: str
    panel: str
    inheritance: str
    n_variants: int
    n_qc_pass: int
    homozygous_candidate: bool
    compound_het_candidate: bool
    verdict: str
    detail: str


def evaluate_gene(df_gene: pd.DataFrame, gene: str, panel: str, inheritance: str) -> tuple[GeneVerdict, list[dict]]:
    cand = df_gene[df_gene["candidate"]]
    hom = cand[cand["zygosity"] == "HOM_ALT"]
    het = cand[cand["zygosity"].isin(["HET", "HET_MULTI"])]

    pairs: list[dict] = []
    for a, b in combinations([r for _, r in het.iterrows()], 2):
        status, why = phase_status(a, b)
        if status == "in_cis":
            continue                                   # same haplotype: not a biallelic genotype
        pairs.append(
            dict(
                gene=gene,
                a_pos=int(a["pos"]), a_ref=a["ref"], a_alt=a["alt"], a_csq=a["most_severe"], a_hgvsp=a.get("hgvsp"),
                b_pos=int(b["pos"]), b_ref=b["ref"], b_alt=b["alt"], b_csq=b["most_severe"], b_hgvsp=b.get("hgvsp"),
                phase=status, phase_note=why,
            )
        )

    if not hom.empty:
        verdict, detail = "homozygous candidate", f"{len(hom)} rare damaging homozygous genotype(s)"
    elif pairs:
        verdict, detail = "compound heterozygous candidate", f"{len(pairs)} candidate pair(s), phase: {pairs[0]['phase']}"
    else:
        n_common = int((~df_gene["is_rare"]).sum())
        n_noncoding = int((~df_gene["is_damaging"]).sum())
        verdict = "negative"
        detail = (
            f"no qualifying biallelic genotype; {n_common} variant(s) too common (AF >= threshold), "
            f"{n_noncoding} without coding/splice impact"
        )

    return (
        GeneVerdict(
            gene=gene, panel=panel, inheritance=inheritance,
            n_variants=len(df_gene), n_qc_pass=int(df_gene["qc_pass"].sum()),
            homozygous_candidate=not hom.empty,
            compound_het_candidate=bool(pairs),
            verdict=verdict, detail=detail,
        ),
        pairs,
    )


def evaluate(df: pd.DataFrame, genes: list[str], panel_of, inheritance_of) -> tuple[pd.DataFrame, pd.DataFrame]:
    verdicts, all_pairs = [], []
    for gene in genes:
        sub = df[df["gene"] == gene]
        if sub.empty:
            continue
        v, pairs = evaluate_gene(sub, gene, panel_of(gene), inheritance_of(gene))
        verdicts.append(asdict(v))
        all_pairs.extend(pairs)
    return pd.DataFrame(verdicts), pd.DataFrame(all_pairs)
