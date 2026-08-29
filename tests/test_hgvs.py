"""Offline unit tests for coordinate handling and HGVS derivation.

The BUB1B expectations are the two submitted candidate alleles, checked against
the CDS exon structure of ENST00000287598 (plus strand, 23 coding exons,
CDS 3,153 nt, BubR1 1,050 aa) and against ClinVar's own title for the
truncating allele, ``NM_001211.6(BUB1B):c.2210T>G (p.Leu737Ter)``.
"""

import json
from pathlib import Path

import pytest

from mva_track1 import hgvs

FIX = json.loads((Path(__file__).parent / "data" / "bub1b_cds.json").read_text())
EXONS = [tuple(e) for e in FIX["exons"]]
STRAND = FIX["strand"]


def test_cds_and_protein_length():
    assert hgvs.cds_length(EXONS) == 3153
    assert hgvs.protein_length(EXONS) == 1050


@pytest.mark.parametrize("gpos,expected_c", [(40209701, 2210), (40220612, 3006)])
def test_coding_position_matches_clinvar_and_vep(gpos, expected_c):
    assert hgvs.coding_position(EXONS, STRAND, gpos) == expected_c


@pytest.mark.parametrize("gpos,expected_p", [(40209701, 737), (40220612, 1002)])
def test_protein_position(gpos, expected_p):
    c = hgvs.coding_position(EXONS, STRAND, gpos)
    assert hgvs.protein_position(c) == expected_p


def test_intronic_position_is_not_in_cds():
    # chr15:40,216,470 A>G, the deep-intronic variant we kept as a low-ranked
    # alternative second allele: it must not resolve to a coding position.
    assert hgvs.coding_position(EXONS, STRAND, 40216470) is None


def test_minus_strand_numbering_runs_the_other_way():
    exons = [(100, 109), (200, 209)]        # 20 nt CDS
    assert hgvs.coding_position(exons, +1, 100) == 1
    assert hgvs.coding_position(exons, +1, 209) == 20
    assert hgvs.coding_position(exons, -1, 209) == 1
    assert hgvs.coding_position(exons, -1, 100) == 20


@pytest.mark.parametrize(
    "pos,ref,alt,region,allele",
    [
        (40209701, "T", "G", "15:40209701-40209701:1", "G"),      # SNV
        (1000, "GT", "G", "15:1001-1001:1", "-"),                 # 1 bp deletion
        (1000, "G", "GT", "15:1001-1000:1", "T"),                 # 1 bp insertion
        (1000, "GAAT", "G", "15:1001-1003:1", "-"),               # 3 bp deletion
        (1000, "AT", "GC", "15:1000-1001:1", "GC"),               # MNV, nothing to trim
    ],
)
def test_vep_region_allele(pos, ref, alt, region, allele):
    assert hgvs.to_vep_region_allele("chr15", pos, ref, alt) == (region, allele)
