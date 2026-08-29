"""Offline unit tests for genotype, zygosity and QC derivation."""

import pandas as pd
import pytest

from mva_track1 import vcf_extract as vx


@pytest.mark.parametrize(
    "gt,expected",
    [("0/1", "HET"), ("1|0", "HET"), ("1/1", "HOM_ALT"), ("0/0", "HOM_REF"),
     ("1/2", "HET_MULTI"), ("./.", "NO_CALL"), (".", "NO_CALL")],
)
def test_zygosity(gt, expected):
    assert vx.zygosity(gt) == expected


def test_allele_balance_counts_only_the_called_alt():
    assert vx.allele_balance("21,25", "0/1") == pytest.approx(25 / 46)
    assert vx.allele_balance("15,13", "0/1") == pytest.approx(13 / 28)
    assert vx.allele_balance("0,40", "1/1") == pytest.approx(1.0)
    assert vx.allele_balance("10,5,20", "0/2") == pytest.approx(20 / 35)
    assert vx.allele_balance("0,0", "0/1") is None
    assert vx.allele_balance("", "0/1") is None


def test_ab_plausibility_gate():
    assert vx.ab_plausible("HET", 0.54)
    assert vx.ab_plausible("HET", 0.46)
    assert not vx.ab_plausible("HET", 0.05)       # likely artefact or contamination
    assert not vx.ab_plausible("HOM_ALT", 0.55)   # not really homozygous
    assert vx.ab_plausible("HOM_ALT", 0.98)
    assert vx.ab_plausible("HET", None)           # missing AD is not a failure


def test_add_genotype_qc_columns():
    df = pd.DataFrame([
        dict(gene="BUB1B", chrom="15", pos=40209701, vcf_id=".", ref="T", alt="G", qual="500",
             filt="PASS", fmt="GT:AD:DP:GQ", smp="0/1:21,25:46:99"),
        dict(gene="BUB1B", chrom="15", pos=1, vcf_id=".", ref="A", alt="C", qual="30",
             filt="PASS", fmt="GT:AD:DP:GQ", smp="0/1:2,1:3:9"),
    ])
    out = vx.add_genotype_qc(df, min_dp=10, min_gq=30)
    assert out.loc[0, "zygosity"] == "HET" and out.loc[0, "qc_pass"]
    assert out.loc[0, "AB"] == pytest.approx(25 / 46, abs=1e-3)
    assert not out.loc[1, "qc_pass"]              # DP 3, GQ 9
