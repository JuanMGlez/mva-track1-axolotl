"""Offline unit tests for the recessive-model logic and EPCR ranking."""

import pandas as pd

from mva_track1 import inheritance as inh
from mva_track1 import submission


def _v(**kw):
    base = dict(gene="BUB1B", chrom="15", pos=1, ref="A", alt="G", zygosity="HET",
                impact="HIGH", most_severe="stop_gained", af_joint=1e-5, qc_pass=True,
                PGT=None, PID=None, hgvsp=None)
    base.update(kw)
    return base


def test_absent_from_gnomad_counts_as_rare():
    df = inh.flag_candidates(pd.DataFrame([_v(af_joint=None)]), 0.01)
    assert bool(df.loc[0, "is_rare"])


def test_common_variant_is_not_a_candidate_even_if_severe():
    # The ATP6V1B1 start-loss, rs11681642: AF 43%, ClinVar Benign 2 stars.
    df = inh.flag_candidates(pd.DataFrame([_v(gene="ATP6V1B1", most_severe="start_lost", af_joint=0.43)]), 0.01)
    assert not bool(df.loc[0, "candidate"])


def test_two_rare_damaging_hets_give_a_compound_het_candidate():
    df = inh.flag_candidates(pd.DataFrame([
        _v(pos=40209701, ref="T", alt="G", af_joint=7.4e-5),
        _v(pos=40220612, ref="T", alt="G", impact="MODERATE", most_severe="missense_variant", af_joint=6.8e-7),
    ]), 0.01)
    v, pairs = inh.evaluate(df, ["BUB1B"], lambda g: "MVA", lambda g: "AR")
    assert v.loc[0, "verdict"] == "compound heterozygous candidate"
    assert len(pairs) == 1
    assert pairs.loc[0, "phase"] == "unphased"
    assert "10,911 bp apart" in pairs.loc[0, "phase_note"]


def test_cis_pair_is_rejected_not_reported():
    df = inh.flag_candidates(pd.DataFrame([
        _v(pos=100, PGT="0|1", PID="ps1"),
        _v(pos=200, PGT="0|1", PID="ps1"),
    ]), 0.01)
    v, pairs = inh.evaluate(df, ["BUB1B"], lambda g: "MVA", lambda g: "AR")
    assert pairs.empty
    assert v.loc[0, "verdict"] == "negative"


def test_trans_pair_is_recognised():
    df = inh.flag_candidates(pd.DataFrame([
        _v(pos=100, PGT="0|1", PID="ps1"),
        _v(pos=200, PGT="1|0", PID="ps1"),
    ]), 0.01)
    _, pairs = inh.evaluate(df, ["BUB1B"], lambda g: "MVA", lambda g: "AR")
    assert pairs.loc[0, "phase"] == "in_trans"


def test_homozygous_takes_precedence_over_pairing():
    df = inh.flag_candidates(pd.DataFrame([_v(zygosity="HOM_ALT"), _v(pos=2)]), 0.01)
    v, _ = inh.evaluate(df, ["BUB1B"], lambda g: "MVA", lambda g: "AR")
    assert v.loc[0, "verdict"] == "homozygous candidate"


def test_heuristic_epcr_ranks_the_submitted_alleles_in_the_right_order():
    truncating = submission.heuristic_epcr(dict(
        clinvar_sig="Pathogenic/Likely pathogenic", clinvar_stars=2, impact="HIGH",
        af_joint=7.4e-5, phenotype_match=True))
    intronic = submission.heuristic_epcr(dict(impact="MODIFIER", af_joint=None))
    assert truncating > intronic
    assert 0 < intronic <= truncating <= 0.99


def test_submission_is_written_in_descending_epcr_order(tmp_path):
    rows = [dict(proband_id="X", chrom_1="chr15", pos_1=1, ref_1="A", alt_1="G", epcr=0.1,
                 finding_type="primary", notes="low"),
            dict(proband_id="X", chrom_1="chr15", pos_1=2, ref_1="A", alt_1="G", epcr=0.85,
                 finding_type="primary", notes="high")]
    p = submission.write_submission(rows, tmp_path / "s.csv")
    lines = p.read_text().strip().split("\n")
    assert lines[0].startswith("proband_id,chrom_1")
    assert ",0.85," in lines[1] and ",0.1," in lines[2]
    assert lines[1].count(",") == lines[0].count(",")     # blank allele-2 columns kept
