"""End-to-end orchestration: VCF in, submission CSV + tables + figure out."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import annotate as ann
from . import hgvs, inheritance, submission, vcf_extract
from .config import Config

CODING_OR_SPLICE = inheritance.DAMAGING_TERMS | {
    "synonymous_variant", "splice_region_variant", "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant", "splice_donor_region_variant",
    "stop_retained_variant", "coding_sequence_variant", "incomplete_terminal_codon_variant",
}


def _log(msg: str) -> None:
    print(f"[mva] {msg}", flush=True)


def gene_coordinates(cfg: Config, ens: ann.EnsemblClient) -> pd.DataFrame:
    rows = []
    for gene in cfg.genes:
        rec = ens.lookup_gene(gene)
        rows.append(dict(
            gene=gene, panel=cfg.panel_of(gene), inheritance=cfg.inheritance_of(gene),
            ensembl_id=rec["id"], chrom=rec["seq_region_name"],
            start=int(rec["start"]), end=int(rec["end"]), strand=int(rec["strand"]),
            canonical_transcript=rec.get("canonical_transcript", "").split(".")[0],
            assembly=rec.get("assembly_name"), description=rec.get("description"),
        ))
    return pd.DataFrame(rows)


def cds_structure(coords: pd.DataFrame, ens: ann.EnsemblClient) -> dict[str, dict]:
    """``{gene: {transcript, strand, exons, cds_len, protein_len}}`` for the
    canonical transcript of each gene."""
    out = {}
    for r in coords.itertuples(index=False):
        feats = ens.cds_intervals(r.chrom, r.start, r.end)
        exons = sorted({(f["start"], f["end"]) for f in feats
                        if str(f.get("Parent", "")).split(".")[0] == r.canonical_transcript})
        if not exons:
            continue
        out[r.gene] = dict(
            transcript=r.canonical_transcript, strand=int(r.strand), exons=exons,
            cds_len=hgvs.cds_length(exons), protein_len=hgvs.protein_length(exons),
        )
    return out


def annotate_variants(df: pd.DataFrame, cfg: Config, ens, gnm, cvr, ucsc, cds: dict) -> pd.DataFrame:
    """VEP + gnomAD + ClinVar + phyloP + arithmetic HGVS, one row per variant."""
    recs = []
    for i, r in enumerate(df.itertuples(index=False), 1):
        if i % 25 == 0:
            _log(f"annotating {i}/{len(df)}")
        alt = str(r.alt).split(",")[0]
        region, allele = hgvs.to_vep_region_allele(r.chrom, r.pos, r.ref, alt)

        vep = ens.vep(region, allele)
        most_severe, impact, aa, ppos, sift, polyphen, hgvsc_vep, hgvsp = (None,) * 8
        if vep:
            v0 = vep[0]
            most_severe = v0.get("most_severe_consequence")
            tcs = v0.get("transcript_consequences") or []
            pick = next((t for t in tcs if t.get("gene_symbol") == r.gene and t.get("canonical")), None)
            pick = pick or next((t for t in tcs if t.get("gene_symbol") == r.gene), None)
            if pick:
                impact = pick.get("impact")
                aa = pick.get("amino_acids")
                ppos = pick.get("protein_start")
                sift = pick.get("sift_prediction")
                polyphen = pick.get("polyphen_prediction")
                hgvsc_vep = pick.get("hgvsc")
                hgvsp = pick.get("hgvsp")

        g = gnm.variant(r.chrom, r.pos, r.ref, alt) or {}
        def _freq(block):
            b = g.get(block) or {}
            ac, an = b.get("ac"), b.get("an")
            return (ac, an, (ac / an) if ac is not None and an else None, b.get("homozygote_count"))
        ac_e, an_e, af_e, hom_e = _freq("exome")
        ac_g, an_g, af_g, hom_g = _freq("genome")
        ac_j, an_j, af_j, hom_j = _freq("joint")

        cv = cvr.classify(r.chrom, r.pos, r.ref, alt)
        phylop = ucsc.phylop(r.chrom, r.pos)

        c = cds.get(r.gene)
        cpos = hgvs.coding_position(c["exons"], c["strand"], r.pos) if c else None

        recs.append(dict(
            most_severe=most_severe, impact=impact, aa=aa, protein_pos=ppos,
            sift=sift, polyphen=polyphen, hgvsc_vep=hgvsc_vep, hgvsp=hgvsp,
            c_pos=cpos, hgvsc_derived=(f"c.{cpos}{r.ref}>{alt}" if cpos else None),
            protein_pos_derived=(hgvs.protein_position(cpos) if cpos else None),
            rsid=(g.get("rsids") or [None])[0],
            ac_exome=ac_e, an_exome=an_e, af_exome=af_e, hom_exome=hom_e,
            ac_genome=ac_g, an_genome=an_g, af_genome=af_g, hom_genome=hom_g,
            ac_joint=ac_j, an_joint=an_j, af_joint=af_j, hom_joint=hom_j,
            phylop100way=phylop, **cv,
        ))
    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(recs)], axis=1)

    # Hard consistency check: our arithmetic HGVS must agree with VEP's own
    # protein position wherever both exist.
    both = out.dropna(subset=["protein_pos", "protein_pos_derived"])
    bad = both[both["protein_pos"].astype(int) != both["protein_pos_derived"].astype(int)]
    if not bad.empty:
        raise AssertionError(
            "derived protein position disagrees with VEP for:\n"
            + bad[["gene", "pos", "ref", "alt", "protein_pos", "protein_pos_derived"]].to_string(index=False)
        )
    return out


def run(vcf: str | Path, cfg: Config, outdir: str | Path, cache_dir: str | Path | None,
        use_submitted_epcr: bool = True, make_figure: bool = True) -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cache = ann.Cache(cache_dir)
    ens, gnm = ann.EnsemblClient(cache), ann.GnomadClient(cache)
    cvr, ucsc = ann.ClinVarClient(cache), ann.UcscClient(cache)
    paths: dict[str, Path] = {}

    _log(f"resolving GRCh38 coordinates for {len(cfg.genes)} genes")
    coords = gene_coordinates(cfg, ens)
    assert (coords["assembly"] == cfg.assembly).all(), f"expected {cfg.assembly}"

    _log("fetching gnomAD gene constraint")
    con = []
    for gene in cfg.genes:
        g = gnm.gene_constraint(gene) or {}
        con.append(dict(gene=gene, gene_id=g.get("gene_id"),
                        canonical_transcript=g.get("canonical_transcript_id"),
                        **(g.get("gnomad_constraint") or {})))
    constraint = pd.DataFrame(con)
    paths["constraint"] = outdir / "gene_constraint.csv"
    constraint.to_csv(paths["constraint"], index=False)

    _log("extracting variants from the VCF")
    ivs = [vcf_extract.Interval(r.chrom, r.start - cfg.padding_bp, r.end + cfg.padding_bp, r.gene)
           for r in coords.itertuples(index=False)]
    raw = vcf_extract.extract(vcf, ivs)
    v = vcf_extract.add_genotype_qc(raw, cfg.min_dp, cfg.min_gq)
    paths["all_variants"] = outdir / "gene_variants_qc.tsv"
    v.to_csv(paths["all_variants"], sep="\t", index=False)
    _log(f"{len(v)} variants extracted, {int(v.qc_pass.sum())} pass QC")

    _log("resolving canonical-transcript CDS structure")
    cds = cds_structure(coords, ens)
    (outdir / "cds_structure.json").write_text(json.dumps(
        {g: {k: val for k, val in d.items() if k != "exons"} | {"n_cds_exons": len(d["exons"])}
         for g, d in cds.items()}, indent=1))

    # Annotate everything in the primary panel(s) flagged annotate_all, plus
    # coding/splice-region variants everywhere else. Cheap first pass with VEP
    # on the whole set would be wasteful, so pre-filter on distance to CDS.
    def _near_cds(row) -> bool:
        c = cds.get(row.gene)
        if not c:
            return False
        w = cfg.splice_window_bp
        return any(s - w <= row.pos <= e + w for s, e in c["exons"])

    keep_all = {g for p in cfg.panels.values() if p.annotate_all for g in p.genes}
    mask = v.apply(lambda r: (r.gene in keep_all) or _near_cds(r), axis=1)
    sel = v[mask & v["filt"].eq("PASS")].reset_index(drop=True)
    _log(f"annotating {len(sel)} variants ({len(keep_all)} gene(s) annotated in full)")

    master = annotate_variants(sel, cfg, ens, gnm, cvr, ucsc, cds)
    master = inheritance.flag_candidates(master, cfg.rarity_af_max)
    master["panel"] = master["gene"].map(cfg.panel_of)
    paths["annotated"] = outdir / "annotated_candidate_variants.csv"
    master.to_csv(paths["annotated"], index=False)

    _log("evaluating inheritance models")
    verdicts, pairs = inheritance.evaluate(master, cfg.genes, cfg.panel_of, cfg.inheritance_of)
    paths["verdicts"] = outdir / "recessive_models.csv"
    verdicts.to_csv(paths["verdicts"], index=False)
    if not pairs.empty:
        paths["pairs"] = outdir / "compound_het_pairs.csv"
        pairs.to_csv(paths["pairs"], index=False)

    _log("writing submission")
    if use_submitted_epcr and cfg.submitted_epcr:
        rows = submission.rows_from_config(cfg.submitted_epcr, cfg.proband_id)
        name = "Axolotl_targeted-panel-vep-gnomad-acmg.csv"
    else:
        rows = []
        for _, r in master[master["candidate"]].iterrows():
            d = r.to_dict()
            rows.append(dict(proband_id=cfg.proband_id, chrom_1=f"chr{r['chrom']}", pos_1=int(r["pos"]),
                             ref_1=r["ref"], alt_1=str(r["alt"]).split(",")[0],
                             epcr=submission.heuristic_epcr(d),
                             finding_type="primary" if cfg.panel_of(r["gene"]) == "MVA" else "secondary",
                             notes=f"{r['gene']} {r['most_severe']} ({r['impact']}); "
                                   f"gnomAD AF {r['af_joint']}; ClinVar {r.get('clinvar_sig')}"))
        name = "heuristic_predictions.csv"
    paths["submission"] = submission.write_submission(rows, outdir / name)

    if make_figure:
        from .figure import make_figure as _mk
        prim = [g for g, p in ((g, cfg.panels[cfg.panel_of(g)]) for g in cfg.genes) if p.annotate_all]
        cand = master[master["candidate"] & master["gene"].isin(prim)].copy()
        if not cand.empty:
            cand["label"] = cand.apply(
                lambda r: (str(r["hgvsp"]).split(":")[-1] if pd.notna(r["hgvsp"])
                           else f"c.{r['c_pos']}"), axis=1)
            gene = cand["gene"].iloc[0]
            plen = cds[gene]["protein_len"]
            paths["figure"] = _mk(master, cand, plen, (766, plen, "kinase domain"),
                                  outdir / "mva_bub1b_candidates.png", tuple(prim))

    _log("done: " + ", ".join(f"{k}={p.name}" for k, p in paths.items()))
    return paths
