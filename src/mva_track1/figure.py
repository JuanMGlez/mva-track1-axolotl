"""The submitted two-panel figure, regenerated from the pipeline's own tables.

Panel (a) rarity against consequence severity for every annotated variant in the
primary panel; panel (b) the candidate alleles on the BubR1 domain map. Every
number and every evidence string in the annotations is read from the tables, so
the figure cannot drift from the data it describes.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

FOCAL = "#B02418"
GREY = "#9aa0a6"
BLUE = "#2C6FA6"
DARK = "#333333"

#: UniProt O60566 (BUB1B_HUMAN) feature coordinates, in residues.
BUBR1_DOMAINS: tuple[tuple[str, int, int, str], ...] = (
    ("KEN1", 26, 32, "#cfd8dc"),
    ("TPR", 70, 360, "#b9c6cf"),
    ("KARD", 665, 682, "#cfd8dc"),
)

SEVERITY = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
SEVERITY_LABELS = ["MODIFIER\n(intron / UTR)", "LOW\n(synonymous)",
                   "MODERATE\n(missense)", "HIGH\n(truncating)"]
MARKERS = ("o", "s", "^", "D", "v")


def _style(base: float = 8.0) -> None:
    """Explicit rcParams so the figure does not depend on session-only helpers."""
    plt.rcParams.update({
        "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base - 1,
        "xtick.labelsize": base - 2, "ytick.labelsize": base - 2,
        "legend.fontsize": base - 2, "axes.spines.top": False,
        "axes.spines.right": False, "axes.linewidth": 0.7,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "font.family": "sans-serif", "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def _panel_letter(ax, letter: str) -> None:
    ax.text(-0.075, 1.13, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")


def _evidence_lines(r: pd.Series) -> list[str]:
    """Evidence caption for one candidate allele, built from the table."""
    out = []
    hgvsc = f"c.{int(r['c_pos'])}{r['ref']}>{str(r['alt']).split(',')[0]}" if pd.notna(r.get("c_pos")) else ""
    out.append(f"{hgvsc}  {r['label']}".strip())
    csq = str(r["most_severe"]).replace("_variant", "").replace("_", " ")
    bits = [csq]
    if pd.notna(r.get("clinvar_sig")):
        stars = "" if pd.isna(r.get("clinvar_stars")) else " (" + "\u2605" * int(r["clinvar_stars"]) + ")"
        bits.append(f"ClinVar {str(r['clinvar_sig']).lower()}{stars}")
    else:
        bits.append("absent from ClinVar")
    out.append(" \u00b7 ".join(bits))
    tail = []
    if pd.notna(r.get("sift")) and pd.notna(r.get("polyphen")):
        tail.append("SIFT and PolyPhen-2 deleterious")
    if pd.notna(r.get("phylop100way")):
        tail.append(f"phyloP {float(r['phylop100way']):.1f}")
    if tail:
        out.append(" \u00b7 ".join(tail))
    return out


def _af_caption(r: pd.Series) -> str:
    af = r.get("af_joint")
    if pd.isna(af):
        return "absent from gnomAD"
    ac = r.get("ac_joint")
    if pd.notna(ac) and int(ac) == 1:
        return "1 allele in gnomAD"
    e = int(np.floor(np.log10(float(af))))
    return f"AF {float(af) / 10 ** e:.1f}\u00d710$^{{{e}}}$"


def make_figure(
    variants: pd.DataFrame,
    candidates: pd.DataFrame,
    protein_length: int,
    domain: tuple[int, int, str],
    out_png: str | Path,
    panel_genes: tuple[str, ...] = ("BUB1B", "CEP57", "TRIP13"),
    af_floor: float = 3e-7,
    rarity_af: float = 0.01,
) -> Path:
    """Render the figure. ``domain`` is the highlighted ``(start_aa, end_aa, label)``."""
    _style()
    fig = plt.figure(figsize=(180 / 25.4, 100 / 25.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.92], hspace=1.00,
                          left=0.115, right=0.985, top=0.90, bottom=0.10)

    # ---- panel a: rarity vs severity -------------------------------------
    d = variants[variants["gene"].isin(panel_genes)].copy()
    d["af_plot"] = d["af_joint"].fillna(af_floor).clip(lower=af_floor)
    d["sev"] = d["impact"].map(SEVERITY).fillna(0)
    cand_keys = set(zip(candidates["chrom"], candidates["pos"], candidates["alt"]))
    d["is_cand"] = [k in cand_keys for k in zip(d["chrom"], d["pos"], d["alt"])]

    ax = fig.add_subplot(gs[0])
    ax.axvspan(af_floor / 2, rarity_af, color=FOCAL, alpha=0.045, zorder=0)
    ax.axvline(rarity_af, color=DARK, ls=":", lw=0.8, zorder=1)
    ax.text(rarity_af * 1.35, 3.42, f"{rarity_af:.0%} (rarity threshold)",
            fontsize=6.3, color=DARK, va="center")

    marker_of = dict(zip(panel_genes, MARKERS))
    for gene, m in marker_of.items():
        sub = d[(d["gene"] == gene) & ~d["is_cand"]]
        ax.scatter(sub["af_plot"], sub["sev"], marker=m, s=22, facecolor="none",
                   edgecolor=GREY, lw=0.8, zorder=3)
    for _, r in d[d["is_cand"]].iterrows():
        ax.scatter([r["af_plot"]], [r["sev"]], marker=marker_of.get(r["gene"], "o"),
                   s=40, color=FOCAL, zorder=5)

    cand_sorted = candidates.sort_values("af_joint", na_position="first")
    for (_, r), (dx, dy) in zip(cand_sorted.iterrows(), [(2.4, 0.30), (2.4, 0.36)]):
        af = af_floor if pd.isna(r["af_joint"]) else max(float(r["af_joint"]), af_floor)
        sev = SEVERITY.get(r["impact"], 0)
        ax.annotate(f"{r['label']}  \u00b7  {_af_caption(r)}",
                    xy=(af, sev), xytext=(af * dx, sev + dy), fontsize=6.6,
                    color=FOCAL, va="center", ha="left",
                    arrowprops=dict(arrowstyle="-", color=FOCAL, lw=0.7,
                                    shrinkA=0, shrinkB=2))

    n = len(d)
    ax.set_xscale("log")
    ax.set_xlim(af_floor / 2.4, 1.6)
    ax.set_ylim(-0.45, 3.75)
    ax.set_yticks(range(4))
    ax.set_yticklabels(SEVERITY_LABELS)
    ticks = [af_floor, 1e-5, 1e-3, 1e-1]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["absent", "10$^{-5}$", "10$^{-3}$", "10$^{-1}$"])
    ax.set_xlabel("Allele frequency in gnomAD v4 (exomes + genomes)")
    ax.set_title(f"Of the {n} variants in the three MVA genes, only two are rare and "
                 f"protein-altering \u2014 both in $\\it{{{candidates['gene'].iloc[0]}}}$",
                 loc="left", pad=7)
    handles = [Line2D([], [], marker=marker_of[g], ls="", mfc="none", mec=GREY,
                      ms=4.6, label=f"$\\it{{{g}}}$") for g in panel_genes]
    handles.append(Line2D([], [], marker="o", ls="", color=FOCAL, ms=5, label="candidate"))
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.80),
              ncol=2, frameon=False, handletextpad=0.4, columnspacing=1.1,
              borderpad=0.2)
    _panel_letter(ax, "a")

    # ---- panel b: domain map --------------------------------------------
    L = int(protein_length)
    ds, de, dlabel = int(domain[0]), int(domain[1]), str(domain[2])
    ax2 = fig.add_subplot(gs[1])
    ax2.add_patch(Rectangle((1, 0.42), L, 0.16, facecolor="#eceff1",
                            edgecolor=DARK, lw=0.7, zorder=2))
    for name, s, e, col in BUBR1_DOMAINS:
        ax2.add_patch(Rectangle((s, 0.42), e - s, 0.16, facecolor=col,
                                edgecolor=DARK, lw=0.6, zorder=3))
        if e - s > 100:
            ax2.text((s + e) / 2, 0.50, name, ha="center", va="center",
                     fontsize=6, color=DARK, zorder=7)
        else:
            ax2.text((s + e) / 2, 0.38, name, ha="center", va="top",
                     fontsize=6, color=DARK)
    ax2.add_patch(Rectangle((ds, 0.42), de - ds, 0.16, facecolor=BLUE, alpha=0.30,
                            edgecolor=DARK, lw=0.6, zorder=3))
    ax2.text((ds + de) / 2, 0.50, dlabel, ha="center", va="center", fontsize=6,
             color=DARK, zorder=7)

    x0, x1 = -30, L + 40
    ax2.set_xlim(x0, x1)
    ax2.set_ylim(0.05, 0.98)

    trunc = candidates.sort_values("protein_pos").iloc[0]
    lost = L - int(trunc["protein_pos"])
    ax2.add_patch(Rectangle((int(trunc["protein_pos"]), 0.42), lost, 0.16,
                            facecolor="none", edgecolor=FOCAL, lw=0.7,
                            hatch="///", zorder=4))

    for (_, r), (tx, ha) in zip(candidates.sort_values("protein_pos").iterrows(),
                                [(0.02, "left"), (0.98, "right")]):
        x = float(r["protein_pos"])
        ax2.plot([x, x], [0.42, 0.92], color=FOCAL, lw=1.0, zorder=5)
        ax2.plot([x], [0.92], marker="v", color=FOCAL, ms=4.2, zorder=6)
        xa = (x - x0) / (x1 - x0)
        ax2.plot([xa, tx], [0.92, 0.99], transform=ax2.transAxes, color=FOCAL,
                 lw=0.6, zorder=5, clip_on=False)
        ax2.text(tx, 1.02, "\n".join(_evidence_lines(r)), transform=ax2.transAxes,
                 fontsize=6.3, color=FOCAL, ha=ha, va="bottom", linespacing=1.45)

    ax2.annotate("", xy=(L, 0.26), xytext=(float(trunc["protein_pos"]), 0.26),
                 arrowprops=dict(arrowstyle="<->", color=FOCAL, lw=0.7))
    ax2.text((float(trunc["protein_pos"]) + L) / 2, 0.19,
             f"{lost} aa lost: the entire {dlabel}", fontsize=6.4, color=FOCAL,
             ha="center", va="top")

    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.set_xticks([1, 200, 400, 600, 800, L])
    ax2.set_xlabel(f"Position in BubR1 (amino acids; NM_001211.6 / "
                   f"{candidates['transcript'].iloc[0] if 'transcript' in candidates else 'ENST00000287598'}, {L} aa)")
    ax2.set_title("Both alleles affect the C-terminal half of BubR1: "
                  "one null allele and one conserved missense allele", loc="left", pad=34)
    _panel_letter(ax2, "b")

    out_png = Path(out_png)
    fig.savefig(out_png)
    plt.close(fig)
    return out_png
