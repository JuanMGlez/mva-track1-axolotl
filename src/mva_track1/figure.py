"""The submitted two-panel figure, regenerated from the pipeline's own tables.

Panel (a) rarity vs. consequence severity for every variant in the primary
panel; panel (b) the candidate positions on the protein.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

FOCAL, GREY, BLUE, DARK = "#B02418", "#9aa0a6", "#2C6FA6", "#333333"

SEVERITY = [
    ("MODIFIER\n(intronic /\nregulatory)", {"MODIFIER"}),
    ("LOW\n(synonymous /\nsplice region)", {"LOW"}),
    ("MODERATE\n(missense)", {"MODERATE"}),
    ("HIGH\n(truncating)", {"HIGH"}),
]


def _style(base: float = 8.0) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": base,
        "axes.titlesize": base + 1,
        "axes.labelsize": base,
        "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1,
        "legend.fontsize": base - 1,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "figure.dpi": 200,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _panel_letter(ax, letter: str) -> None:
    ax.text(-0.085, 1.06, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def make_figure(
    variants: pd.DataFrame,
    candidates: pd.DataFrame,
    protein_length: int,
    domain: tuple[int, int, str],
    out_png: str | Path,
    panel_genes: tuple[str, ...] = ("BUB1B", "CEP57", "TRIP13"),
    af_floor: float = 3e-7,
) -> Path:
    """Render the figure. ``domain`` is ``(start_aa, end_aa, label)``."""
    _style()
    fig = plt.figure(figsize=(180 / 25.4, 92 / 25.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.80], hspace=0.72,
                          left=0.10, right=0.985, top=0.905, bottom=0.095)

    # ---- panel a: rarity vs consequence severity -------------------------
    ax = fig.add_subplot(gs[0])
    d = variants[variants["gene"].isin(panel_genes)].copy()
    d["y"] = [next((i for i, (_, s) in enumerate(SEVERITY) if str(imp) in s), 0) for imp in d["impact"]]
    d["af_plot"] = d["af_joint"].fillna(af_floor).clip(lower=af_floor)

    markers = {"BUB1B": "o", "CEP57": "s", "TRIP13": "^"}
    rng = np.random.default_rng(0)
    for gene, m in markers.items():
        sub = d[d["gene"] == gene]
        if sub.empty:
            continue
        jitter = rng.uniform(-0.13, 0.13, len(sub))
        is_cand = sub["candidate"].to_numpy(dtype=bool)
        ax.scatter(sub["af_plot"][~is_cand], sub["y"][~is_cand] + jitter[~is_cand],
                   marker=m, s=22, facecolor="none", edgecolor=GREY, linewidth=0.8,
                   label=gene, zorder=3)
        ax.scatter(sub["af_plot"][is_cand], sub["y"][is_cand] + jitter[is_cand],
                   marker=m, s=46, color=FOCAL, edgecolor="white", linewidth=0.6, zorder=5)

    ax.axvline(0.01, color=DARK, ls=(0, (4, 3)), lw=0.8, zorder=1)
    ax.text(0.01, 3.62, " rarity threshold, AF = 1%", fontsize=6.5, color=DARK, va="top", ha="left")
    ax.set_xscale("log")
    ax.set_xlim(af_floor / 2, 1.6)
    ax.set_ylim(-0.55, 3.75)
    ax.set_yticks(range(len(SEVERITY)))
    ax.set_yticklabels([s[0] for s in SEVERITY], fontsize=6.2, linespacing=0.95)
    ax.set_xlabel("gnomAD v4 allele frequency (joint exomes + genomes; leftmost column = absent)")
    ax.set_title("Only two of the variants in the three MVA genes are both rare and protein-altering — both in $\\it{BUB1B}$",
                 loc="left", pad=8)
    ax.legend(loc="upper left", frameon=False, handletextpad=0.25, borderpad=0.1,
              bbox_to_anchor=(0.005, 1.005), labelspacing=0.25)
    _panel_letter(ax, "a")

    for _, r in candidates.iterrows():
        af = r["af_joint"] if pd.notna(r.get("af_joint")) else af_floor
        y = next((i for i, (_, s) in enumerate(SEVERITY) if str(r["impact"]) in s), 0)
        ax.annotate(r["label"], xy=(max(af, af_floor), y), xytext=(0, -14),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=6.8, color=FOCAL,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=FOCAL, shrinkA=1, shrinkB=3))

    # ---- panel b: protein schematic --------------------------------------
    ax2 = fig.add_subplot(gs[1])
    L = protein_length
    ax2.add_patch(Rectangle((1, 0.42), L, 0.16, facecolor="#eceff1",
                            edgecolor=DARK, lw=0.7, zorder=2))
    ds, de, dlabel = domain
    ax2.add_patch(Rectangle((ds, 0.42), de - ds, 0.16, facecolor=BLUE, alpha=0.30,
                            edgecolor=BLUE, lw=0.7, hatch="///", zorder=3))
    ax2.text((ds + de) / 2, 0.34, dlabel, fontsize=6.5, color=BLUE, ha="center", va="top")

    for _, r in candidates.iterrows():
        if pd.isna(r.get("protein_pos")):
            continue
        x = float(r["protein_pos"])
        ax2.plot([x, x], [0.42, 0.72], color=FOCAL, lw=1.0, zorder=4)
        ax2.plot([x], [0.72], marker="v", color=FOCAL, ms=4.5, zorder=5)
        ax2.text(x, 0.78, r["label"], fontsize=6.8, color=FOCAL, ha="center", va="bottom")

    ax2.annotate("", xy=(L, 0.24), xytext=(float(candidates["protein_pos"].min()), 0.24),
                 arrowprops=dict(arrowstyle="<->", lw=0.6, color=DARK))
    lost = int(L - candidates["protein_pos"].min())
    ax2.text(L, 0.17, f"{lost} residues lost by the truncating allele", fontsize=6.5,
             color=DARK, ha="right", va="top")

    ax2.set_xlim(-30, L + 40)
    ax2.set_ylim(0.05, 1.02)
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.set_xlabel("BubR1 residue")
    ax2.set_title("The two alleles hit the same essential domain from opposite directions", loc="left", pad=6)
    _panel_letter(ax2, "b")

    out_png = Path(out_png)
    fig.savefig(out_png)
    plt.close(fig)
    return out_png
