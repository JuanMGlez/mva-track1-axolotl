"""Configuration loading: gene panels, thresholds, submitted EPCR values."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class Panel:
    name: str
    label: str
    genes: tuple[str, ...]
    inheritance: dict[str, str]
    #: If True, every variant in the gene interval is annotated, not just
    #: coding/splice-region ones. Used for the primary-hypothesis panel so that
    #: deep-intronic rare variants are still visible.
    annotate_all: bool = False


@dataclass(frozen=True)
class Config:
    assembly: str
    proband_id: str
    padding_bp: int
    rarity_af_max: float
    min_dp: int
    min_gq: int
    splice_window_bp: int
    hpo_terms: dict[str, str]
    panels: dict[str, Panel]
    submitted_epcr: list[dict[str, Any]] = field(default_factory=list)

    @property
    def genes(self) -> list[str]:
        out: list[str] = []
        for p in self.panels.values():
            out.extend(p.genes)
        return out

    def panel_of(self, gene: str) -> str:
        for name, p in self.panels.items():
            if gene in p.genes:
                return name
        raise KeyError(gene)

    def inheritance_of(self, gene: str) -> str:
        p = self.panels[self.panel_of(gene)]
        return p.inheritance.get(gene, "AR")


def load_config(path: str | Path) -> Config:
    raw = yaml.safe_load(Path(path).read_text())
    panels = {}
    for name, blk in raw["panels"].items():
        panels[name] = Panel(
            name=name,
            label=blk["label"],
            genes=tuple(blk["genes"]),
            inheritance=dict(blk.get("inheritance", {})),
            annotate_all=bool(blk.get("annotate_all", False)),
        )
    return Config(
        assembly=raw["assembly"],
        proband_id=raw["proband_id"],
        padding_bp=int(raw["padding_bp"]),
        rarity_af_max=float(raw["rarity_af_max"]),
        min_dp=int(raw["qc"]["min_dp"]),
        min_gq=int(raw["qc"]["min_gq"]),
        splice_window_bp=int(raw.get("splice_window_bp", 8)),
        hpo_terms=dict(raw.get("hpo_terms", {})),
        panels=panels,
        submitted_epcr=list(raw.get("submitted_epcr", [])),
    )
