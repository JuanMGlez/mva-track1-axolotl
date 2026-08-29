"""Coordinate normalisation and HGVS derivation.

Two things live here that we deliberately do NOT delegate to an annotation
service:

* :func:`to_vep_region_allele` converts a VCF-style ``(pos, ref, alt)`` into the
  ``region`` / ``allele`` pair the Ensembl VEP REST endpoint expects, using
  standard VCF left-trimming.
* :func:`coding_position` recomputes the HGVS ``c.`` position arithmetically
  from the canonical transcript's own CDS exon structure. We cross-check the
  implied protein position against the one VEP returns independently; a
  mismatch is a hard error rather than a warning.
"""

from __future__ import annotations

Interval = tuple[int, int]


def trim(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """VCF left-trim. Returns ``(pos, ref, alt)`` with ``''`` for a pure
    insertion/deletion allele (the anchor base removed)."""
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:      # deletion
        return pos + 1, ref[1:], ""
    if len(alt) > 1 and len(ref) == 1 and ref[0] == alt[0]:      # insertion
        return pos + 1, "", alt[1:]
    return pos, ref, alt


def to_vep_region_allele(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, str]:
    """Ensembl VEP ``/vep/human/region/{region}/{allele}`` representation."""
    c = chrom.replace("chr", "")
    p, r, a = trim(pos, ref, alt)
    if r and not a:                                  # deletion
        return f"{c}:{p}-{p + len(r) - 1}:1", "-"
    if a and not r:                                  # insertion
        return f"{c}:{p}-{p - 1}:1", a
    return f"{c}:{p}-{p + len(r) - 1}:1", a          # substitution / MNV / complex


def cds_length(exons: list[Interval]) -> int:
    return sum(e - s + 1 for s, e in exons)


def coding_position(exons: list[Interval], strand: int, gpos: int) -> int | None:
    """HGVS ``c.`` position of a genomic coordinate, or ``None`` if not in CDS.

    ``exons`` are genomic CDS intervals (any order, 1-based inclusive) of a
    single transcript; ``strand`` is ``+1`` or ``-1``.
    """
    ex = sorted(exons, reverse=strand < 0)
    cum = 0
    for s, e in ex:
        if s <= gpos <= e:
            return cum + (gpos - s + 1 if strand > 0 else e - gpos + 1)
        cum += e - s + 1
    return None


def protein_position(cpos: int) -> int:
    return (cpos - 1) // 3 + 1


def protein_length(exons: list[Interval]) -> int:
    """Residues excluding the stop codon."""
    return cds_length(exons) // 3 - 1
