"""Clients for the four public annotation sources used by the pipeline.

All of them are free and keyless. Every response is cached on disk under
``--cache-dir`` so a re-run costs no network calls, which is what makes the
submission reproducible rather than merely repeatable.

* **Ensembl REST** — gene coordinates, canonical transcript, CDS exon structure,
  VEP consequence with SIFT and PolyPhen-2.
* **gnomAD v4 GraphQL** — per-allele AC/AN/homozygote count for exomes, genomes
  and the joint callset; gene-level constraint.
* **NCBI E-utilities (ClinVar)** — variation ID, germline classification, review
  status; star rating derived from the review status with the published mapping.
* **UCSC Genome Browser API** — phyloP100way basewise conservation.

Set ``NCBI_EMAIL`` (and optionally ``NCBI_API_KEY``) to identify yourself to
E-utilities as NCBI requests; the client omits the field when unset rather than
sending a placeholder.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import requests

ENSEMBL = "https://rest.ensembl.org"
GNOMAD = "https://gnomad.broadinstitute.org/api"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UCSC = "https://api.genome.ucsc.edu"

#: ClinVar review status -> gold stars.
#: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
}


class Cache:
    def __init__(self, path: str | Path | None):
        self.dir = Path(path) if path else None
        if self.dir:
            self.dir.mkdir(parents=True, exist_ok=True)

    def _p(self, key: str) -> Path | None:
        if not self.dir:
            return None
        return self.dir / (hashlib.sha1(key.encode()).hexdigest() + ".json")

    def get(self, key: str) -> Any | None:
        p = self._p(key)
        if p and p.exists():
            return json.loads(p.read_text())
        return None

    def put(self, key: str, value: Any) -> None:
        p = self._p(key)
        if p:
            p.write_text(json.dumps(value))


class Client:
    """Shared HTTP behaviour: caching, retry with backoff, polite pacing."""

    def __init__(self, cache: Cache, min_interval_s: float = 0.15, retries: int = 4, timeout: int = 60):
        self.cache, self.retries, self.timeout = cache, retries, timeout
        self.min_interval_s = min_interval_s
        self._last = 0.0
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "mva-track1-pipeline/1.0"

    def _pace(self) -> None:
        wait = self.min_interval_s - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()

    def _request(self, key: str, method: str, url: str, **kw) -> Any:
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        last: Exception | None = None
        for attempt in range(self.retries):
            try:
                self._pace()
                r = self.session.request(method, url, timeout=self.timeout, **kw)
                if r.status_code in (429, 500, 502, 503, 504):
                    raise requests.HTTPError(f"{r.status_code} from {url}")
                r.raise_for_status()
                out = r.json()
                self.cache.put(key, out)
                return out
            except Exception as exc:                       # noqa: BLE001 - retried below
                last = exc
                time.sleep(2 ** attempt)
        raise RuntimeError(f"request failed after {self.retries} attempts: {url}") from last


class EnsemblClient(Client):
    def lookup_gene(self, symbol: str, species: str = "homo_sapiens") -> dict:
        url = f"{ENSEMBL}/lookup/symbol/{species}/{symbol}"
        return self._request(
            f"ens:lookup:{species}:{symbol}", "GET", url,
            params={"expand": 0}, headers={"Content-Type": "application/json"},
        )

    def cds_intervals(self, chrom: str, start: int, end: int) -> list[dict]:
        url = f"{ENSEMBL}/overlap/region/homo_sapiens/{chrom}:{start}-{end}"
        return self._request(
            f"ens:cds:{chrom}:{start}:{end}", "GET", url,
            params={"feature": "cds"}, headers={"Content-Type": "application/json"},
        )

    def vep(self, region: str, allele: str) -> list[dict]:
        url = f"{ENSEMBL}/vep/human/region/{region}/{allele}"
        return self._request(
            f"ens:vep:{region}:{allele}", "GET", url,
            params={"canonical": 1, "hgvs": 1, "numbers": 1, "domains": 0,
                    "variant_class": 1, "sift": "b", "polyphen": "b"},
            headers={"Content-Type": "application/json"},
        )


class GnomadClient(Client):
    VARIANT_Q = """
    query V($v: String!, $d: DatasetId!) {
      variant(variantId: $v, dataset: $d) {
        variant_id rsids
        exome  { ac an homozygote_count }
        genome { ac an homozygote_count }
        joint  { ac an homozygote_count }
      }
    }"""

    GENE_Q = """
    query G($g: String!) {
      gene(gene_symbol: $g, reference_genome: GRCh38) {
        gene_id canonical_transcript_id chrom start stop
        gnomad_constraint { pLI oe_lof oe_lof_upper oe_mis obs_lof exp_lof mis_z lof_z }
      }
    }"""

    #: gnomAD reports an absent variant as a GraphQL *error* with a null data
    #: field. For us that is a meaningful result (absence from 1.6 M alleles is
    #: rarity evidence), not a failure, so it is filtered out here rather than
    #: raised. Any other error message still raises.
    NOT_FOUND = "variant not found"

    def _gql(self, key: str, query: str, variables: dict) -> dict:
        body = self._request(key, "POST", GNOMAD, json={"query": query, "variables": variables})
        errors = [e for e in (body.get("errors") or [])
                  if self.NOT_FOUND not in str(e.get("message", "")).lower()]
        if errors:
            raise RuntimeError(f"gnomAD GraphQL error: {errors}")
        return body.get("data") or {}

    def variant(self, chrom: str, pos: int, ref: str, alt: str, dataset: str = "gnomad_r4") -> dict | None:
        """Allele record, or ``None`` when the allele is absent from gnomAD."""
        vid = f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"
        data = self._gql(f"gnomad:var:{dataset}:{vid}", self.VARIANT_Q, {"v": vid, "d": dataset})
        return data.get("variant")

    def gene_constraint(self, symbol: str) -> dict | None:
        data = self._gql(f"gnomad:gene:{symbol}", self.GENE_Q, {"g": symbol})
        return data.get("gene")


class ClinVarClient(Client):
    def _eutils_params(self, extra: dict) -> dict:
        p = {"db": "clinvar", "retmode": "json", **extra}
        if os.environ.get("NCBI_EMAIL"):
            p["email"] = os.environ["NCBI_EMAIL"]
        if os.environ.get("NCBI_API_KEY"):
            p["api_key"] = os.environ["NCBI_API_KEY"]
        return p

    def uids_at_position(self, chrom: str, pos: int) -> list[str]:
        term = f"{chrom.replace('chr', '')}[chr] AND {pos}[chrpos38]"
        body = self._request(
            f"cv:search:{chrom}:{pos}", "GET", f"{EUTILS}/esearch.fcgi",
            params=self._eutils_params({"term": term, "retmax": 20}),
        )
        return list(body.get("esearchresult", {}).get("idlist", []))

    def summary(self, uid: str) -> dict:
        body = self._request(
            f"cv:sum:{uid}", "GET", f"{EUTILS}/esummary.fcgi",
            params=self._eutils_params({"id": uid}),
        )
        return body.get("result", {}).get(uid, {})

    def classify(self, chrom: str, pos: int, ref: str, alt: str) -> dict:
        """Best-matching ClinVar record for one allele, or empty dict."""
        for uid in self.uids_at_position(chrom, pos):
            rec = self.summary(uid)
            title = rec.get("title", "")
            if alt and f">{alt}" not in title and alt not in title:
                continue                                  # different allele at same position
            germ = rec.get("germline_classification", {}) or {}
            review = (germ.get("review_status") or "").strip().lower()
            return dict(
                clinvar_uid=uid,
                clinvar_accession=rec.get("accession"),
                clinvar_title=title,
                clinvar_sig=germ.get("description"),
                clinvar_review=review or None,
                clinvar_stars=REVIEW_STARS.get(review),
                clinvar_last_evaluated=germ.get("last_evaluated"),
            )
        return {}


class UcscClient(Client):
    def phylop(self, chrom: str, pos: int, track: str = "phyloP100way", genome: str = "hg38") -> float | None:
        c = chrom if chrom.startswith("chr") else f"chr{chrom}"
        body = self._request(
            f"ucsc:{track}:{c}:{pos}", "GET", f"{UCSC}/getData/track",
            params={"genome": genome, "track": track, "chrom": c, "start": pos - 1, "end": pos},
        )
        rows = body.get(track) or body.get(c) or []
        if isinstance(rows, dict):
            rows = rows.get(c, [])
        vals = [row.get("value") for row in rows if isinstance(row, dict) and row.get("value") is not None]
        return float(sum(vals) / len(vals)) if vals else None
