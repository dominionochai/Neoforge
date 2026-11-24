"""Extract per-cell-line expression rows from the DepMap public release.

Provenance chain (files kept under refs/ccle/, checksums pinned in
fetch/manifests.csv):

1. DepMap Public 23Q2 (CC BY 4.0), figshare article 22765112
   - `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (441,143,067 bytes,
     sha256 0d031e2e691e7bd7db195b28a9af92e1f7fefefb27a0e7dd11c56fb9a7d36ca3):
     matrix values are log2(TPM+1) per the file's own column header; rows are
     DepMap model IDs; columns are 'gene_name (entrez_id)';
   - `Model.csv` verifies the ID/cell-line mapping (ACH-000699 == HCC1395);
2. Ensembl BioMart (martservice, hsapiens_gene_ensembl, release-aligned
   ensembl_gene_id -> entrezgene_id): `biomart_ensembl_entrez.tsv` — the
   Ensembl->Entrez crosswalk the pipeline needs to join the ENSG gene ids on
   the GENCODE v46 annotation to DepMap columns (36,824 pairs with Entrez).

This script verifies the pinned checksums, then outputs one verbatim row per
requested model. Values stay on the source's log2(TPM+1) scale; the pipeline
applies the single documented transform tpm = 2^(v) - 1 and never imputes
missing genes.
"""
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

REFS = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/workspace/refs/ccle")
# ACH-000699 == HCC1395 per Model.csv. The matched normal HCC1395BL
# (ACH-002382) is listed in Model.csv but absent from the 23Q2 RNA matrix —
# noted as unavailable rather than substituted.
MODELS = {
    "ACH-000699": "hcc1395_expression_tpm_log2p1.csv",
}
PINNED = {
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv":
        "0d031e2e691e7bd7db195b28a9af92e1f7fefefb27a0e7dd11c56fb9a7d36ca3",
    "Model.csv":
        "aef4502b249f88a6e26cd02bbd1c04b6905365b542225ccbbbf09f650a92d4d4",
    "biomart_ensembl_entrez.tsv":
        "1bd390e19458f0223ef019285410ad5639805182802474ac63efc9fff331a4d8",
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    for name, digest in PINNED.items():
        path = REFS / name
        if not path.exists():
            raise SystemExit(f"missing source {path}")
        if sha256_of(path) != digest:
            raise SystemExit(f"checksum mismatch for {path}")

    found: set[str] = set()
    with (REFS / "OmicsExpressionProteinCodingGenesTPMLogp1.csv").open(
            newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if row and row[0] in MODELS:
                out = REFS / MODELS[row[0]]
                with out.open("w", newline="", encoding="utf-8") as g:
                    writer = csv.writer(g)
                    writer.writerow(["depmap_id"] + header[1:])
                    writer.writerow([row[0]] + row[1:])
                print(f"extracted {row[0]} -> {out.name}")
                found.add(row[0])
            if found == set(MODELS):
                break
    missing = set(MODELS) - found
    if missing:
        raise SystemExit(f"models not found in expression matrix: {sorted(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
