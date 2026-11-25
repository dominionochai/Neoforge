"""Measured expression (L1): DepMap TPM for one cell-line model.

Reads the checksum-pinned per-model row extracted from
`OmicsExpressionProteinCodingGenesTPMLogp1.csv` (columns 'gene_name
(entrez_id)', values log2(TPM+1)) and the Ensembl->Entrez crosswalk from
Ensembl BioMart, and produces a gene->TPM map keyed by *unversioned* Ensembl
gene ids to match GENCODE v46 annotation ids.

Also reads GDC `rna_seq.augmented_star_gene_counts.tsv` (one patient's
STAR alignment; column `tpm_unstranded` is TPM as published by GDC — no
transform applied, the value is copied verbatim) keyed by unversioned ENSG.

Honesty rules: the single documented transform is tpm = 2^(v) - 1 (the exact
inverse of the source's stated log2(TPM+1)); an entrez id that maps to
multiple Ensembl genes is ambiguous and that column is dropped and counted —
never guessed; a missing value stays missing.
"""
from __future__ import annotations

import csv
import re

_COLUMN_ID = re.compile(r"^(?P<symbol>[^()]+)\((?P<entrez>[0-9]+)\)$")


def crosswalk(path) -> tuple[dict[str, str], int]:
    ensg_by_entrez: dict[str, str] = {}
    dropped: set[str] = set()
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 2 or not row[1]:
                continue
            ensg, entrez = row[0], row[1]
            existing = ensg_by_entrez.get(entrez)
            if existing and existing != ensg:
                dropped.add(entrez)
                del ensg_by_entrez[entrez]
            elif existing is None:
                ensg_by_entrez[entrez] = ensg
    return ensg_by_entrez, len(dropped)


def gene_tpm(row_path, ensg_by_entrez: dict[str, str]) -> tuple[dict[str, float], dict]:
    with open(row_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        _header = next(reader)  # "depmap_id, gene (entrez), ..."
        row = next(reader)
        depmap_id, values = row[0], row
    tpm_by_gene: dict[str, float] = {}
    counts = {"columns": 0, "unmapped_entrez": 0, "value_gaps": 0}
    for name, value in zip(_header[1:], values[1:]):
        counts["columns"] += 1
        m = _COLUMN_ID.match(name)
        if not m:
            counts["unmapped_entrez"] += 1
            continue
        ensg = ensg_by_entrez.get(m.group("entrez"))
        if ensg is None:
            counts["unmapped_entrez"] += 1
            continue
        if value == "":
            counts["value_gaps"] += 1
            continue
        tpm_by_gene.setdefault(ensg, 2.0 ** float(value) - 1.0)
    return depmap_id, tpm_by_gene, counts


def gene_tpm_gdc_star(tsv_path) -> tuple[dict[str, float], dict]:
    tpm_by_gene: dict[str, float] = {}
    counts = {"rows": 0, "value_gaps": 0, "duplicate_gene_id": 0}
    with open(tsv_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] == "gene_id":
                header = cols
                tpm_idx = header.index("tpm_unstranded")
                continue
            counts["rows"] += 1
            value = cols[tpm_idx]
            if value == "":
                counts["value_gaps"] += 1
                continue
            ensg = cols[0].split(".", 1)[0]
            if ensg in tpm_by_gene:
                counts["duplicate_gene_id"] += 1
                continue
            tpm_by_gene[ensg] = float(value)
    return tpm_by_gene, counts
