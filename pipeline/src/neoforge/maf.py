"""GDC masked somatic MAF ingestion (open tier).

Reads the per-case `Masked Somatic Mutation` MAF (ensemble-callset, masked by
the GDC's open-access rules) into the same column-major variant table shape
the VCF path produces, so every downstream stage stays format-agnostic.

Integrity rules enforced here:
- the file must carry exactly one tumour sample barcode (project-level masking
  exists to protect this property), with every row's normal barcode also
  singular and distinct — the same-patient guarantee upstream layers rely on;
- `Variant_Type` must be one of SNP/INS/DEL (dispatch table, unknown types
  counted and surfaced, never silently dropped);
- tumour VAF does not exist as a column in the masked MAF: it is *computed*
  as t_alt_count / t_depth — the single documented transform, recorded in the
  returned table's `tvaf_source`; rows where either field is absent or zero
  depth stay `None` (unavailable, never zero).
"""
from __future__ import annotations

import gzip
from pathlib import Path

_REQUIRED = ("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
             "Variant_Type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
             "t_depth", "t_alt_count", "n_depth")

_TABLE_COLUMNS = ["chrom", "pos", "ref", "alt", "filters", "tvaf", "nvaf", "tumor_depth"]


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def read_gdc_maf(path: str | Path) -> dict:
    path = Path(path)
    with _open(path) as handle:
        rows = [line for line in handle if not line.startswith("#")]
    if not rows:
        raise ValueError(f"{path}: no data rows after '#' comments")
    header = rows[0].rstrip("\n").split("\t")
    missing = [c for c in _REQUIRED if c not in header]
    if missing:
        raise ValueError(f"{path}: MAF header missing columns {missing}")
    idx = {c: header.index(c) for c in _REQUIRED}

    tumour_barcodes: set[str] = set()
    normal_barcodes: set[str] = set()
    table = {c: [] for c in _TABLE_COLUMNS}
    counters = {"records": 0, "snp": 0, "ins": 0, "del": 0, "other_type": 0,
                "computeable_tvaf": 0, "non_int_coordinates": 0,
                "empty_normal_barcode": 0}

    def cell(row: list[str], col: str) -> str:
        v = row[idx[col]]
        return "" if v in ("NA", "[Not Available]") else v

    for line_no, line in enumerate(rows[1:], start=2):
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < len(header):
            short = len(fields) - len(header)
            raise ValueError(f"{path}: line {line_no} short by {-short} columns")
        counters["records"] += 1
        vtype = cell(fields, "Variant_Type")
        if vtype == "SNP":
            counters["snp"] += 1
        elif vtype == "INS":
            counters["ins"] += 1
        elif vtype == "DEL":
            counters["del"] += 1
        else:
            counters["other_type"] += 1

        tumour_barcodes.add(cell(fields, "Tumor_Sample_Barcode"))
        normal = cell(fields, "Matched_Norm_Sample_Barcode")
        if normal:
            normal_barcodes.add(normal)
        else:
            counters["empty_normal_barcode"] += 1

        kind = "snp" if vtype in ("SNP", "INS", "DEL") else "other"
        if kind != "other":
            try:
                pos = int(cell(fields, "Start_Position"))
                t_depth = int(cell(fields, "t_depth")) if cell(fields, "t_depth") else None
                t_alt = int(cell(fields, "t_alt_count")) if cell(fields, "t_alt_count") else None
                if cell(fields, "n_depth"):
                    int(cell(fields, "n_depth"))
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_no} non-integer coordinate or depth ({exc})") from exc
            tvaf = t_alt / t_depth if (t_depth and t_alt is not None) else None
            if tvaf is not None:
                counters["computeable_tvaf"] += 1
            table["chrom"].append(cell(fields, "Chromosome"))
            table["pos"].append(pos)
            table["ref"].append(cell(fields, "Reference_Allele"))
            table["alt"].append(cell(fields, "Tumor_Seq_Allele2"))
            table["filters"].append(vtype)
            table["tvaf"].append(tvaf)
            table["nvaf"].append(None)
            table["tumor_depth"].append(t_depth)

    if len(tumour_barcodes) != 1:
        raise ValueError(f"{path}: expected exactly one tumour sample barcode, found {sorted(tumour_barcodes)}")
    if len(normal_barcodes) > 1:
        raise ValueError(f"{path}: expected at most one matched normal barcode, found {sorted(normal_barcodes)}")
    if normal_barcodes and normal_barcodes == tumour_barcodes:
        raise ValueError(f"{path}: tumour and normal barcodes are identical")

    return {
        "table": {**table, "columns": _TABLE_COLUMNS,
                  "null_encoding": "null = value unavailable/computed-input absent, not zero",
                  "tvaf_source": "computed: t_alt_count / t_depth (masked MAF carries no VAF column)"},
        "counters": counters,
        "tumour_sample": next(iter(tumour_barcodes), None),
        "normal_sample": next(iter(normal_barcodes), None),
    }
