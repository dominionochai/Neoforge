"""HLA-I genotype handling with explicit resolution semantics.

Variant strings are either well-formed HLA class-I alleles or carry an
explicit failure reason. Whether a well-formed allele is *supported by the
predictor* is decided in the presentation stage, at the moment a prediction
is actually attempted — never assumed here. Nothing is silently coerced."""
from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

_ALLELE_RE = re.compile(r"^HLA-[A-C]\*\d{2}:\d{2,}$")


@dataclass(frozen=True)
class AlleleStatus:
    raw: str
    normalised: str | None
    locus: str | None
    note: str | None

    @property
    def wellformed(self) -> bool:
        return self.normalised is not None


def normalise_allele(raw: str) -> tuple[str | None, str | None]:
    text = raw.strip().upper()
    if not text.startswith("HLA-"):
        text = f"HLA-{text}"
    if _ALLELE_RE.match(text):
        return text, text.split("-")[1].split("*")[0]
    return None, f"not a well-formed HLA class-I allele string: {raw!r}"


def normalise_alleles(raw_alleles: list[str]) -> list[AlleleStatus]:
    statuses: list[AlleleStatus] = []
    seen: set[str] = set()
    for raw in raw_alleles:
        normalised, locus_or_note = normalise_allele(raw)
        if normalised is None:
            statuses.append(AlleleStatus(raw=raw, normalised=None, locus=None, note=locus_or_note))
            continue
        if normalised in seen:
            statuses.append(AlleleStatus(raw=raw, normalised=None, locus=None,
                                         note=f"duplicate allele {normalised} in genotype"))
            continue
        seen.add(normalised)
        statuses.append(AlleleStatus(raw=raw, normalised=normalised,
                                     locus=locus_or_note, note=None))
    return statuses


@dataclass(frozen=True)
class HlaSource:
    cell_line: str
    alleles: tuple[str, ...]
    source_file: Path
    sha256: str


def load_cell_line_alleles(csv_path: str | Path, cell_line: str) -> HlaSource:
    """Read one cell line's class-I genotype from the OpenVax/TCLP CSV
    (columns: cell_line, semicolon-separated HLA-I alleles)."""
    path = Path(csv_path)
    with path.open(newline="") as handle:
        rows = {row[0]: row[1] for row in csv.reader(handle) if row}
    raw = rows.get(cell_line)
    if raw is None:
        raise KeyError(f"cell line {cell_line!r} absent from {path.name}")
    alleles = tuple(a.strip() for a in raw.split(";") if a.strip())
    return HlaSource(cell_line=cell_line, alleles=alleles,
                     source_file=path, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
