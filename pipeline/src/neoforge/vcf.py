"""Dependency-free reader for the VCF subset used by somatic call-set artifacts.

Reads only fixed columns plus FILTER tokens and INFO key=dictionary entries —
the fields the L1 layer reasons about. Malformed rows raise with the line
number instead of being skipped: a truth-set that half-parses must fail loudly.
"""
from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

_INFO_BARE = object()


@dataclass(frozen=True)
class VcfMeta:
    file_format: str | None
    reference: str | None
    contigs: tuple[str, ...]
    filter_definitions: tuple[str, ...]
    n_header_lines: int


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    vid: str
    ref: str
    alts: tuple[str, ...]
    filters: tuple[str, ...]
    info: dict[str, object]

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alts) == 1 and len(self.alts[0]) == 1


class VcfFormatError(ValueError):
    pass


def _parse_info_field(field_value: str) -> tuple[str, object]:
    if "=" not in field_value:
        return field_value, _INFO_BARE
    key, value = field_value.split("=", 1)
    return key, value


class VcfReader:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(self.path)

    def _text_lines(self) -> Iterator[tuple[int, str]]:
        opener = gzip.open if self.path.suffix == ".gz" else open
        with opener(self.path, "rt", encoding="ascii") as handle:
            for lineno, line in enumerate(handle, start=1):
                yield lineno, line.rstrip("\n")

    def read_metadata(self) -> VcfMeta:
        file_format: str | None = None
        reference: str | None = None
        contigs: list[str] = []
        filters: list[str] = []
        header_lines = 0
        for lineno, line in self._text_lines():
            if not line.startswith("#"):
                break
            header_lines += 1
            if line.startswith("##fileformat="):
                file_format = line.split("=", 1)[1]
            elif line.startswith("##reference="):
                reference = line.split("=", 1)[1]
            elif line.startswith("##contig=<ID="):
                contigs.append(line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">"))
            elif line.startswith("##FILTER=<ID="):
                filters.append(line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">"))
        return VcfMeta(
            file_format=file_format,
            reference=reference,
            contigs=tuple(contigs),
            filter_definitions=tuple(filters),
            n_header_lines=header_lines,
        )

    def records(self) -> Iterator[VcfRecord]:
        for lineno, line in self._text_lines():
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                raise VcfFormatError(f"{self.path}:{lineno}: expected >=8 columns, got {len(cols)}")
            chrom, pos, vid, ref, _qual, filt, _info_unused = cols[0], cols[1], cols[2], cols[3], cols[5], cols[6], cols[7]
            info_field = cols[7]
            try:
                pos_i = int(pos)
            except ValueError as exc:
                raise VcfFormatError(f"{self.path}:{lineno}: non-integer POS {pos!r}") from exc
            for allele in (ref, *cols[4].split(",")):
                if not allele or set(allele) - set("ACGTN"):
                    raise VcfFormatError(f"{self.path}:{lineno}: non-nucleotide allele {allele!r}")
            filters_parsed = tuple(tok for tok in filt.split(";") if tok and tok != ".")
            info: dict[str, object] = {}
            if info_field not in (".", ""):
                parsed = (_parse_info_field(part) for part in info_field.split(";"))
                info = {key: value for key, value in parsed if key}
            yield VcfRecord(
                chrom=chrom,
                pos=pos_i,
                vid=vid,
                ref=ref,
                alts=tuple(cols[4].split(",")),
                filters=filters_parsed,
                info=info,
            )
