"""GENCODE GTF transcript index.

Parses CDS blocks of protein-coding transcripts once and answers two questions
fast for tens of thousands of variants:

1. which transcripts have a CDS exon covering genomic position p;
2. for a transcript, the CDS-relative (translation-direction) offset of a
   genomic position, so codons can be assembled exactly across exon junctions
   on either strand — never derived from string manipulation of HGVSp.
"""
from __future__ import annotations

import gzip
from pathlib import Path

BIN_SHIFT = 16
BIN_SIZE = 1 << BIN_SHIFT


class GtfFormatError(ValueError):
    pass


class TranscriptCds:
    __slots__ = ("transcript_id", "chrom", "strand", "exons", "first_frame", "_cums")

    def __init__(self, transcript_id: str, chrom: str, strand: str,
                 exons: tuple[tuple[int, int], ...], first_frame: int) -> None:
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.exons = exons
        self.first_frame = first_frame
        cums = [0]
        for start, end in exons:
            cums.append(cums[-1] + (end - start + 1))
        self._cums = tuple(cums)

    @property
    def n_codons(self) -> int:
        return self._cums[-1] // 3

    def offset_of(self, pos: int) -> int | None:
        """Translation-direction offset (0-based) of 1-based genomic `pos`."""
        for i, (start, end) in enumerate(self.exons):
            if start <= pos <= end:
                local = pos - start + 1 if self.strand == "+" else end - pos + 1
                return self._cums[i] + local - 1
        return None

    def positions_for_offsets(self, offsets: list[int]) -> list[int]:
        """Genomic 1-based positions of translation-direction offsets. A codon
        spans at most two blocks and callers pass short ascending runs, so a
        linear forward walk over blocks is the right granularity."""
        out: list[int] = []
        idx = 0
        for o in offsets:
            while self._cums[idx + 1] <= o:
                idx += 1
            within = o - self._cums[idx]
            start, end = self.exons[idx]
            out.append(start + within if self.strand == "+" else end - within)
        return out


class GtfTranscriptIndex:
    def __init__(self, transcripts: dict[str, TranscriptCds]) -> None:
        self.transcripts = transcripts
        bins: dict[tuple[str, int], list[tuple[int, int, str]]] = {}
        for tid, t in transcripts.items():
            for start, end in t.exons:
                for b in range(start >> BIN_SHIFT, (end >> BIN_SHIFT) + 1):
                    bins.setdefault((t.chrom, b), []).append((start, end, tid))
        self._bins = bins

    @classmethod
    def from_gtf(cls, path: str | Path) -> "GtfTranscriptIndex":
        cds_blocks: dict[str, list[tuple[int, int, int]]] = {}
        meta: dict[str, tuple[str, str]] = {}
        opener = gzip.open if Path(path).suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for lineno, line in enumerate(handle, start=1):
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] != "CDS":
                    continue
                attrs = _parse_attrs(cols[8])
                tid = attrs.get("transcript_id")
                if not tid:
                    continue
                if attrs.get("transcript_type", "protein_coding") != "protein_coding":
                    continue
                try:
                    start, end = int(cols[3]), int(cols[4])
                except ValueError as exc:
                    raise GtfFormatError(f"{path}:{lineno}: non-integer CDS coordinates") from exc
                if end < start:
                    raise GtfFormatError(f"{path}:{lineno}: inverted block")
                frame = 0 if cols[7] == "." else int(cols[7])
                cds_blocks.setdefault(tid, []).append((start, end, frame))
                meta[tid] = (cols[0], cols[6])

        transcripts: dict[str, TranscriptCds] = {}
        for tid, blocks in cds_blocks.items():
            chrom, strand = meta[tid]
            ordered = tuple(sorted(blocks, reverse=(strand == "-")))
            coordinates = tuple((s, e) for s, e, _f in ordered)
            first_frame = ordered[0][2]
            transcripts[tid] = TranscriptCds(tid, chrom, strand, coordinates, first_frame)
        return cls(transcripts)

    def covering_transcripts(self, chrom: str, pos: int) -> list[str]:
        entries = self._bins.get((chrom, pos >> BIN_SHIFT))
        if not entries:
            return []
        return sorted({tid for start, end, tid in entries if start <= pos <= end})

    def n_transcripts(self) -> int:
        return len(self.transcripts)


def _parse_attrs(field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in field.strip().rstrip(";").split(";"):
        part = part.strip()
        if " " not in part:
            continue
        key, value = part.split(" ", 1)
        value = value.strip('"')
        existing = attrs.get(key)
        attrs[key] = existing + "," + value if existing and existing != value else value
    return attrs
