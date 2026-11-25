"""Sorted interval index over a BED file (0-based half-open), with containment
and set-of-positions queries. Used to independently re-verify SEQC2 claims —
e.g. that every v1.2.1 truth call lies inside High-Confidence_Regions_v1.2.bed —
rather than trusting the release notes."""
from __future__ import annotations

import gzip
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path


class BedFormatError(ValueError):
    pass


@dataclass(frozen=True)
class IntervalIndex:
    """chrom -> (starts, ends), both ascending with starts[i] < ends[i]."""

    intervals: dict[str, tuple[list[int], list[int]]]

    @classmethod
    def from_bed(cls, path: str | Path) -> "IntervalIndex":
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        opener = gzip.open if Path(path).suffix == ".gz" else open
        with opener(path, "rt", encoding="ascii") as handle:
            for lineno, line in enumerate(handle, start=1):
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                cols = line.split()
                if len(cols) < 3:
                    raise BedFormatError(f"{path}:{lineno}: expected >=3 BED columns")
                chrom = cols[0]
                try:
                    start, end = int(cols[1]), int(cols[2])
                except ValueError as exc:
                    raise BedFormatError(f"{path}:{lineno}: non-integer interval") from exc
                if end <= start:
                    raise BedFormatError(f"{path}:{lineno}: empty interval {start}-{end}")
                by_chrom.setdefault(chrom, []).append((start, end))
        merged: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, spans in by_chrom.items():
            spans.sort()
            starts: list[int] = []
            ends: list[int] = []
            for s, e in spans:
                if starts and s <= ends[-1]:
                    ends[-1] = max(ends[-1], e)
                else:
                    starts.append(s)
                    ends.append(e)
            merged[chrom] = (starts, ends)
        return cls(intervals=merged)

    def contains(self, chrom: str, start: int, end: int | None = None) -> bool:
        """True iff 1-based inclusive bases [start, end] (or just `start`)
        fall inside one single interval."""
        spans = self.intervals.get(chrom)
        if spans is None:
            return False
        starts, ends = spans
        probe = start - 1  # BED stays 0-based
        idx = bisect_right(starts, probe) - 1
        if idx < 0 or probe >= ends[idx]:
            return False
        if end is None:
            return True
        return ends[idx] >= end

    @property
    def n_intervals(self) -> int:
        return sum(len(ends) for _starts, ends in self.intervals.values())
