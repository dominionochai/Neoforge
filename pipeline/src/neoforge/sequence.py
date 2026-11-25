"""Random-access sequence provider over the reference genome.

The translate stage depends only on a small `fetch(chrom, start0, end0)`
protocol, so production uses a pysam-indexed FASTA and tests can inject an
in-memory genome. Plain-gzip GENCODE assemblies are decompressed once so a
random-access index can be built."""
from __future__ import annotations

import gzip
import shutil
from pathlib import Path


class GenomeHandle:
    """Protocol: fetch(chrom, start0, end0) -> str, pysam-compatible."""

    def fetch(self, chrom: str, start0: int, end0: int) -> str:  # pragma: no cover - protocol
        raise NotImplementedError


class IndexedFasta(GenomeHandle):
    def __init__(self, fa_path: str | Path) -> None:
        import pysam  # noqa: PLC0415 — keeps pysam optional outside translate

        self._fa = pysam.FastaFile(str(fa_path))
        self._names = frozenset(self._fa.references)

    @classmethod
    def decompress_if_needed(cls, gz_path: str | Path, out_dir: Path) -> Path:
        gz_path, out_dir = Path(gz_path), Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        name = gz_path.name[:-3] if gz_path.name.endswith(".gz") else gz_path.name
        fa = out_dir / name
        if fa.exists() and fa.stat().st_size > 0:
            return fa
        tmp = fa.with_name(name + ".tmp")
        with gzip.open(gz_path, "rb") as src, tmp.open("wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 22)
        tmp.rename(fa)
        return fa

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        if chrom not in self._names:
            raise KeyError(f"chromosome {chrom!r} absent from reference FASTA")
        return self._fa.fetch(chrom, start0, end0).upper()

    def close(self) -> None:
        self._fa.close()


class InMemoryGenome(GenomeHandle):
    """Test double: chromosome name -> full sequence string."""

    def __init__(self, chromosomes: dict[str, str]) -> None:
        self._chroms = {k: v.upper() for k, v in chromosomes.items()}

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        return self._chroms[chrom][start0:end0]
