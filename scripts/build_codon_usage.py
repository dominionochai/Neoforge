"""Build the pinned human codon-usage table from the pinned GENCODE v46 CDS.

No external codon-usage download: the reference transcripts the translate stage
already uses are the source, so the codon-choice statistics and the pipeline
share one documented, checksum-pinned input. Counts (not fractions) are emitted
so the table is exact; determinism follows from the pinned genome/GTF.

Usage:
  python scripts/build_codon_usage.py --refs-dir /workspace/refs \
      --out pipeline/src/neoforge/data/codon_usage_human_GRCh38_GENCODE_v46_CDS.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline" / "src"))

from neoforge.gtf import GtfTranscriptIndex          # noqa: E402
from neoforge.sequence import IndexedFasta           # noqa: E402
from neoforge.stages.translate import CODON_TABLE    # noqa: E402

OUT_NAME = "codon_usage_human_GRCh38_GENCODE_v46_CDS.tsv"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refs-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    gtf = args.refs_dir / "gencode/gencode.v46.annotation.gtf.gz"
    fasta = args.refs_dir / "gencode/GRCh38.primary_assembly.genome.fa.gz"
    gtf_sha, fasta_sha = sha256_of(gtf), sha256_of(fasta)

    fasta_plain = IndexedFasta.decompress_if_needed(fasta, args.refs_dir / "gencode")
    index = GtfTranscriptIndex.from_gtf(gtf)
    genome = IndexedFasta(fasta_plain)

    counts: dict[str, int] = {}
    stats = {"transcripts_used": 0, "incomplete_5prime_excluded": 0,
             "length_not_multiple_of_3_excluded": 0, "codons_with_ambiguous_base_skipped": 0,
             "stop_codons_excluded": 0}
    for tid, t in sorted(index.transcripts.items()):
        if t.first_frame != 0:
            stats["incomplete_5prime_excluded"] += 1
            continue
        parts = []
        for start, end in t.exons:
            part = genome.fetch(t.chrom, start - 1, end)
            if t.strand == "+":
                parts.append(part)
            else:
                parts.append(part.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1])
        cds = "".join(parts).upper()
        whole, tail = divmod(len(cds), 3)
        if tail:
            stats["length_not_multiple_of_3_excluded"] += 1
        for i in range(whole - 1):
            codon = cds[i * 3:i * 3 + 3]
            if any(b not in "ACGT" for b in codon):
                stats["codons_with_ambiguous_base_skipped"] += 1
                continue
            counts[codon] = counts.get(codon, 0) + 1
        last = cds[(whole - 1) * 3:whole * 3] if whole else ""
        if CODON_TABLE.get(last) == "*":
            stats["stop_codons_excluded"] += 1
        else:
            if last and not any(b not in "ACGT" for b in last):
                counts[last] = counts.get(last, 0) + 1
        stats["transcripts_used"] += 1

    lines = [
        f"# {OUT_NAME[:-4]}",
        "# Human codon usage: raw CDS codon counts, derived from the pinned GENCODE v46",
        "# GRCh38 protein-coding transcript set (the same reference the translate stage",
        "# uses for peptide derivation), so codon choice statistics and the pipeline",
        "# share one checksum-pinned source.",
        "# Method (scripts/build_codon_usage.py, deterministic): for every",
        "# transcript_type=protein_coding CDS in the GTF, the CDS was spliced from the",
        "# GRCh38 primary-assembly genome strand-aware in translation direction, in",
        "# frame from offset 0; transcripts with a non-zero first CDS frame (5'-",
        "# incomplete) were excluded entirely; the final stop codon was excluded;",
        "# codons containing non-ACGT bases were skipped and counted.",
        f"# input gencode.v46.annotation.gtf.gz sha256={gtf_sha}",
        f"# input GRCh38.primary_assembly.genome.fa.gz sha256={fasta_sha}",
        "# licence: GENCODE is distributed under Apache 2.0 (gencodegenes.org terms).",
        "# codon\taa\tcount",
    ]
    for codon in sorted(counts):
        lines.append(f"{codon}\t{CODON_TABLE[codon]}\t{counts[codon]}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    total = sum(counts.values())
    print(f"{OUT_NAME[:-4]}: {total} codons over {stats['transcripts_used']} transcripts; {len(counts)} of 64 codons observed")
    for k, v in stats.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
