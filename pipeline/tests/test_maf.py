"""GDC masked MAF reader: table shape, computed TVAF, sample-barcode guards."""
from __future__ import annotations

import gzip

import pytest

from neoforge.maf import read_gdc_maf

HEADER = ("\t".join(["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele",
                     "Tumor_Seq_Allele2", "Variant_Type", "Tumor_Sample_Barcode",
                     "Matched_Norm_Sample_Barcode", "t_depth", "t_alt_count", "n_depth",
                     "HGVSp", "Transcript_ID"]))


def row(chrom, pos, ref, alt, vtype="SNP", t_depth="27", t_alt="6",
        tumour="CT-A-06A", normal="CT-A-10A"):
    return "\t".join(["GENE", chrom, pos, ref, alt, vtype, tumour, normal,
                      t_depth, t_alt, "41", "p.Pro1Ser", "ENST00000000000"])


def write_maf(path, lines, compress=False):
    text = HEADER + "\n" + "\n".join(lines) + "\n"
    if compress:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    return path


def test_snv_rows_and_computed_tvaf(tmp_path):
    path = write_maf(tmp_path / "a.maf.gz", [row("chr1", "100", "G", "A")], compress=True)
    out = read_gdc_maf(path)
    t = out["table"]
    assert out["tumour_sample"] == "CT-A-06A" and out["normal_sample"] == "CT-A-10A"
    assert t["chrom"] == ["chr1"] and t["pos"] == [100] and t["ref"] == ["G"] and t["alt"] == ["A"]
    assert t["tvaf"] == [6 / 27] and t["tumor_depth"] == [27]
    assert out["counters"]["snp"] == 1 and out["counters"]["ins"] == 0
    assert out["counters"]["computeable_tvaf"] == 1


def test_ins_del_counted_but_separated_by_ingest_contract(tmp_path):
    path = write_maf(tmp_path / "a.maf", [
        row("chr1", "100", "G", "A"), row("chr2", "200", "C", "CT", "INS"), row("chr3", "300", "GA", "G", "DEL")])
    out = read_gdc_maf(path)
    assert out["counters"]["snp"] == 1 and out["counters"]["ins"] == 1 and out["counters"]["del"] == 1
    assert out["table"]["filters"] == ["SNP", "INS", "DEL"]


def test_multiple_tumour_barcodes_refused(tmp_path):
    path = write_maf(tmp_path / "a.maf", [row("chr1", "100", "G", "A"), row("chr2", "200", "C", "T", tumour="CT-B-06A")])
    with pytest.raises(ValueError, match="tumour sample barcode"):
        read_gdc_maf(path)


def test_identical_tumour_and_normal_refused(tmp_path):
    path = write_maf(tmp_path / "a.maf", [row("chr1", "100", "G", "A", normal="CT-A-06A")])
    with pytest.raises(ValueError, match="identical"):
        read_gdc_maf(path)


def test_zero_or_missing_depth_stays_none(tmp_path):
    path = write_maf(tmp_path / "a.maf", [row("chr1", "100", "G", "A", t_depth="0", t_alt="0"),
                                          row("chr2", "200", "C", "T", t_depth="NA", t_alt="NA")])
    out = read_gdc_maf(path)
    assert out["table"]["tvaf"] == [None, None]
    assert out["counters"]["computeable_tvaf"] == 0


def test_missing_required_column_refused(tmp_path):
    cells = ["GENE", "chr1", "100", "G", "A", "SNP", "CT-A-06A", "CT-A-10A", "27", "6", "p.Pro1Ser"]
    header = "\t".join(["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele",
                        "Tumor_Seq_Allele2", "Variant_Type", "Tumor_Sample_Barcode",
                        "Matched_Norm_Sample_Barcode", "t_depth", "t_alt_count", "HGVSp"])
    path = tmp_path / "a.maf"
    path.write_text(header + "\n" + "\t".join(cells) + "\n")
    with pytest.raises(ValueError, match="missing columns"):
        read_gdc_maf(path)
