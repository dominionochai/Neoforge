# Data Provenance

All sandbox verifications below were executed and recorded during the research
phase. Reference file copies live outside the repo (downloadable via
`fetch/fetch_manifests.py`); their SHA-256 checksums are pinned here.

## 1. Primary demonstration case — HCC1395 / HCC1395BL (SEQC2)

**Same-biological-pair guarantee.** This case is a tumour cell line (HCC1395,
triple-negative breast cancer) and its matched normal (HCC1395BL, EBV B-LCL)
from the same donor, characterised by the SEQC2 consortium as the reference
tumour–normal pair. Every input traces to that pair; no cross-patient merging.

| Input | Source | URL | SHA-256 | Licence |
|---|---|---|---|---|
| SNV truth set (39,447 SNVs, v1.2.1) | NCBI ReferenceSamples / SEQC2 Somatic Mutation WG | `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` | `0c49f7e1519fc7860b773f79bdcf07453fe4def8303707f34542471d92494c9d` | US Govt public data |
| INDEL truth set (v1.2.1) | same directory | `.../high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz` | `db7dc8be2ec04bab830ff1411146c7d538c217ff684e874b85d98a3dbbda91cf` | US Govt public data |
| High-confidence regions BED | same directory | `.../High-Confidence_Regions_v1.2.bed` | `dedaa36a956d6121fb8fccd9ae1891a4fefd608ff31ab277cb55d7c96bb4a6c4` | US Govt public data |
| HLA-I genotype (5 alleles; HLA-A LOH documented) | OpenVax/TCLP `cell-line-class1-mhc-alleles.csv` (HCC1395 row) | `https://github.com/openvax/cancer-cell-line-mhc-alleles/blob/master/cell-line-class1-mhc-alleles.csv` | `940e1d1a56be958c470332ef5ce4ba746241c062c95047e1f30a85e12c739860` | public repository (cite TCLP) |
| Expression (TPM) | CCLE / DepMap, HCC1395 | DepMap portal (CC BY 4.0) | pinned at fetch time | CC BY 4.0 |
| Reference genome + gene model | GENCODE v46 | `gencode.v46.annotation.gtf.gz` (51,006,804 bytes; sha256 `82b0e52bc4956486528152e3850660b3c41767ab103d87949e41d38f92ef8b38`, MD5 verified vs GENCODE MD5SUMS) + `GRCh38.primary_assembly.genome.fa.gz` (845,635,028 bytes; sha256 `ada8bcedcd6ed9359f56778f2430375eddd5580cc4f35b6e581b53422f7035a8`, MD5 verified) + `gencode.v46.pc_translations.fa.gz` (sha256 `80a066d4d278f36dd3971e2c0f3fb9d1c74f8ae8fd5b59321ef8d66060bc039e`, MD5 verified) | verified at fetch time | GENCODE free-use |

**VAF availability (corrected after fetching the files):** the truth VCFs
do carry measured tumour allele fraction in the `TVAF=` INFO field (plus
`NVAF` and `bwaDP` depths). An earlier version of this document claimed "no
tumour VAF"; that was wrong and the pipeline now reads measured `TVAF` directly
(reported exactly as in the source; `None` stays *unavailable*, never imputed).

## 2. Context view — TCGA-SKCM (open tier, single real patients, never merged)

Per-case open harmonised somatic MAFs: chromosome, ref/alt, `t_depth`,
`t_alt_count`, `n_depth`, HGVSp, transcript IDs — verified on file
`02de9cd7-94e2-470c-b2ac-9fe6285b50db` (case `4f7e7a0e-4aa2-4b6a-936f-6429c7d6a2a7`,
submitter `TCGA-D9-A3Z1`, 599 mutations). Open STAR RNA-seq gene counts exist
per case. Context use is limited to cohort-level landscape visualisation.

**Closed path, documented:** patient-level HLA-I genotypes for TCGA (OptiType
calls from the Immune Landscape of Cancer, GDC UUID `2cf05dd5-...`) return
HTTP 403 without a GDC data-access token; patient-level personal neoantigen
analysis for TCGA cases is therefore **not presented**, rather than combining
HLA from one source with mutations from another; the ingest harness enforces
this mechanically (a bundle without `hla_csv` must declare
`hla.status="unavailable"` with its reason, and the present stage refuses to
run on an empty HLA block). GDC token access would unlock it (documented
path), but honesty first.


### 2a. Executed — SKCM case `skcm_d9a3z1` fetch records (2026 session)

| GDC file UUID | File | Content | SHA-256 | Licence/access |
|---|---|---|---|---|
| `02de9cd7-94e2-470c-b2ac-9fe6285b50db` | `tcga/d9a3z1.wxs.aliquot_ensemble_masked.maf.gz` | open Masked Somatic Mutation (WXS ensemble, 591 records; single tumour barcode `TCGA-D9-A3Z1-06A-11D-A23B-08`, matched normal `TCGA-D9-A3Z1-10A-01D-A23B-08`) | `bc37ed2ab8c2ccd1e4f905b75cbc9ce48b8a98eaee8366c8f7c435a912818b01` | GDC open tier |
| `631ef22b-0162-47a7-9dd0-c1d33f9a1de5` | `tcga/d9a3z1.rna_seq.augmented_star_gene_counts.tsv` | same patient's tumour RNA aliquot (-06A-01R): STAR gene counts incl. `tpm_unstranded` (60,666 rows; GENCODE v36 model) | `a1ca53c4e2b4907627daf5c9192aabdd3f5d835dca4b605e1d8c51c5a8686ab1` | GDC open tier |

Both are fetched verbatim via `fetch/manifests.csv`; the case bundle pins their
SHA-256 and the ingest harness re-verifies them at run time (same-patient
barcode check refuses files carrying a second tumour sample).

---

## 3. TESLA external benchmark (ranking-rule audit executed see EVALUATION_PLAN.md; candidate-level immunogenicity: not claimed)

| File | Content | SHA-256 |
|---|---|---|
| `mmc1.xlsx` | Data dictionary for all tables | `e4b230a0b09b908dc87714ba509beacff15d7e1c77c5699f51b1d1e838baab28` |
| `mmc4.xlsx` | **608 candidate pMHC, 6 melanoma patients, 37 experimentally validated immunogenic** (`VALIDATED` column; tetramer/T-cell assays; per-pMHC features) | `5cd384753464e981712af15dffa3b8cb5a3c2c5026e8cd80894a17f9df6236fe` |
| `mmc5.xlsx` | Per-team evaluation (AUPRC/TTIF baselines) | `0e502d9d2fc8f0f366bb8ad1612e9685de739a8169d13c5de7ee95f937dd0aba` |
| `mmc7.xlsx` | Independent validation cohort: 310 peptides (TIL), 4 validated | `7a1475a88449457473993de53e347f219a3c2df194b8de4580278cf516a09b6a` |

Base URL: `https://ars.els-cdn.com/content/image/1-s2.0-S0092867420311569-mmc{n}.xlsx`
(Wells DK et al., *Cell* 183:1556–1569.e17, 2020; academic benchmark; Synapse
project `syn21048999` also openly lists raw WES/RNA FASTQs for patients).
Licensing: supplementary materials of an open Cell resource article, made
publicly available by the authors as an open benchmark; used with citation;
no redistribution of the binaries inside this repository — pinned checksums +
fetch script instead.

## 4. Software provenance

- MHCflurry 2.4.0 (models release 2.3.0, `Class1PresentationPredictor`,
  `models_class1_presentation.20260928`), used unmodified; licence BSD-3.
- pysam 0.24.1, Python 3.12.14 recorded in manifests.
- NetMHCpan: **excluded** (registration/academic licence); stated as limitation.

## 5. Failure conduct

Any missing input (VAF, expression, an unresolvable HLA allele, an unavailable
source) is surfaced in artifacts and UI as `missing` with a reason — never
imputed, never silently dropped.
