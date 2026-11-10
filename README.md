<div align="center">

# MITOCLIN

### Automated Mitochondrial Genome Analysis and NUMT-Purified Variant Calling Framework

**Laboratory of Genome Informatics (LGI)**  
**Centre for DNA Fingerprinting and Diagnostics (BRIC-CDFD), Hyderabad, India**

[![Release](https://img.shields.io/badge/Release-v2.1.0-blue.svg?style=flat-square)](https://github.com/lgi-cdfd/mitoclin/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-teal.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.9%2B-informational.svg?style=flat-square)](https://www.python.org/)
[![Conda](https://img.shields.io/badge/Environment-Conda%20%7C%20Bioconda-brightgreen.svg?style=flat-square)](https://bioconda.github.io/)
[![Reference](https://img.shields.io/badge/Reference-rCRS%20(NC__012920.1)-purple.svg?style=flat-square)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)
[![NUMT Decoys](https://img.shields.io/badge/NUMT%20Decoys-766%20Contigs-orange.svg?style=flat-square)](https://genome.ucsc.edu/)
[![PhyloTree](https://img.shields.io/badge/PhyloTree-Build%2017-forestgreen.svg?style=flat-square)](https://www.phylotree.org/)
[![Institution](https://img.shields.io/badge/Institution-BRIC--CDFD%20%7C%20LGI-crimson.svg?style=flat-square)](https://cdfd.org.in/)

</div>

---

## 1. Overview

MITOCLIN is an open-source, reproducible command-line workflow and computational toolkit engineered for high-throughput mitochondrial genomics, heteroplasmy quantification, and clinical variant prioritization from paired-end high-throughput sequencing (NGS) data.

Accurate identification of mitochondrial DNA (mtDNA) variants is routinely confounded by Nuclear Mitochondrial DNA Segments (NUMTs)—ancient pseudogene fragments transposed from the ancestral mitochondrial genome into the human nuclear autosomes. When standard alignment protocols map sequencing reads exclusively against the revised Cambridge Reference Sequence (rCRS, NC_012920.1), nuclear pseudogene fragments align non-specifically to mtDNA, generating rampant false-positive low-frequency heteroplasmic variants (typically clustered between 5% and 20% Variant Allele Fraction).

MITOCLIN addresses this fundamental challenge through a **competitive dual-alignment strategy**: sequencing reads are simultaneously aligned against a composite reference comprising the canonical human mitochondrial genome (NC_012920.1) and 766 nuclear NUMT decoy contigs (RHNumtS.2 / UCSC human NUMT catalog). Reads of nuclear origin preferentially anchor to their cognate autosomal contigs and are depleted away from the mitochondrial alignment. MITOCLIN subsequently performs non-diploid somatic/mitochondrial variant calling via GATK Mutect2, evaluates phylogenetic haplogroups and cross-sample contamination thresholds using HaploCheck (PhyloTree Build 17), and integrates offline clinical pathogenicity annotations against MITOMAP and gnomAD v3.1 nuclear pseudogene (`npg`) reference datasets.

---

## 2. Pipeline Architecture & Workflow

```
Raw Paired-End FASTQs (R1 & R2)
               │
               ▼
   [Step 0: Read Pair Synchronization] (repair_fastq_pairs.py)
               │
   ┌───────────┴───────────────────────┐
   ▼                                   ▼
[Step 1: FastQC Quality Profiling]   [Step 2: Adapter & Quality Trimming]
                                                (Trim Galore / Cutadapt)
                                       │
                                       ▼
                     [Step 3: Competitive Alignment (BWA-MEM)]
                     rCRS (NC_012920.1) + 766 Nuclear NUMT Decoys
                                       │
                                       ▼
                     [Step 4: NUMT Depletion & Coordinate Normalization]
                     (numt_engine.py -> Extract authentic mtDNA BAM)
                                       │
                                       ▼
                     [Step 5: Duplicate Marking (GATK MarkDuplicates)]
                                       │
                     ┌─────────────────┴─────────────────┐
                     ▼                                   ▼
          [Step 6: BAM QC & Stats]            [Step 7: Per-Base Coverage]
          (samtools flagstat/idxstats)        (samtools depth, mean coverage)
                     │                                   │
                     └─────────────────┬─────────────────┘
                                       │
                                       ▼
                     [Step 8: Variant Discovery (GATK Mutect2)]
                     (--mitochondria-mode, native PairHMM)
                                       │
                                       ▼
                     [Step 9-14: Variant Filtering & Allelic Splitting]
                     (GATK FilterMutectCalls, BCFtools norm, VAF >= 1%)
                                       │
                     ┌─────────────────┴─────────────────┐
                     ▼                                   ▼
          [Step 15: Haplogroup & Contam]     [Step 16-17: Pathogenicity Curation]
          (HaploCheck, PhyloTree B17)        (MITOMAP local offline DB)
                     │                                   │
                     └─────────────────┬─────────────────┘
                                       │
                                       ▼
                     [Step 18: Variant NUMT & gnomAD Annotation]
                     (gnomAD v3.1 npg filter + UCSC NUMT homology)
                                       │
                                       ▼
                     Final Standardized Clinical TSV & VCF
```

### Detailed Execution Steps

1. **Read Pair Synchronization:** Validates that paired FASTQ files contain identical, synchronized read identifiers, repairing interleaved discrepancies before execution.
2. **Quality Control & Adapter Trimming:** Computes base quality profiles (FastQC) and strips adapter sequences and low-quality 3' ends with a Phred score cutoff of Q30 (Trim Galore / Cutadapt).
3. **Competitive Decoy Alignment:** Aligns trimmed reads using BWA-MEM against the composite `rCRS_plus_NUMTs` reference. Reads exhibiting higher homology to nuclear paralogs align to the decoy contigs, preventing artifactual misattribution.
4. **NUMT Depletion & Filtering:** Extracts authentic NC_012920.1 mitochondrial alignments, sorts coordinates, indexes BAMs, and exports granular read absorption statistics.
5. **Duplicate Marking:** Identifies and tags optical and PCR duplicate read pairs using GATK MarkDuplicates.
6. **Coverage Depth Profiling:** Generates comprehensive per-nucleotide depth profiles across all 16,569 bp, calculating overall mean coverage depth and regional uniformity metrics.
7. **Mitochondrial Variant Discovery:** Calls non-diploid SNPs and INDELs using GATK Mutect2 in `--mitochondria-mode` with strand-bias modeling.
8. **Variant Quality & Heteroplasmy Filtering:** Executes GATK FilterMutectCalls, splits multiallelic sites using BCFtools, and derives exact per-allele Variant Allele Fractions (VAF = Alternate Depth / Total Depth).
9. **Haplogroup & Contamination Assessment:** Leverages HaploCheck against PhyloTree Build 17 using Lowest Common Ancestor (LCA) distance, estimating major and minor haplogroup lineages and detecting sample-mix contamination down to 1.8% minor fraction.
10. **Clinical Pathogenicity Annotation:** Annotates variants against a curated local MITOMAP database (coding control mutations, rRNA/tRNA alterations, and confirmed pathogenic mutations) without requiring external web requests.
11. **NUMT Knowledgebase Cross-Referencing:** Tags variant calls against Broad Institute gnomAD v3.1 `npg` (nuclear pseudogene) flag sites and coordinates intersecting known human nuclear NUMT intervals.

---

## 3. System Requirements & Dependencies

### Hardware Requirements
* **Operating System:** Linux (CentOS 7+, RHEL 8+, Ubuntu 18.04+, Debian 10+) or macOS.
* **Processor:** 64-bit x86_64 architecture, minimum 4 CPU cores (8+ cores recommended for high throughput).
* **Memory (RAM):** 8 GB minimum, 16 GB recommended for high-depth coverage libraries.
* **Storage:** 10 GB free disk space for references, databases, and temporary alignment files.

### Bioinformatic Tool Dependencies
MITOCLIN requires the following standard command-line tools in your environment:
* **BWA** (`>= 0.7.17`)
* **SAMtools** (`>= 1.15`)
* **BCFtools** (`>= 1.15`)
* **GATK4** (`>= 4.2.0.0`)
* **FastQC** (`>= 0.11.9`)
* **Trim Galore** (`>= 0.6.7`) / **Cutadapt** (`>= 3.4`)
* **Java Runtime Environment (JRE)** (`>= 17`)
* **Python** (`>= 3.9`) with packages: `numpy`, `pandas`, `matplotlib`, `reportlab`

---

## 4. Installation

### Method A: Conda / Mamba Environment (Recommended)

The most reliable and reproducible installation method is via Conda:

```bash
# 1. Clone the MITOCLIN repository
git clone https://github.com/lgi-cdfd/mitoclin.git
cd mitoclin

# 2. Create and activate the conda environment
conda env create -f environment.yml
conda activate mitoclin

# 3. Verify environment configuration
./mitoclin --check-env
```

### Method B: Manual Installation in an Existing Environment

If your system or cluster module already contains BWA, SAMtools, GATK4, and Python:

```bash
# Clone the repository
git clone https://github.com/lgi-cdfd/mitoclin.git
cd mitoclin

# Install required Python dependencies
pip install -r requirements.txt

# Run the environment pre-flight validator
./mitoclin --check-env
```

---

## 5. Quick Start Tutorial

To verify that MITOCLIN is operating correctly on your system, execute the bundled test script:

```bash
cd mitoclin
bash example/run_example.sh
```

This runs the complete 18-step pipeline on a bundled paired-end test FASTQ dataset in under one minute and produces an output directory (`example_output/`) with indexed BAMs, filtered VCFs, coverage metrics, and the finalized clinical TSV report.

---

## 6. Command-Line Usage

```bash
./mitoclin [options] -s <SAMPLE_ID> -1 <R1.fastq.gz> -2 <R2.fastq.gz>
```

### Argument Reference

| Option | Flag | Description | Default |
| :--- | :--- | :--- | :--- |
| `--sample` | `-s` | Unique sample identifier (alphanumeric, dashes, underscores) | Required |
| `--r1` | `-1` | Path to forward paired-end FASTQ (`.fastq`, `.fq.gz`) | Required |
| `--r2` | `-2` | Path to reverse paired-end FASTQ (`.fastq`, `.fq.gz`) | Required |
| `--outdir` | `-o` | Output directory destination | `./output_<SAMPLE_ID>` |
| `--reference` | `-r` | Reference genome FASTA (`numt`, `rcrs`, or custom path) | `rCRS_plus_NUMTs.fasta` |
| `--threads` | `-t` | Number of CPU cores for multi-threaded alignment & calling | `4` |
| `--min-vaf` | | Minimum Variant Allele Fraction cutoff | `0.01` (1% VAF) |
| `--check-env` | | Perform pre-flight dependency check and exit | — |
| `--help` | `-h` | Display usage instructions and exit | — |

### Execution Examples

#### 1. Standard Clinical Run (NUMT-Aware Competitive Alignment)
```bash
./mitoclin \
    --sample Patient01 \
    --r1 data/Patient01_R1.fastq.gz \
    --r2 data/Patient01_R2.fastq.gz \
    --outdir results/Patient01 \
    --threads 8
```

#### 2. Fast Screening Mode (rCRS-Only Alignment)
To run in standard rCRS-only mode without nuclear decoy contigs:
```bash
./mitoclin \
    --sample Patient01 \
    --r1 data/Patient01_R1.fastq.gz \
    --r2 data/Patient01_R2.fastq.gz \
    --reference rcrs \
    --threads 4
```

#### 3. Positional Shorthand Syntax
```bash
./mitoclin run Patient01 data/Patient01_R1.fq.gz data/Patient01_R2.fq.gz results/Patient01
```

---

## 7. Directory Structure and Output Files

Following execution, the output directory is structured as follows:

```
output_<SAMPLE_ID>/
├── 01_raw_fastqc/              # Raw sequence quality reports (HTML, ZIP)
├── 02_trimmed/                 # Adapter-trimmed and synchronized FASTQ reads
├── 03_mapping/
│   ├── <SAMPLE>_sorted.bam     # Coordinated-sorted mitochondrial BAM
│   ├── <SAMPLE>_sorted.bam.bai # BAM index
│   ├── <SAMPLE>_dedup.bam      # Deduplicated mitochondrial BAM
│   ├── <SAMPLE>_dedup.bam.bai  # Deduplicated BAM index
│   └── <SAMPLE>_numt_depletion.json # Quantitative NUMT read depletion metrics
├── 04_variant_calling/
│   ├── <SAMPLE>_unfiltered.vcf # Raw candidate variant calls (Mutect2)
│   ├── <SAMPLE>_filtered.vcf   # Filtered variants passing quality criteria
│   └── <SAMPLE>_snps.vcf       # Cleaned, biallelic mitochondrial SNPs
├── 05_bam_stats/               # Alignment statistics (flagstat, idxstats, stats)
├── 06_coverage/
│   ├── <SAMPLE>_coverage.txt   # Per-base read depth across 16,569 bp
│   └── <SAMPLE>_mean_coverage.txt # Mean cohort coverage summary
├── 07_haplocheck/
│   └── <SAMPLE>_haplocheck.tsv # Major/minor haplogroups and contamination %
├── 09_mitomaster/              # Parsed mutational functional annotations
├── 10_final_report_data/
│   └── <SAMPLE>_final_report.tsv # Primary finalized clinical report table
└── logs/                       # Comprehensive step-by-step stdout/stderr runtime logs
```

### Primary Output Table Columns (`_final_report.tsv`)

The primary deliverable is a tab-delimited file containing prioritized, clinically curated variant information:

| Column Header | Description |
| :--- | :--- |
| `Sample` | Sample identifier |
| `POS` | Genomic coordinate on rCRS (1..16,569) |
| `REF` | Reference nucleotide allele |
| `ALT` | Alternate nucleotide allele |
| `Mutation` | Standard mitochondrial HGVS nomenclature (e.g., `m.11778G>A`) |
| `Gene` | Mitochondrial gene symbol (*MT-ND4*, *MT-ATP6*, *MT-TL1*, etc.) |
| `Protein_Change` | Amino acid alteration (e.g., `p.R340H`) |
| `Disease` | Curated disease association (LHON, Leigh, MELAS, MERRF) |
| `Pathogenicity` | Clinical classification (Confirmed Pathogenic, Likely Pathogenic, Polymorphism) |
| `Coverage` | Read depth at the specific variant site ($DP$) |
| `Mean_Coverage` | Overall whole-mitochondrial mean coverage |
| `Major_Haplogroup` | Primary phylogenetic lineage assigned by PhyloTree Build 17 |
| `Minor_Haplogroup` | Secondary lineage detected (indicative of mixture) |
| `Contamination_Status` | Sample contamination verdict (`YES` / `NO`) |
| `Contamination_Pct` | Quantified exogenous mixture fraction (%) |
| `Variant_Allele_Fraction` | Heteroplasmy percentage ($\text{VAF} = \frac{AD}{DP} \times 100\%$) |
| `Origin` | Biological categorization: `Germline` (homoplasmic) vs. `Somatic` (heteroplasmic) |
| `NUMT_Homology` | Overlapping nuclear NUMT contig identifier and pseudogene locus |
| `gnomAD_npg_Status` | Broad Institute gnomAD v3.1 nuclear pseudogene artifact flag |

---

## 8. Benchmark Validation & Analytical Performance

MITOCLIN’s analytical accuracy and NUMT-depletion capabilities were benchmarked across a cohort of **128 clinical patient mitochondrial genomes** sequenced at ultra-deep coverage (mean coverage: **13,787×** depth) at the Centre for DNA Fingerprinting and Diagnostics (CDFD):

* **Orthogonal Clinical Concordance:** Evaluated against confirmed diagnostic cases (LHON $m.11778G>A$, $m.3460G>A$, $m.14484T>C$; Leigh Syndrome $m.8993T>G/C$; MERRF $m.8344A>G$) verified by diagnostic capillary sequencing, achieving 100% concordance.
* **NUMT Pseudogene Purging:** Competitive dual-alignment against 766 nuclear NUMT contigs absorbed **63,230,940 nuclear reads** away from the mitochondrial alignment, purging **527 false NUMT-derived heteroplasmies** (a **55.7% reduction** in false low-frequency variants in the 5%–20% VAF window) while retaining 100% of authentic biological homoplasmies ($\ge 95\%$ VAF).
* **Contamination Sensitivity:** Integrated HaploCheck engine with PhyloTree Build 17 LCA distance reliably detects cross-sample contamination down to a **1.8%** minor haplotype fraction ($F_{\text{minor}} = 0.018$).
* **Limit of Detection (LOD):** Under ultra-deep sequencing coverage ($>5,000\times$), dual-strand read support provides a reliable Limit of Detection of **1.00% VAF**.

---

## 9. Citation & Institutional Affiliation

If you use MITOCLIN in your research, please cite:

> **MITOCLIN: an integrated and reproducible computational framework for comprehensive mitochondrial genome analysis and clinical interpretation from next-generation sequencing data**  
> Laboratory of Genome Informatics (LGI),  
> Centre for DNA Fingerprinting and Diagnostics (BRIC-CDFD), Hyderabad, Telangana, India.  
> Repository: [https://github.com/lgi-cdfd/mitoclin](https://github.com/lgi-cdfd/mitoclin)

For inquiries, clinical collaboration, or bug reports, please contact:
* **Dr. Ajay Kumar Mahato** (Principal Investigator / Corresponding Author): `akmahato@cdfd.org.in`
* **Laboratory of Genome Informatics (LGI)**, BRIC-CDFD, Hyderabad, India.

---

## 10. License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for complete details.
