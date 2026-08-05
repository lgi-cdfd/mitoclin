<div align="center">

# MITOCLIN
### An Integrated Web Server & Automated Pipeline for High-Resolution Mitochondrial Genome Analysis, Heteroplasmy Quantification, and Clinical Pathogenicity Annotation

[![License: MIT](https://img.shields.io/badge/License-MIT-teal.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Flask 3.0](https://img.shields.io/badge/Framework-Flask_3.0-green.svg)](https://flask.palletsprojects.com/)
[![Reference rCRS](https://img.shields.io/badge/Reference-rCRS_NC__012920.1-purple.svg)](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)
[![Institution](https://img.shields.io/badge/BRIC--CDFD-Hyderabad_India-orange.svg)](https://cdfd.org.in/)

</div>

---

## Table of Contents
- [Overview](#-overview)
- [Key Features & Scientific Capabilities](#-key-features--scientific-capabilities)
- [Platform Architecture](#-platform-architecture)
- [Installation & Setup](#-installation--setup)
- [Usage Modes](#-usage-modes)
  - [1. Web Server Mode (GUI)](#1-web-server-mode-gui)
  - [2. Command-Line Interface (CLI / Batch Mode)](#2-command-line-interface-cli--batch-mode)
  - [3. Standalone Bioinformatic Pipeline (Bash)](#3-standalone-bioinformatic-pipeline-bash)
- [Directory & File Structure](#-directory--file-structure)
- [Output Reports & Visualizations](#-output-reports--visualizations)
- [Biological Validation & Literature Framework](#-biological-validation--literature-framework)
- [License](#-license)
- [Citation & Contact](#-citation--contact)

---

## Overview

**MITOCLIN** is an end-to-end, automated web application and computational pipeline developed at the **Laboratory of Human Molecular Genetics (LHMG), Centre for DNA Fingerprinting and Diagnostics (CDFD), Hyderabad, India**.

It bridges the gap between raw high-throughput Next-Generation Sequencing (NGS) data and actionable clinical insights by integrating sequence quality control, reference alignment, high-sensitivity variant calling, heteroplasmy quantification, contamination assessment, variant origin classification, disease annotation, and interactive clinical reporting into a unified workflow.

---

## Key Features & Scientific Capabilities

1. **Dual Somatic vs. Germline Variant Origin Engine**:
   - **Germline Origin ($\text{VAF} \ge 95\%$ or Population Polymorphisms)**: Classifies inherited homoplasmic/high-grade alleles defining maternal haplogroups or ancestral traits.
   - **Somatic Origin ($5\% \le \text{VAF} < 95\%$)**: Identifies post-zygotic heteroplasmic mutations acquired during tissue aging, oxidative stress, or tumorigenesis.

2. **High-Sensitivity Variant Calling**:
   - Powered by **GATK Mutect2** in mitochondrial mode (`--mitochondria-mode`), enabling detection of low-frequency heteroplasmy down to **5% Allele Fraction (VAF)** without homoplasmy bias.

3. **LCA-Based Contamination Detection**:
   - Integrated **HaploCheck (v1.3.3)** identifies sample contamination using Lowest Common Ancestor (LCA) phylogenetic distance across PhyloTree Build 17, discriminating true heteroplasmy from exogenous DNA mixtures.

4. **Offline Clinical Database Annotation**:
   - Direct matching against a locally curated **MITOMAP** database to report gene names, nucleotide/protein modifications, disease associations, and pathogenicity tiers (Pathogenic, Likely Pathogenic, VUS, Benign).

5. **100% Parity HTML & PDF Diagnostic Reports**:
   - Generates publication-ready, multi-page HTML dashboards and PDF clinical reports featuring institutional branding, executive stat cards, depth plots, heteroplasmy histograms, and 12-column variant summary tables.

---

## Platform Architecture

```
[Raw FASTQ / VCF Input]
          │
          ▼
   ┌──────────────┐
   │ FastQC 0.11  │  --> Read Quality Check
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ Trim Galore  │  --> Adapter Trimming & Quality Filtering (Q >= 20)
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ BWA-MEM      │  --> rCRS NC_012920.1 Alignment & MarkDuplicates
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ GATK Mutect2 │  --> Mitochondria-Mode Variant Calling & Filtering
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ BCFtools     │  --> Heteroplasmy Fraction: VAF = AD / (RD + AD)
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ HaploCheck   │  --> Major/Minor Haplogroups & LCA Contamination
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ MITOMAP Local│  --> Disease Annotation & Somatic/Germline Origin
   └──────┬───────┘
          ▼
   ┌──────────────┐
   │ HTML & PDF   │  --> Executive Diagnostic Reports & VCF Outputs
   └──────────────┘
```

---

## Installation & Setup

### Prerequisites
- **Operating System**: Linux (Ubuntu 20.04/22.04 LTS recommended)
- **Python**: Version 3.10 or higher
- **Bioinformatics Command-Line Tools**: `fastqc`, `trim_galore`, `bwa`, `samtools`, `bcftools`, `gatk` (v4.2.6+), `haplocheck` (v1.3.3)

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/lgi-cdfd/mitoclin.git
cd mitoclin/

# Install Python dependencies
pip install -r requirements.txt
```

---

## Usage Modes

### 1. Web Server Mode (GUI)
Host the full interactive web platform for researchers and clinicians:

```bash
# Direct Python start
python3 mitoclin_app.py --host 0.0.0.0 --port 5005

# Or via shell wrapper script
cd mitoclin_server/
bash run_server.sh 5005
```
Access the web dashboard in your browser at `http://<SERVER_IP>:5005`.

---

### 2. Command-Line Interface (CLI / Batch Mode)
Generate publication-ready HTML and PDF clinical reports directly from a VCF or pipeline TSV file without launching the web server:

```bash
# Generate reports from a VCF file:
python3 mitoclin_app.py --report-only \
  --sample SAMPLE_01 \
  --vcf /path/to/sample.vcf \
  --outdir /path/to/output_dir/

# Generate reports from pipeline final_report.tsv:
python3 mitoclin_app.py --report-only \
  --sample SAMPLE_01 \
  --csv /path/to/sample_final_report.tsv \
  --outdir /path/to/output_dir/
```

---

### 3. Standalone Bioinformatic Pipeline (Bash)
Execute the complete bioinformatic workflow on raw FASTQ files via command line:

```bash
bash mtdna_analysis.sh \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -s SAMPLE_01 \
  -o output_SAMPLE_01/
```

---

## Directory & File Structure

```
mitoclin/
├── mitoclin_app.py            # Flask web application & HTML/PDF report engine
├── mtdna_analysis.sh          # Core bioinformatic execution pipeline
├── requirements.txt           # Python dependencies
├── run_server.sh              # Web server startup script
├── draft_mitoclin_paper.docx  # Publication manuscript document
├── database/                  # Curated MITOMAP & polymorphism databases
├── reference/                 # rCRS FASTA reference & index files
├── figures/                   # High-resolution 300 DPI publication figures
└── output_*/                  # Sample-specific analytical output directories
```

---

## Output Reports & Visualizations

Every completed analysis produces:
1. **Interactive HTML Clinical Report** (`<sample>_mitoclin_report.html`)
2. **Publication-Grade PDF Clinical Report** (`<sample>_mitoclin_report.pdf`)
3. **Filtered Variant Call Format** (`<sample>_final_variants.vcf`)
4. **Annotated Summary Table** (`<sample>_final_report.tsv`)

### Executive Report Sections:
- **Header & Institutional Branding**: BRIC-CDFD logo and metadata bar.
- **Executive Stat Cards**: Mean Coverage, Pathogenic, Likely Pathogenic, VUS, Benign, Germline, Somatic, and Total Variants.
- **Quality Control Summary**: Coverage depth, HaploCheck contamination score, major/minor haplogroups, heteroplasmy cutoff (5% AF), and GATK Mutect2 status.
- **Clinical Interpretation Summary**: Haplogroup badge and automated clinical text narrative.
- **Visual Plots**: 16,569 bp rCRS coverage depth chart, VAF frequency distribution histogram, and pathogenicity pie chart.
- **12-Column Variant Summary Table**: POS, REF, ALT, Mutation, Type, Gene, Protein/RNA, Disease, HET %, Pathogenicity, Origin, Coverage.

---

## Biological Validation & Literature Framework

MITOCLIN algorithms adhere strictly to published clinical guidelines:
- **ACMG/AMP Guidelines**: McCormick, E. M., et al. (2020). *Specifications of the ACMG/AMP variant interpretation guidelines for mitochondrial DNA variants*. Human Mutation, 41(12), 2028–2057.
- **HaploCheck Contamination**: Weissensteiner, H., et al. (2021). *HaploCheck: contamination detection in sequencing data using the mitochondrial phylogenetic tree*. Genome Medicine, 13(1), 56.
- **MITOMAP Database**: Lott, M. T., et al. (2013). *mtDNA Variation and Analysis Using MITOMAP and MITOMASTER*. Current Protocols in Bioinformatics, 44(1), 1.23.1–1.23.26.

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details. Free for academic, clinical, and research use.

---

## Citation & Contact

If you use **MITOCLIN** in your research or clinical diagnostics, please cite:

> Eerapagula, R., Yadav, A., & Mahato, A. K. (2026). *MITOCLIN: An Integrated Web Server and Automated Pipeline for High-Resolution Mitochondrial Genome Analysis, Heteroplasmy Quantification, and Clinical Pathogenicity Annotation*. Laboratory of Genome Informatics, Centre for DNA Fingerprinting and Diagnostics (CDFD), Hyderabad, India.

- **Institution**: Centre for DNA Fingerprinting and Diagnostics (CDFD), Hyderabad, Telangana, India
- **Repository**: [github.com/lgi-cdfd/mitoclin](https://github.com/lgi-cdfd/mitoclin)
- **Contact**: Dr. Ajay Kumar Mahato (`akmlgi2021@gmail.com`)
