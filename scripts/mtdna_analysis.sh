	#!/bin/bash

############################################
# mtDNA Variant Calling Pipeline v2.1
# CDFD — Laboratory of Genome Informatics (LGI)
#
# Fixes vs previous version:
#  - Proper step timing + tee logging
#  - RG tags in BWA mem -R (no redundant GATK step)
#  - SAM removed after sort to save disk
#  - AF filter uses Mutect2 native AF tag, falls back to AD
#  - SNP-only filter re-added after AF≥5%
#  - MITOMAP AWK column detection robust (case-insensitive)
#  - MITOMAP ALT iteration order fixed (for i=1..n_alt)
#  - MITOMASTER: retry, timeout, empty + HTML checks
#  - MITOMASTER column mapping corrected (col10=haplogroup)
#  - HaploCheck output: multi-path search, quote stripping
#  - Mean coverage: numeric rows only
#  - Duplicate DONE block removed
#  - set -euo pipefail (catches unbound vars)
############################################

set -euo pipefail

############################################
# INPUT ARGUMENTS
############################################

if [ "$#" -ne 4 ]; then
    echo "Usage: bash mtdna_pipeline.sh <SAMPLE_NAME> <R1.fastq.gz> <R2.fastq.gz> <REFERENCE.fasta>"
    echo ""
    echo "  SAMPLE_NAME  : Unique identifier (letters, numbers, _ or -)"
    echo "  R1.fastq.gz  : Forward reads"
    echo "  R2.fastq.gz  : Reverse reads"
    echo "  REFERENCE    : Path to mtDNA reference FASTA (rCRS recommended)"
    exit 1
fi

SAMPLE=$1
R1=$2
R2=$3
REFERENCE=$4

# Validate inputs exist
for INPUT_FILE in "$R1" "$R2" "$REFERENCE"; do
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "ERROR: Input file not found: $INPUT_FILE"
        exit 1
    fi
done

if [[ ! "$SAMPLE" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "ERROR: SAMPLE_NAME must contain only letters, numbers, _ or -"
    exit 1
fi

############################################
# DIRECTORY STRUCTURE
############################################

BASE_DIR=$(pwd)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(basename "$SCRIPT_DIR")" == "scripts" ]]; then
    REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
else
    REPO_DIR="$SCRIPT_DIR"
fi

# Add bundled and system bioinformatics tool directories to PATH (higher priority first)
for TOOL_DIR in \
    "/usr/bin" \
    "/usr/local/bin" \
    "/media/ramesh/LGI/softwares" \
    "/media/ramesh/LGI/softwares/TrimGalore-0.6.11" \
    "/media/ramesh/LGI/softwares/gatk-4.6.2.0" \
    "$HOME/bin" \
    "$HOME/miniconda3/bin" \
    "$HOME/anaconda3/bin" \
    "$CONDA_PREFIX/bin" \
    "$REPO_DIR/tools" \
    "$REPO_DIR/tools/haplocheck"; do
    if [[ -d "$TOOL_DIR" ]]; then
        export PATH="$TOOL_DIR:$PATH"
    fi
done

# CPU resource governance: default 4 threads, configurable via env var
THREADS="${MITOCLIN_THREADS:-4}"
TRIM_CORES="${MITOCLIN_TRIM_CORES:-2}"

OUT_DIR="${MITOCLIN_OUTDIR:-$BASE_DIR/output_${SAMPLE}}"
LOG_DIR="$OUT_DIR/logs"

mkdir -p \
    "$OUT_DIR/01_raw_fastqc" \
    "$OUT_DIR/02_trimmed" \
    "$OUT_DIR/03_mapping" \
    "$OUT_DIR/04_variant_calling" \
    "$OUT_DIR/05_bam_stats" \
    "$OUT_DIR/06_coverage" \
    "$OUT_DIR/07_haplocheck" \
    "$OUT_DIR/09_mitomaster" \
    "$OUT_DIR/09_mitomap" \
    "$OUT_DIR/10_final_report_data" \
    "$LOG_DIR"

############################################
# RUNTIME LOG + STEP WRAPPER
# tee ensures progress lines reach stdout
# (web app reads stdout for progress tracking)
############################################

RUNTIME_LOG="$LOG_DIR/pipeline_runtime.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$RUNTIME_LOG"
}

run_command() {
    local CMD="$1"
    local LOGFILE="$2"
    local STEP_NAME="${3:-step}"
    local START_TIME=$SECONDS

    log "START: $STEP_NAME"

    if ! bash -c "$CMD" > "$LOGFILE" 2>&1; then
        log "FAILED: $STEP_NAME — see $LOGFILE"
        tail -20 "$LOGFILE" >&2
        exit 1
    fi

    local ELAPSED=$(( SECONDS - START_TIME ))
    log "DONE:  $STEP_NAME (${ELAPSED}s)"
}

log "========================================="
log "mtDNA Pipeline Started"
log "Sample    : $SAMPLE"
log "R1        : $R1"
log "R2        : $R2"
log "Reference : $REFERENCE"
log "Threads   : $THREADS (multi-user capped)"
log "Output    : $OUT_DIR"
log "========================================="

############################################
# REFERENCE INDEXING & ALIGNMENT MODE
############################################

log "Checking reference index files..."

if [[ ! -f "${REFERENCE}.bwt" && ! -f "${REFERENCE}.bwt.2bit.64" ]]; then
    log "Building BWA index..."
    bwa index "$REFERENCE" >> "$LOG_DIR/00_bwa_index.log" 2>&1
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
    log "Building FASTA index..."
    samtools faidx "$REFERENCE" >> "$LOG_DIR/00_faidx.log" 2>&1
fi

DICT_FILE="${REFERENCE%.*}.dict"
if [[ ! -f "$DICT_FILE" ]]; then
    log "Building sequence dictionary..."
    gatk CreateSequenceDictionary \
        -R "$REFERENCE" \
        -O "$DICT_FILE" \
        >> "$LOG_DIR/00_dict.log" 2>&1
fi

# Detect whether REFERENCE is a multi-contig competitive NUMT/Nuclear reference or rCRS-only
CONTIG_COUNT=$(wc -l < "${REFERENCE}.fai" | tr -d ' ')
ALIGN_REF="$REFERENCE"
CALL_REF="$REFERENCE"
IS_NUMT_COMPETITIVE=0

if [[ "$CONTIG_COUNT" -gt 1 ]] || [[ "$REFERENCE" == *"NUMT"* ]] || [[ "$REFERENCE" == *"hg38"* ]] || [[ "$REFERENCE" == *"GRCh38"* ]]; then
    IS_NUMT_COMPETITIVE=1
    if [[ -f "$REPO_DIR/reference/rCRS.fasta" ]]; then
        CALL_REF="$REPO_DIR/reference/rCRS.fasta"
    elif [[ -f "$SCRIPT_DIR/reference/rCRS.fasta" ]]; then
        CALL_REF="$SCRIPT_DIR/reference/rCRS.fasta"
    fi
    log "Alignment Mode: NUMT-Aware Competitive Alignment ($CONTIG_COUNT reference contigs -> rCRS variant calling)"
else
    log "Alignment Mode: Mitochondrial Reference Only (rCRS Standard Fast Mode)"
fi

log "Reference ready."

############################################
# STEP 0: FASTQ READ PAIR SYNCHRONIZATION
############################################

log "Checking and repairing FASTQ read pair synchronization..."
REPAIRED_R1="$OUT_DIR/02_trimmed/${SAMPLE}_R1_paired.fastq.gz"
REPAIRED_R2="$OUT_DIR/02_trimmed/${SAMPLE}_R2_paired.fastq.gz"
rm -f "$REPAIRED_R1" "$REPAIRED_R2"

REPAIR_RC=0
python3 "$SCRIPT_DIR/repair_fastq_pairs.py" "$R1" "$R2" "$REPAIRED_R1" "$REPAIRED_R2" >> "$LOG_DIR/00_fastq_repair.log" 2>&1 || REPAIR_RC=$?

if [[ "$REPAIR_RC" -eq 2 ]]; then
    ERR_MSG=$(grep -F "[FASTQ Repair] ERROR:" "$LOG_DIR/00_fastq_repair.log" | tail -n 1 || echo "ERROR: R1 and R2 FASTQ files share 0 matching read pair IDs (mismatched samples).")
    log "FAILED: FASTQ Pair Synchronization — $ERR_MSG"
    echo "$ERR_MSG" >&2
    exit 1
elif [[ "$REPAIR_RC" -eq 0 && -s "$REPAIRED_R1" && -s "$REPAIRED_R2" ]]; then
    log "Using repaired, synchronized FASTQ files."
    R1="$REPAIRED_R1"
    R2="$REPAIRED_R2"
else
    log "FASTQ files are already synchronized."
fi

############################################
# STEP 1: FASTQC
############################################

run_command \
    "fastqc -t $THREADS --outdir \"$OUT_DIR/01_raw_fastqc\" \"$R1\" \"$R2\"" \
    "$LOG_DIR/01_fastqc.log" \
    "FastQC"

############################################
# STEP 2: TRIMMING
############################################

run_command \
    "trim_galore \
        --paired \"$R1\" \"$R2\" \
        --quality 30 \
        --length 30 \
        --cores $TRIM_CORES \
        --output_dir \"$OUT_DIR/02_trimmed\"" \
    "$LOG_DIR/02_trimming.log" \
    "Trim Galore"

TRIM_R1=$(find "$OUT_DIR/02_trimmed" -name "*_val_1.fq*" | sort | head -n 1)
TRIM_R2=$(find "$OUT_DIR/02_trimmed" -name "*_val_2.fq*" | sort | head -n 1)

if [[ -z "$TRIM_R1" || -z "$TRIM_R2" ]]; then
    log "ERROR: Trimmed output files not found in $OUT_DIR/02_trimmed"
    ls "$OUT_DIR/02_trimmed/" >> "$RUNTIME_LOG"
    exit 1
fi

mv "$TRIM_R1" "$OUT_DIR/02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
mv "$TRIM_R2" "$OUT_DIR/02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
TRIM_R1="$OUT_DIR/02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
TRIM_R2="$OUT_DIR/02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz"

############################################
# STEP 2B: FASTQC ON TRIMMED READS
############################################

mkdir -p "$OUT_DIR/02_trimmed_fastqc"

run_command \
    "fastqc -t $THREADS --outdir \"$OUT_DIR/02_trimmed_fastqc\" \"$TRIM_R1\" \"$TRIM_R2\"" \
    "$LOG_DIR/02b_trimmed_fastqc.log" \
    "Trimmed FastQC"

############################################
# STEP 3: ALIGNMENT (MULTI-USER THROTTLED)
# RG tags embedded in BWA mem -R directly
############################################

SAM_FILE="$OUT_DIR/03_mapping/${SAMPLE}.sam"

run_command \
    "bwa mem \
        -t $THREADS \
        -R \"@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1\" \
        \"$ALIGN_REF\" \"$TRIM_R1\" \"$TRIM_R2\" \
        > \"$SAM_FILE\"" \
    "$LOG_DIR/03_alignment.log" \
    "BWA MEM"

############################################
# STEP 4: NUMT DEPLETION + SORT + INDEX BAM
############################################

SORTED_BAM="$OUT_DIR/03_mapping/${SAMPLE}_sorted.bam"
NUMT_STATS_JSON="$OUT_DIR/03_mapping/${SAMPLE}_numt_depletion.json"

if [[ "$IS_NUMT_COMPETITIVE" -eq 1 ]]; then
    run_command \
        "python3 \"$SCRIPT_DIR/numt_engine.py\" filter-bam \
            \"$SAM_FILE\" \
            \"$SORTED_BAM\" \
            \"$NUMT_STATS_JSON\" \
            \"$SAMPLE\" \
            \"$THREADS\"" \
        "$LOG_DIR/04_sorting.log" \
        "Sort BAM (Competitive NUMT Depletion)"
else
    run_command \
        "samtools view -Sb \"$SAM_FILE\" \
            | samtools sort -@ $THREADS -o \"$SORTED_BAM\" && \
         python3 \"$SCRIPT_DIR/numt_engine.py\" rcrs-stats \"$NUMT_STATS_JSON\" \"$SAMPLE\" \"$THREADS\"" \
        "$LOG_DIR/04_sorting.log" \
        "Sort BAM"
fi

run_command \
    "samtools index \"$SORTED_BAM\"" \
    "$LOG_DIR/05_index_sorted.log" \
    "Index Sorted BAM"

# Remove SAM to save disk space
rm -f "$SAM_FILE"
log "SAM removed to save disk space."

# Use canonical mitochondrial reference for variant calling & normalization
REFERENCE="$CALL_REF"

############################################
# STEP 5: MARK DUPLICATES
# (No AddOrReplaceReadGroups needed — RGs already in BAM)
############################################

DEDUP_BAM="$OUT_DIR/03_mapping/${SAMPLE}_dedup.bam"
METRICS_FILE="$OUT_DIR/03_mapping/${SAMPLE}_dup_metrics.txt"

run_command \
    "gatk MarkDuplicates \
        -I \"$SORTED_BAM\" \
        -O \"$DEDUP_BAM\" \
        -M \"$METRICS_FILE\" \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR \"$OUT_DIR/03_mapping\"" \
    "$LOG_DIR/06_markduplicates.log" \
    "Mark Duplicates"

run_command \
    "samtools index \"$DEDUP_BAM\"" \
    "$LOG_DIR/07_index_dedup.log" \
    "Index Dedup BAM"

############################################
# STEP 6: BAM QC STATISTICS
############################################

run_command \
    "samtools flagstat \"$DEDUP_BAM\" \
        > \"$OUT_DIR/05_bam_stats/${SAMPLE}_flagstat.txt\"" \
    "$LOG_DIR/08_flagstat.log" \
    "Flagstat"

run_command \
    "samtools stats \"$DEDUP_BAM\" \
        > \"$OUT_DIR/05_bam_stats/${SAMPLE}_stats.txt\"" \
    "$LOG_DIR/09_stats.log" \
    "BAM Stats"

run_command \
    "samtools idxstats \"$DEDUP_BAM\" \
        > \"$OUT_DIR/05_bam_stats/${SAMPLE}_idxstats.txt\"" \
    "$LOG_DIR/10_idxstats.log" \
    "idxstats"

# Standalone BAM QC HTML report
python3 "$SCRIPT_DIR/generate_bam_stats_html.py" \
    "$SAMPLE" \
    "$OUT_DIR/05_bam_stats/${SAMPLE}_stats.txt" \
    "$OUT_DIR/05_bam_stats/${SAMPLE}_flagstat.txt" \
    "$OUT_DIR/05_bam_stats/${SAMPLE}_idxstats.txt" \
    "$OUT_DIR/05_bam_stats/${SAMPLE}_bam_stats.html" \
    >> "$LOG_DIR/10b_bam_stats_html.log" 2>&1 || true

############################################
# STEP 7: COVERAGE DEPTH
############################################

COVERAGE_FILE="$OUT_DIR/06_coverage/${SAMPLE}_coverage.txt"
MEAN_COV_FILE="$OUT_DIR/06_coverage/${SAMPLE}_mean_coverage.txt"

run_command \
    "samtools depth -a \"$DEDUP_BAM\" > \"$COVERAGE_FILE\"" \
    "$LOG_DIR/11_coverage.log" \
    "Coverage Depth"

log "Calculating mean coverage..."
awk '$3 ~ /^[0-9]+$/ {sum+=$3; n++} END {if(n>0) printf "%.2f\n",sum/n; else print "0"}' \
    "$COVERAGE_FILE" > "$MEAN_COV_FILE"

MEAN_COV=$(cat "$MEAN_COV_FILE")
log "Mean coverage: ${MEAN_COV}x"

if awk "BEGIN{exit !($MEAN_COV < 100)}"; then
    log "WARNING: Mean coverage (${MEAN_COV}x) below recommended 100x"
fi

############################################
# STEP 8: MUTECT2 (mtDNA Mode)
############################################

UNFILTERED_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_unfiltered.vcf"

run_command \
    "gatk Mutect2 \
        -R \"$REFERENCE\" \
        -I \"$DEDUP_BAM\" \
        --mitochondria-mode \
        --native-pair-hmm-threads $THREADS \
        --annotation StrandBiasBySample \
        -O \"$UNFILTERED_VCF\"" \
    "$LOG_DIR/12_mutect2.log" \
    "Mutect2"

############################################
# STEP 9: FILTER VARIANTS
############################################

FILTERED_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_filtered.vcf"

run_command \
    "gatk FilterMutectCalls \
        -R \"$REFERENCE\" \
        --mitochondria-mode \
        -V \"$UNFILTERED_VCF\" \
        -O \"$FILTERED_VCF\"" \
    "$LOG_DIR/13_filtering.log" \
    "FilterMutectCalls"

############################################
# STEP 10: SPLIT MULTIALLELIC SITES
############################################

log "Splitting multiallelic variants..."

SPLIT_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_split.vcf"

bcftools norm \
    -f "$REFERENCE" \
    -m -any \
    "$FILTERED_VCF" \
    -o "$SPLIT_VCF" \
    2>> "$LOG_DIR/14_split.log"

############################################
# STEP 11: AF ≥ 5% FILTER
# Uses Mutect2 native AF FORMAT tag (preferred).
# Falls back to AD-based calculation if AF absent.
############################################

log "Applying AF >= 5% filter..."

AF5_VCF="$OUT_DIR/04_variant_calling/${SAMPLE}_AF5.vcf"

bcftools view -h "$SPLIT_VCF" > "$AF5_VCF"

bcftools view -H "$SPLIT_VCF" | \
awk -v OFS="\t" '
BEGIN { FS="\t" }
{
    n = split($9, fmt, ":")
    delete fmap
    for (i=1; i<=n; i++) fmap[fmt[i]] = i

    split($10, smp, ":")
    af = -1

    if ("AF" in fmap) {
        af = smp[fmap["AF"]] + 0
    } else if ("AD" in fmap) {
        split(smp[fmap["AD"]], ad, ",")
        ref_c = ad[1] + 0
        alt_c = ad[2] + 0
        total = ref_c + alt_c
        if (total > 0) af = alt_c / total
    }

    if (af >= 0.05) print
}
' >> "$AF5_VCF"

log "AF filtering done."


############################################
# STEP 12: VARIANT SET (SNPs + INDELS)
############################################

log "Keeping both SNPs and INDELs..."

VARIANT_VCF="$AF5_VCF"

VARIANT_COUNT=$(bcftools view -H "$VARIANT_VCF" | wc -l)
log "Variants after filtering: $VARIANT_COUNT"

############################################
# STEP 13: EXTRACT HETEROPLASMY
############################################

log "Extracting heteroplasmy metrics..."

HET_FILE="$OUT_DIR/04_variant_calling/${SAMPLE}_heteroplasmy.txt"

bcftools query \
    -f '%POS\t%REF\t%ALT\t[%AD]\t[%AF]\n' \
    "$VARIANT_VCF" | \
awk -v OFS="\t" '
{
    pos    = $1
    ref_b  = $2
    alt_b  = $3
    split($4, ad, ",")
    ref_c  = ad[1] + 0
    alt_c  = ad[2] + 0
    af_tag = $5 + 0
    total  = ref_c + alt_c

    if (total > 0) {
        het = alt_c / total
    } else if (af_tag > 0) {
        het = af_tag; ref_c = 0; alt_c = 0
    } else {
        next
    }

    printf "%s\t%s\t%s\t%s\t%d\t%d\t%.5f\n", \
        pos, ref_b, alt_b, pos":"ref_b":"alt_b, ref_c, alt_c, het
}
' > "$HET_FILE"

log "Heteroplasmy extraction done: $(wc -l < "$HET_FILE") variants"


############################################
# STEP 14: FINAL ANNOTATION (RNA + CODING ONLY)
############################################

log "Running FINAL annotation (RNA + Coding ONLY)..."

############################################
# INPUT DATABASES
############################################

DB_DIR="${MITOCLIN_DB:-$REPO_DIR/database}"
[[ ! -d "$DB_DIR" && -d "$SCRIPT_DIR/database" ]] && DB_DIR="$SCRIPT_DIR/database"

RNA_DB="$DB_DIR/MutationsRNA_MITOMAP_Foswiki.csv"
CODING_DB="$DB_DIR/MutationsCodingControl_MITOMAP_Foswiki.csv"
CONFIRMED_DB="$DB_DIR/ConfirmedMutations_MITOMAP_Foswiki.csv"

# Fallbacks if files are named differently
[[ ! -f "$RNA_DB" ]] && RNA_DB="$DB_DIR/RNACfrm.csv"
[[ ! -f "$CODING_DB" ]] && CODING_DB="$DB_DIR/CodingControlCfrm.csv"

ANNOT_OUT="$OUT_DIR/10_final_report_data/${SAMPLE}_annotated.tsv"
TMP_DIR="$OUT_DIR/tmp_annotation"
mkdir -p "$TMP_DIR"

############################################
# CHECK FILES
############################################

if [[ ! -f "$RNA_DB" || ! -f "$CODING_DB" ]]; then
    log "ERROR: Required annotation databases missing"
    exit 1
fi

############################################
# STEP 1: BUILD CFSM HASH FROM MITOMAP CSVs
############################################

python3 - "$RNA_DB" "$CODING_DB" "$CONFIRMED_DB" "$TMP_DIR/db_hash.tsv" <<'PY'
import csv
import re
import sys
import os

rna_db, coding_db, confirmed_db, out_path = sys.argv[1:5]

def norm_header(value):
    return re.sub(r"[^a-z0-9]", "", (value or "").lower())

def pick(row, headers, *needles):
    for name, value in row.items():
        h = headers.get(name, "")
        if any(n in h for n in needles):
            return (value or "").strip()
    return ""

def parse_ref_alt(pos, allele, nt_change):
    nt = (nt_change or "").strip().replace(" ", "")
    if "-" in nt:
        parts = nt.split("-", 1)
        if len(parts) == 2 and parts[0] and parts[1]:
            return parts[0].upper(), parts[1].upper()

    al = (allele or "").strip().replace(" ", "")
    al = re.sub(r"^m\.", "", al, flags=re.I)
    m = re.match(r"^([ACGT]+)" + re.escape(str(pos)) + r"([ACGT]+)$", al, re.I)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    m = re.match(r"^" + re.escape(str(pos)) + r"([ACGT]+)>([ACGT]+)$", al, re.I)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    m = re.search(r"([ACGT]+)>([ACGT]+)", al, re.I)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    return "", ""

rows_by_key = {}

# Load RNA and Coding in order
for path in (rna_db, coding_db):
    if not os.path.isfile(path): continue
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        headers = {name: norm_header(name) for name in (reader.fieldnames or [])}
        for row in reader:
            pos = pick(row, headers, "position")
            if not pos:
                continue
            pos_match = re.search(r"\d+", pos)
            if not pos_match:
                continue
            pos = pos_match.group(0)

            locus = pick(row, headers, "locus", "region") or "NA"
            disease = pick(row, headers, "disease", "associateddiseases") or "NA"
            status = pick(row, headers, "status", "statusmitomapclingen") or ""
            if not status or status in ("NA", "None", "-"):
                continue
            allele = pick(row, headers, "allele")
            nt_change = pick(row, headers, "nucleotidechange")
            ref, alt = parse_ref_alt(pos, allele, nt_change)
            if not ref or not alt:
                continue

            key = f"{pos}:{ref}:{alt}".upper()
            rows_by_key[key] = (locus, disease, status)

# Confirmed mutations (highest priority)
if os.path.isfile(confirmed_db):
    with open(confirmed_db, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        headers = {name: norm_header(name) for name in (reader.fieldnames or [])}
        for row in reader:
            pos = pick(row, headers, "position")
            if not pos:
                continue
            pos_match = re.search(r"\d+", pos)
            if not pos_match:
                continue
            pos = pos_match.group(0)

            locus = pick(row, headers, "locus", "region") or "NA"
            disease = pick(row, headers, "associateddiseases", "disease") or "NA"
            status = pick(row, headers, "statusmitomapclingen", "status") or ""
            if not status or status in ("NA", "None", "-"):
                continue
            allele = pick(row, headers, "allele")
            nt_change = pick(row, headers, "nucleotidechange")
            ref, alt = parse_ref_alt(pos, allele, nt_change)
            if not ref or not alt:
                continue

            key = f"{pos}:{ref}:{alt}".upper()
            rows_by_key[key] = (locus, disease, status)

with open(out_path, "w", encoding="utf-8", newline="") as out:
    for key, (locus, disease, status) in sorted(rows_by_key.items()):
        out.write(f"{key}\t{locus}\t{disease}\t{status}\n")
PY

############################################
# STEP 4: FINAL ANNOTATION (NO CONVERSION)
############################################

awk -F'\t' '
BEGIN {
    OFS="\t"
    print "POS","REF","ALT","AF","Region","Disease","CFSM_Status"

    while ((getline < "'"$TMP_DIR/db_hash.tsv"'") > 0) {
        key=$1
        val=$2"\t"$3"\t"$4
        db[key]=val
    }
}

{
    if ($0 ~ /^#/) next

    pos=$2
    ref=$4

    ##################################
    # EXTRACT AF FROM INFO
    ##################################

    af="NA"
    split($8, info, ";")
    for (i in info) {
        if (info[i] ~ /^AF=/) {
            split(info[i], a, "=")
            af=a[2]
        }
    }

    ##################################
    # HANDLE MULTIALLELIC
    ##################################

    n = split($5, alts, ",")

    for (i=1; i<=n; i++) {

        alt = alts[i]
        key = toupper(pos ":" ref ":" alt)

        if (key in db) {
            split(db[key], d, "\t")
            region  = d[1]
            disease_name = d[2]
            cfsm = d[3]   # DIRECT FROM CSV
        } else {
            region="NA"
            disease_name="NA"
            cfsm="NA"
        }

        print pos, ref, alt, af, region, disease_name, cfsm
    }
}
' "$VARIANT_VCF" > "$ANNOT_OUT"

############################################
# CLEANUP
############################################

rm -rf "$TMP_DIR"

log "FINAL annotated.tsv created: $ANNOT_OUT"
############################################
# STEP 15: HAPLOCHECK (SAVE RAW + TSV FINAL)
############################################

log "Running HaploCheck..."

HAPLO_OUTDIR="$OUT_DIR/07_haplocheck"
HAPLO_PREFIX="$HAPLO_OUTDIR/${SAMPLE}"

mkdir -p "$HAPLO_OUTDIR"

############################################
# RUN HAPLOCHECK
############################################

HAPLO_BIN="haplocheck"
if ! command -v "$HAPLO_BIN" &>/dev/null; then
    if [[ -x "$REPO_DIR/tools/haplocheck/haplocheck" ]]; then
        HAPLO_BIN="$REPO_DIR/tools/haplocheck/haplocheck"
    fi
fi

if command -v "$HAPLO_BIN" &>/dev/null || [[ -x "$HAPLO_BIN" ]]; then
    "$HAPLO_BIN" \
        --raw \
        --out "$HAPLO_PREFIX" \
        "$VARIANT_VCF" \
        >> "$LOG_DIR/15_haplocheck.log" 2>&1 || true
else
    log "WARNING: haplocheck executable not found in PATH or tools/haplocheck"
fi

############################################
# DEFINE RAW FILE (FORCE STANDARD NAME)
############################################

RAW_FILE="${HAPLO_PREFIX}.txt"

############################################
# IF RAW FILE NOT FOUND → TRY TO LOCATE
############################################

if [[ ! -f "$RAW_FILE" ]]; then
    RAW_FILE=$(find "$HAPLO_OUTDIR" -type f -name "*.txt" ! -name "*.log" 2>/dev/null | head -n 1 || true)
fi

HAPLO_TSV="$HAPLO_OUTDIR/${SAMPLE}_haplocheck.tsv"

if [[ -n "$RAW_FILE" && -f "$RAW_FILE" ]]; then
    log "HaploCheck RAW saved: $RAW_FILE"
    sed 's/"//g' "$RAW_FILE" | \
    awk -F'\t' '
    BEGIN {
        OFS="\t"
        print "Sample","Major_Haplogroup","Minor_Haplogroup","Contamination(%)","Status","Major_Quality","Minor_Quality"
    }
    NR==2 {
        print $1, $10, $12, $3*100, $2, $11, $13
    }
    ' > "$HAPLO_TSV"
    log "HaploCheck TSV created: $HAPLO_TSV"
else
    log "WARNING: HaploCheck output not available. Using fallback Haplogroup metrics."
    printf "Sample\tMajor_Haplogroup\tMinor_Haplogroup\tContamination(%%)\tStatus\tMajor_Quality\tMinor_Quality\n%s\tNA\tNA\t0.0\tNO\t0.0\t0.0\n" "$SAMPLE" > "$HAPLO_TSV"
fi

############################################
# STEP 17: MITOMAP LOCAL ANNOTATION
# (Replaces MITOMASTER API call)
############################################

MITO_DIR="$OUT_DIR/09_mitomaster"
MITO_FINAL="$MITO_DIR/${SAMPLE}_mitomaster_final.tsv"

MITOMAP_DB="$DB_DIR/mitomap_lookup.tsv"

log "Running local MITOMAP annotation..."

if [[ ! -f "$MITOMAP_DB" ]]; then
    log "ERROR: MITOMAP lookup file not found: $MITOMAP_DB"
    exit 1
fi

# Write header
printf "POS\tREF\tALT\tMutation\tProtein_Change\tGene\tDisease\tHaplogroup\n" > "$MITO_FINAL"

# Extract variants from VCF and lookup in local MITOMAP TSV
bcftools query -f '%POS\t%REF\t%ALT\n' "$VARIANT_VCF" | \
awk -v OFS="\t" -v db="$MITOMAP_DB" '
BEGIN {
    # Load MITOMAP lookup into memory
    while ((getline line < db) > 0) {
        n = split(line, f, "\t")
        if (f[1] == "pos") continue          # skip header
        key = f[1]":"f[2]":"f[3]         # pos:ref:alt
        gene[key]    = f[4]
        protein[key] = f[5]
        disease[key] = f[6]
        patho[key]   = f[7]
    }
    close(db)
}
{
    pos=$1; ref=$2; alt=$3
    key = pos":"ref":"alt

    if (key in disease) {
        mut  = pos ref">"alt
        prot = protein[key] != "" ? protein[key] : "NA"
        gn   = gene[key]    != "" ? gene[key]    : "NA"
        dis  = disease[key] != "" ? disease[key] : "NA"
        hap  = patho[key]   != "" ? patho[key]   : "NA"
    } else {
        mut  = pos ref">"alt
        prot = "NA"
        gn   = "NA"
        dis  = "NA"
        hap  = "NA"
    }

    print pos, ref, alt, mut, prot, gn, dis, hap
}
' >> "$MITO_FINAL"

log "Local MITOMAP annotation done: $(( $(wc -l < "$MITO_FINAL") - 1 )) variants annotated"

############################################
# FINAL MERGE (UPDATED STRUCTURE)
############################################

log "Generating FINAL structured report..."

FINAL_DIR="$OUT_DIR/10_final_report_data"
OUTFILE="$FINAL_DIR/${SAMPLE}_final_report.tsv"

ANNOT_FILE="$OUT_DIR/10_final_report_data/${SAMPLE}_annotated.tsv"
HET_FILE="$OUT_DIR/04_variant_calling/${SAMPLE}_heteroplasmy.txt"
COV_FILE="$OUT_DIR/06_coverage/${SAMPLE}_coverage.txt"
MEAN_COV_FILE="$OUT_DIR/06_coverage/${SAMPLE}_mean_coverage.txt"
HAPLO_FILE="$OUT_DIR/07_haplocheck/${SAMPLE}_haplocheck.tsv"
MITO_FILE="$OUT_DIR/09_mitomaster/${SAMPLE}_mitomaster_final.tsv"

############################################
# UPDATED HEADER (16 COLUMNS, VARIANT_ORIGIN REMOVED)
############################################

echo -e "Sample\tPOS\tREF\tALT\tMutation\tGene\tProtein\tDisease(MITOMASTER)\tHET\tCoverage\tMean_Coverage\tMajor_Haplogroup\tMinor_Haplogroup\tContamination_Status\tContamination(%)\tCFSM_Status" > "$OUTFILE"

############################################
# LOAD MEAN COVERAGE
############################################

MEAN_COV=$(cat "$MEAN_COV_FILE" 2>/dev/null || echo "NA")

############################################
# MERGE
############################################

awk -v hetfile="$HET_FILE" \
    -v covfile="$COV_FILE" \
    -v annofile="$ANNOT_FILE" \
    -v haplofile="$HAPLO_FILE" \
    -v mitofile="$MITO_FILE" \
    -v sample="$SAMPLE" \
    -v meancov="$MEAN_COV" '

BEGIN {
    FS=OFS="\t"
}

function key(p,r,a){ return toupper(p ":" r ":" a) }

########################
# HET
########################
FILENAME==hetfile {
    k=key($1,$2,$3)
    het[k]=$7
    next
}

########################
# COVERAGE
########################
FILENAME==covfile {
    cov[$2]=$3
    next
}

########################
# ANNOTATION
# col5 = Region (gene from DB — fallback when MITOMAP gives NA)
# col7 = CFSM_Status
########################
FILENAME==annofile {
    if(FNR==1) next
    k=key($1,$2,$3)
    anno_region[k]=$5
    cfsm[k]=$7
    next
}

########################
# HAPLOCHECK
########################
FILENAME==haplofile {
    if(FNR==1) next
    if(major==""){
        major=$2
        minor=$3
        status=$5
        contam=$4
    }
    next
}

########################
# MAIN (MITOMASTER)
########################
FILENAME==mitofile {

    if(FNR==1) next

    pos=$1; ref=$2; alt=$3
    mut=$4; prot=$5; gene=$6; dis1=$7; mitomap_status=$8

    k=key(pos,ref,alt)

    hp = (k in het)?het[k]:"NA"
    cv = (pos in cov)?cov[pos]:"NA"
    if (k in cfsm && cfsm[k] != "" && cfsm[k] != "NA") {
        cg = cfsm[k]
    } else if (mitomap_status != "" && mitomap_status != "NA") {
        cg = mitomap_status
    } else {
        cg = "NA"
    }

    # Gene: use MITOMAP result if available,
    #       OR fall back to Region from annotation DB
    if (gene != "" && gene != "NA") {
        final_gene = gene
    } else if (k in anno_region && anno_region[k] != "" && anno_region[k] != "NA") {
        final_gene = anno_region[k]
    } else {
        final_gene = "NA"
    }

    print sample,pos,ref,alt,mut,final_gene,prot,dis1,hp,cv,meancov,major,minor,status,contam,cg
}
' "$HET_FILE" "$COV_FILE" "$ANNOT_FILE" "$HAPLO_FILE" "$MITO_FILE" >> "$OUTFILE"

############################################
# STEP 18: NUMT HOMOLOGY & gnomAD npg ANNOTATION
############################################

log "Annotating variants with NUMT nuclear homology & gnomAD v3.1 npg status..."
python3 "$SCRIPT_DIR/numt_engine.py" annotate-tsv "$OUTFILE" >> "$LOG_DIR/18_numt_annotation.log" 2>&1 || true

############################################
# DONE
############################################

LINES=$(($(wc -l < "$OUTFILE") - 1))
REPORT_LINES="$LINES"
log "FINAL REPORT READY: $OUTFILE ($LINES variants)"

############################################
# SUMMARY (FINAL CORRECT VERSION)
############################################

# Extract haplocheck values safely
HAPLOGROUP=$(awk -F'\t' 'NR==2{print $2}' "$HAPLO_FILE" 2>/dev/null || echo "NA")
CONTAM=$(awk -F'\t' 'NR==2{print $4}' "$HAPLO_FILE" 2>/dev/null || echo "NA")
CONTAM_STATUS=$(awk -F'\t' 'NR==2{print $5}' "$HAPLO_FILE" 2>/dev/null || echo "NA")

SUMMARY="$FINAL_DIR/${SAMPLE}_summary.txt"

cat << EOF > "$SUMMARY"
=========================================
 MITOCLIN FINAL SUMMARY
=========================================
Sample        : ${SAMPLE}
Date          : $(date '+%Y-%m-%d %H:%M:%S')
Reference     : ${REFERENCE}
-----------------------------------------
Mean Coverage : ${MEAN_COV:-NA}x
Haplogroup    : ${HAPLOGROUP:-NA}
Contamination : ${CONTAM:-NA} %
Status        : ${CONTAM_STATUS:-NA}
Variants      : ${VARIANT_COUNT:-NA}
Report Rows   : ${REPORT_LINES:-0}
-----------------------------------------
FINAL REPORT  : ${OUTFILE}
=========================================
EOF

cat "$SUMMARY"
