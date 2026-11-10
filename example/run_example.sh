#!/usr/bin/env bash

# ==============================================================================
# MITOCLIN Example Test Runner
# Executes the pipeline on bundled sample paired-end FASTQ reads
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SAMPLE="TestSample"
R1="$SCRIPT_DIR/sample_R1.fastq.gz"
R2="$SCRIPT_DIR/sample_R2.fastq.gz"
OUTDIR="$REPO_DIR/example_output"

echo "================================================================================"
echo "Running MITOCLIN Validation Test on Example FASTQ Dataset"
echo "Sample ID : $SAMPLE"
echo "R1        : $R1"
echo "R2        : $R2"
echo "Output    : $OUTDIR"
echo "================================================================================"

"$REPO_DIR/mitoclin" \
    --sample "$SAMPLE" \
    --r1 "$R1" \
    --r2 "$R2" \
    --outdir "$OUTDIR" \
    --threads 2

echo ""
echo "================================================================================"
echo "Example run completed successfully."
echo "Results are stored in: $OUTDIR"
echo "Final annotated variants table: $OUTDIR/10_final_report_data/${SAMPLE}_final_report.tsv"
echo "================================================================================"
