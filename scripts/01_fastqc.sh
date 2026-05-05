#!/usr/bin/env bash
# Step 1 — Raw read QC (FastQC + MultiQC aggregation)
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/data"
RESULTS="$PROJECT/results"

RAW="$DATA/student_dataset"
QC_RAW="$RESULTS/qc/fastqc_raw"
QC_OUT="$RESULTS/qc"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 1: Raw read QC ==="

mkdir -p "$QC_RAW" "$DATA/trimmed" "$DATA/aligned" "$RESULTS/peaks"

# build FASTQ list (R1 + R2 for each assigned sample)
FASTQS=()
for sample in "${SAMPLES[@]}"; do
  FASTQS+=("$RAW/${sample}_R1_1M.fastq.gz" "$RAW/${sample}_R2_1M.fastq.gz")
done

echo ">>> Running FastQC on ${#FASTQS[@]} files"
fastqc -t 8 -o "$QC_RAW" "${FASTQS[@]}"

echo ">>> Aggregating with MultiQC → $QC_OUT/multiqc_report.html"
multiqc "$QC_RAW" -o "$QC_OUT" --filename multiqc_report.html --force

echo ">>> Cleaning per-sample HTML/ZIP (kept by multiqc)"
rm -f "$QC_RAW"/*_fastqc.html "$QC_RAW"/*_fastqc.zip "$QC_RAW"/fastqc_run.log

echo "=== Step 1 complete ==="
echo "Report: $QC_OUT/multiqc_report.html"
