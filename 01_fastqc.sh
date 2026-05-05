#!/usr/bin/env bash
# Step 1 — Raw read QC (FastQC + MultiQC aggregation)
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/student_dataset"
QC_OUT="$PROJECT/qc/fastqc_raw"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 1: Raw read QC ==="

mkdir -p "$QC_OUT" "$PROJECT/trimmed" "$PROJECT/aligned" "$PROJECT/peaks"

# build FASTQ list (R1 + R2 for each assigned sample)
FASTQS=()
for sample in "${SAMPLES[@]}"; do
  FASTQS+=("$DATA/${sample}_R1_1M.fastq.gz" "$DATA/${sample}_R2_1M.fastq.gz")
done

echo ">>> Running FastQC on ${#FASTQS[@]} files"
fastqc -t 8 -o "$QC_OUT" "${FASTQS[@]}"

echo ">>> Aggregating with MultiQC"
multiqc "$QC_OUT" -o "$QC_OUT" --filename multiqc_report.html --force

echo ">>> Cleaning per-sample HTML/ZIP (kept by multiqc)"
rm -f "$QC_OUT"/*_fastqc.html "$QC_OUT"/*_fastqc.zip "$QC_OUT"/fastqc_run.log

echo "=== Step 1 complete ==="
echo "Report: $QC_OUT/multiqc_report.html"
