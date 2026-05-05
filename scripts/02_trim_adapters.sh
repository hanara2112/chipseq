#!/usr/bin/env bash
# Step 2 — Adapter + quality trimming (Trim Galore, Q20, paired)
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/data"
RESULTS="$PROJECT/results"

RAW="$DATA/student_dataset"
OUT="$DATA/trimmed"
mkdir -p "$OUT"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 2: Adapter trimming (Q20, paired-end) ==="

for sample in "${SAMPLES[@]}"; do
  R1="$RAW/${sample}_R1_1M.fastq.gz"
  R2="$RAW/${sample}_R2_1M.fastq.gz"

  echo ">>> Trimming $sample"
  trim_galore \
    --paired \
    --cores 4 \
    --quality 20 \
    --fastqc \
    -o "$OUT" \
    "$R1" "$R2"
done

echo "=== Step 2 complete ==="
echo "Outputs: $OUT/*_val_1.fq.gz, *_val_2.fq.gz, *_trimming_report.txt"
