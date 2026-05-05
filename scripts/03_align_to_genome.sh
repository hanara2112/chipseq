#!/usr/bin/env bash
# Step 3 — Align trimmed reads to mm10 with Bowtie2, sort + index
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/data"
RESULTS="$PROJECT/results"
TRIMMED="$DATA/trimmed"
OUT="$DATA/aligned"
INDEX="$DATA/genome/mm10"
mkdir -p "$OUT"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 3: Bowtie2 alignment to mm10 ==="

for sample in "${SAMPLES[@]}"; do
  R1="$TRIMMED/${sample}_R1_1M_val_1.fq.gz"
  R2="$TRIMMED/${sample}_R2_1M_val_2.fq.gz"
  BAM="$OUT/${sample}.bam"
  LOG="$OUT/${sample}_bowtie2.log"

  echo ">>> Aligning $sample"
  bowtie2 \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    --threads 8 \
    -X 2000 \
    --no-mixed \
    --no-discordant \
    2> "$LOG" \
  | samtools sort -@ 4 -o "$BAM"

  samtools index "$BAM"
  echo "  done: $BAM"
  echo "  log:  $LOG"
done

echo "=== Step 3 complete ==="
