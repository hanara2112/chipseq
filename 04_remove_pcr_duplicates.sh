#!/usr/bin/env bash
# Step 4 — Remove PCR duplicates (samtools fixmate + samtools markdup)
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
ALIGNED="$PROJECT/aligned"
QC="$PROJECT/qc"
mkdir -p "$QC"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 4: PCR duplicate removal ==="

for sample in "${SAMPLES[@]}"; do
  IN="$ALIGNED/${sample}.bam"
  NSORT="$ALIGNED/${sample}.nsort.bam"
  FIXED="$ALIGNED/${sample}.fixmate.bam"
  CSORT="$ALIGNED/${sample}.fixmate.sorted.bam"
  DEDUP="$ALIGNED/${sample}_deduplicated.bam"
  METRICS="$QC/${sample}_dup_metrics.txt"

  echo ">>> Processing $sample"

  echo "  [1/4] name-sorting BAM"
  samtools sort -n -@ 4 -o "$NSORT" "$IN"

  echo "  [2/4] fixmate (adds MC/ms tags)"
  samtools fixmate -m -@ 4 "$NSORT" "$FIXED"

  echo "  [3/4] coordinate-sorting"
  samtools sort -@ 4 -o "$CSORT" "$FIXED"

  echo "  [4/4] samtools markdup (-r removes duplicates, -s writes stats)"
  samtools markdup -r -s -@ 4 \
    -f "$METRICS" \
    "$CSORT" \
    "$DEDUP"

  samtools index "$DEDUP"

  echo "  cleaning intermediates"
  rm -f "$NSORT" "$FIXED" "$CSORT"

  echo "  done: $DEDUP"
  echo "  metrics: $METRICS"
done

echo "=== Step 4 complete ==="
