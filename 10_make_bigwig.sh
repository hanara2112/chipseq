#!/usr/bin/env bash
# Step 10 — Generate RPKM-normalised bigWig coverage tracks for IGV
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
ALIGNED="$PROJECT/aligned"
QC="$PROJECT/qc"
mkdir -p "$QC"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 10: bigWig generation (deepTools bamCoverage) ==="

for sample in "${SAMPLES[@]}"; do
  IN="$ALIGNED/${sample}_filtered_MAPQ20.bam"
  OUT="$QC/${sample}.bw"

  echo ">>> $sample"
  bamCoverage \
    -b "$IN" \
    -o "$OUT" \
    --normalizeUsing RPKM \
    --binSize 10 \
    --extendReads \
    -p 8
  echo "  done: $OUT"
done

echo "=== Step 10 complete ==="
echo ""
echo "Load in IGV:"
echo "  1. Genome: mm10"
echo "  2. File → Load from File:"
echo "       $QC/IL7-H_H3k27ac_S295.bw"
echo "       $QC/IL7-H_igG_S286.bw"
echo "       $PROJECT/peaks/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak"
echo "  3. Navigate to a gene (e.g. Actb, Gapdh, Cd3e) and screenshot."
