#!/usr/bin/env bash
# Step 10 — Generate RPKM-normalised bigWig coverage tracks for IGV
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/data"
RESULTS="$PROJECT/results"

ALIGNED="$DATA/aligned"
COV="$RESULTS/coverage"
PEAKS="$RESULTS/peaks"
mkdir -p "$COV"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)

echo "=== Step 10: bigWig generation (deepTools bamCoverage) ==="

for sample in "${SAMPLES[@]}"; do
  IN="$ALIGNED/${sample}_filtered_MAPQ20.bam"
  OUT="$COV/${sample}.bw"

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
echo "       $COV/IL7-H_H3k27ac_S295.bw"
echo "       $COV/IL7-H_igG_S286.bw"
echo "       $PEAKS/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak"
echo "  3. Navigate to a gene (e.g. Actb, Gapdh, Cd3e) and screenshot."
