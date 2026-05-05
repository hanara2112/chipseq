#!/usr/bin/env bash
# Step 6 — MACS2 peak calling
#   (a) H3K27ac (broad mark) using IgG as control
#   (b) IgG self-call (no control) → noise baseline for Step 7
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
ALIGNED="$PROJECT/aligned"
PEAKS="$PROJECT/peaks"
mkdir -p "$PEAKS/H3K27ac" "$PEAKS/IgG"

CHIP="IL7-H_H3k27ac_S295"
CTRL="IL7-H_igG_S286"

CHIP_BAM="$ALIGNED/${CHIP}_filtered_MAPQ20.bam"
CTRL_BAM="$ALIGNED/${CTRL}_filtered_MAPQ20.bam"

echo "=== Step 6a: MACS2 — H3K27ac (broad) vs IgG ==="
macs2 callpeak \
  -t "$CHIP_BAM" \
  -c "$CTRL_BAM" \
  -f BAMPE \
  -g mm \
  --broad \
  -q 0.05 \
  -n "$CHIP" \
  --outdir "$PEAKS/H3K27ac" \
  2> "$PEAKS/H3K27ac/${CHIP}_macs2.log"

echo "=== Step 6b: MACS2 — IgG self-call (no control) ==="
macs2 callpeak \
  -t "$CTRL_BAM" \
  -f BAMPE \
  -g mm \
  -q 0.05 \
  -n "$CTRL" \
  --outdir "$PEAKS/IgG" \
  2> "$PEAKS/IgG/${CTRL}_macs2.log"

echo "=== Step 6 complete ==="
echo "H3K27ac peaks: $PEAKS/H3K27ac/${CHIP}_peaks.broadPeak"
echo "IgG peaks:     $PEAKS/IgG/${CTRL}_peaks.narrowPeak"
