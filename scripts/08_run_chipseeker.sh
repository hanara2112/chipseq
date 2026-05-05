#!/usr/bin/env bash
# Step 8 — wrapper that runs the ChIPseeker R script
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
RESULTS="$PROJECT/results"
mkdir -p "$RESULTS/qc"

echo "=== Step 8: ChIPseeker peak annotation ==="
Rscript "$PROJECT/scripts/08_chipseeker_annotate.R"
echo "=== Step 8 complete ==="
echo "Outputs:"
echo "  $RESULTS/peaks/H3K27ac/H3K27ac_annotated.csv"
echo "  $RESULTS/qc/chipseeker_H3K27ac.pdf"
echo "  $RESULTS/qc/chipseeker_summary.txt"
