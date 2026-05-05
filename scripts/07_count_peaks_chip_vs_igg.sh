#!/usr/bin/env bash
# Step 7 — Count peaks (ChIP vs IgG) and write summary table
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/data"
RESULTS="$PROJECT/results"
PEAKS="$RESULTS/peaks"
QC="$RESULTS/qc"
mkdir -p "$QC"

CHIP_PEAKS="$PEAKS/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak"
IGG_PEAKS="$PEAKS/IgG/IL7-H_igG_S286_peaks.narrowPeak"

# count peaks (handle missing IgG file: MACS2 may produce no narrowPeak if 0 peaks)
N_CHIP=$( [ -s "$CHIP_PEAKS" ] && wc -l < "$CHIP_PEAKS" || echo 0 )
N_IGG=$(  [ -s "$IGG_PEAKS"  ] && wc -l < "$IGG_PEAKS"  || echo 0 )

# enrichment ratio (avoid div-by-zero)
if [ "$N_IGG" -eq 0 ]; then
  RATIO="∞ (no IgG peaks)"
else
  RATIO=$(awk -v c="$N_CHIP" -v i="$N_IGG" 'BEGIN{printf "%.1fx", c/i}')
fi

OUT="$RESULTS/peak_count_summary.tsv"
{
  printf "Sample\tCondition\tMark\tPeak_File\tN_Peaks\n"
  printf "S295\tIL7-High\tH3K27ac\t%s\t%s\n" "$(basename "$CHIP_PEAKS")" "$N_CHIP"
  printf "S286\tIL7-High\tIgG\t%s\t%s\n"     "$(basename "$IGG_PEAKS")"  "$N_IGG"
} > "$OUT"

echo "=== Peak Count Summary ==="
column -t -s $'\t' "$OUT"
echo ""
echo "ChIP / IgG enrichment: $RATIO"
echo ""
echo "Saved: $OUT"
