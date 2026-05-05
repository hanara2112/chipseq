#!/usr/bin/env bash
# Step 5 — Filter low-MAPQ + keep only properly-paired reads
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
ALIGNED="$PROJECT/aligned"
QC="$PROJECT/qc"

SAMPLES=(IL7-H_H3k27ac_S295 IL7-H_igG_S286)
MAPQ=20

echo "=== Step 5: MAPQ>=${MAPQ} + properly-paired filter ==="

for sample in "${SAMPLES[@]}"; do
  IN="$ALIGNED/${sample}_deduplicated.bam"
  OUT="$ALIGNED/${sample}_filtered_MAPQ20.bam"
  STATS="$QC/${sample}_filter_stats.txt"

  echo ">>> Filtering $sample"

  # -b BAM out, -q 20 MAPQ filter, -f 2 properly paired, -F 1804 drops unmapped/secondary/dup/QC-fail
  samtools view -b -@ 4 -q "$MAPQ" -f 2 -F 1804 "$IN" -o "$OUT"
  samtools index "$OUT"

  {
    echo "Sample:           $sample"
    echo "Input BAM:        $IN"
    echo "Output BAM:       $OUT"
    echo "MAPQ threshold:   >=$MAPQ"
    echo "Flag require (-f): 2  (properly paired)"
    echo "Flag exclude (-F): 1804  (unmapped/mate-unmapped/secondary/QCfail/dup)"
    echo "--- read counts ---"
    echo "before: $(samtools view -c "$IN")"
    echo "after:  $(samtools view -c "$OUT")"
  } > "$STATS"

  echo "  done: $OUT"
  echo "  stats: $STATS"
done

echo "=== Step 5 complete ==="
