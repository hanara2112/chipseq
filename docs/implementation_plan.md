# ChIP-seq Assignment — Implementation Plan

**Assignment:** Histone modification ChIP-seq analysis (Mouse mm10)  
**Your marks:** H3K4me3 (promoter) + H3K27ac (enhancer) across 4 IL7 conditions  
**Conda env:** `chipseq-env` | **Project root:** `~/IIITH/S26/mb/project/`

> **Run every bash script with:** `conda activate chipseq-env && bash stepN_xxx.sh`

---

## Pipeline Overview

```
Raw FASTQs → FastQC → Trim Galore → Bowtie2 → SAMtools → Picard → MACS2 → ChIPseeker → IGV
  (Step 1)   (Step 1)  (Step 2)    (Step 3)   (Step 4)  (Step 4)  (Step 6)  (Step 8)   (Step 10)
```

---

## Step 1 — Raw Read Quality Assessment ✅ DONE

**Script:** `step1_fastqc.sh` | **Output:** `qc/fastqc_raw/multiqc_report.html`

### Assignment Questions & Answers

| # | Question | Answer |
|---|----------|--------|
| 1 | Average per-base sequence quality? | All 24 samples show Phred ≥ 30 across all 151 bp positions — **excellent quality** |
| 2 | Any positions where quality decreases significantly? | No — quality is uniformly high, no 3'-end drop observed |
| 3 | Evidence of adapter contamination? | **Yes** — all samples fail the Adapter Content module. Illumina TruSeq adapters detected in all 24 files |
| 4 | Why is QC important before downstream analysis? | Adapter sequences can misalign or reduce alignment rate; low-quality bases increase mismatches; QC ensures only clean data enters the pipeline, preventing false-positive peaks |

---

## Step 2 — Adapter Trimming

**Tool:** Trim Galore (wrapper for Cutadapt)  
**Script to create:** `step2_trim.sh`

### Commands
```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
DATA="$PROJECT/student_dataset"
OUT="$PROJECT/trimmed"

# Paired-end trimming for all 12 sample pairs
# Trim Galore auto-detects Illumina adapters
SAMPLES=(
  "IL7-H_H3k4me3_S294"       "IL7-H_H3k4me3_S294"
  "IL7-H_H3k27ac_S295"       "IL7-H_H3k27ac_S295"
  "IL7-H_igG_S286"           "IL7-H_igG_S286"
  "Il7_H400_H3k4me3_S296"    "Il7_H400_H3k4me3_S296"
  "IL7-H400_H3k27ac_S297"    "IL7-H400_H3k27ac_S297"
  "IL7-H400_IgG_S287"        "IL7-H400_IgG_S287"
  "IL7-Low_H3k4me3_S298"     "IL7-Low_H3k4me3_S298"
  "IL7-_Low_H3k27ac_S299"    "IL7-_Low_H3k27ac_S299"
  "IL7-Low_IgG_S288"         "IL7-Low_IgG_S288"
  "IL7-L-400_H3k4me3_S300"   "IL7-L-400_H3k4me3_S300"
  "IL7-L-400_H3k27ac_S301"   "IL7-L-400_H3k27ac_S301"
  "IL7-L-400IgG_S289"        "IL7-L-400IgG_S289"
)

for sample in IL7-H_H3k4me3_S294 IL7-H_H3k27ac_S295 IL7-H_igG_S286 \
              Il7_H400_H3k4me3_S296 IL7-H400_H3k27ac_S297 IL7-H400_IgG_S287 \
              IL7-Low_H3k4me3_S298 IL7-_Low_H3k27ac_S299 IL7-Low_IgG_S288 \
              IL7-L-400_H3k4me3_S300 IL7-L-400_H3k27ac_S301 IL7-L-400IgG_S289; do
  trim_galore \
    --paired \
    --cores 4 \
    --quality 20 \
    --fastqc \
    -o "$OUT" \
    "$DATA/${sample}_R1_1M.fastq.gz" \
    "$DATA/${sample}_R2_1M.fastq.gz"
done
```

### Expected Output
```
trimmed/
├── *_R1_1M_val_1.fq.gz     # 12 trimmed R1 files
├── *_R2_1M_val_2.fq.gz     # 12 trimmed R2 files
└── *_trimming_report.txt   # per-sample trim stats (keep these)
```

**Keep:** `*_val_1.fq.gz`, `*_val_2.fq.gz`, `*_trimming_report.txt`  
**Delete after alignment:** the trimmed `.fq.gz` files (large, no longer needed)

### Assignment Questions to Answer

| # | Question | Answer (fill after running) |
|---|----------|-----------------------------|
| 1 | Why must adapter sequences be removed before alignment? | Adapter bases are synthetic sequences not in the genome — they cause misalignments, reduce mapping rates, and introduce artefactual signals near read ends |
| 2 | What could happen if adapters are not removed? | Reads with adapters fail to align or align incorrectly → reduced coverage, false-negative peak calls, inflated background |
| 3 | What quality threshold was used? | Q20 (Phred 20 = 99% base call accuracy), standard for ChIP-seq |

---

## Step 3 — Alignment to Reference Genome (mm10)

**Tool:** Bowtie2  
**Script to create:** `step3_align.sh`  

> ⚠️ **First: Download mm10 index** — ~3.5 GB, do once:
> ```bash
> mkdir -p ~/IIITH/S26/mb/project/genome
> # Option A: Download pre-built Bowtie2 index from NCBI/iGenomes
> wget -P ~/IIITH/S26/mb/project/genome \
>   https://genome-idx.s3.amazonaws.com/bt/mm10.zip
> unzip ~/IIITH/S26/mb/project/genome/mm10.zip -d ~/IIITH/S26/mb/project/genome/
> ```

### Commands
```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
TRIMMED="$PROJECT/trimmed"
OUT="$PROJECT/aligned"
INDEX="$PROJECT/genome/mm10"   # path to Bowtie2 index prefix

for sample in IL7-H_H3k4me3_S294 IL7-H_H3k27ac_S295 IL7-H_igG_S286 \
              Il7_H400_H3k4me3_S296 IL7-H400_H3k27ac_S297 IL7-H400_IgG_S287 \
              IL7-Low_H3k4me3_S298 IL7-_Low_H3k27ac_S299 IL7-Low_IgG_S288 \
              IL7-L-400_H3k4me3_S300 IL7-L-400_H3k27ac_S301 IL7-L-400IgG_S289; do

  R1=$(ls "$TRIMMED/${sample}_R1_1M_val_1.fq.gz")
  R2=$(ls "$TRIMMED/${sample}_R2_1M_val_2.fq.gz")

  bowtie2 \
    -x "$INDEX" \
    -1 "$R1" \
    -2 "$R2" \
    --threads 8 \
    -X 2000 \
    --no-mixed \
    --no-discordant \
    2> "$OUT/${sample}_bowtie2.log" \
  | samtools sort -@ 4 -o "$OUT/${sample}.bam"

  samtools index "$OUT/${sample}.bam"
  echo "Done: $sample"
done
```

### Expected Output
```
aligned/
├── *.bam                    # sorted BAM files (keep)
├── *.bam.bai                # BAM indexes (keep)
└── *_bowtie2.log            # alignment stats (keep — fill table below)
```

### Assignment Questions to Answer

| # | Question | Answer (fill after running — check `*_bowtie2.log`) |
|---|----------|------------------------------------------------------|
| 1 | % reads aligned successfully? | **[TODO]** — check bowtie2 log, expect >70% for ChIP, ~90% for IgG |
| 2 | What does "properly paired reads" mean? | Both R1 and R2 align to same chromosome, correct orientation, within expected fragment size |
| 3 | Why specify `-X` (max fragment length)? | ChIP fragments are size-selected (~200bp); allowing very large inserts captures random ligation artefacts |
| 4 | Why might reads fail alignment? | Repetitive regions (multi-mapping), low-quality reads, reads derived from adapter dimers, or contamination from other organisms |

---

## Step 4 — BAM Processing & Duplicate Removal

**Tools:** SAMtools + Picard  
**Script to create:** `step4_process.sh`

### Commands
```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
ALIGNED="$PROJECT/aligned"

for sample in IL7-H_H3k4me3_S294 IL7-H_H3k27ac_S295 IL7-H_igG_S286 \
              Il7_H400_H3k4me3_S296 IL7-H400_H3k27ac_S297 IL7-H400_IgG_S287 \
              IL7-Low_H3k4me3_S298 IL7-_Low_H3k27ac_S299 IL7-Low_IgG_S288 \
              IL7-L-400_H3k4me3_S300 IL7-L-400_H3k27ac_S301 IL7-L-400IgG_S289; do

  BAM="$ALIGNED/${sample}.bam"

  # 1. Fix mate information (required before dedup)
  samtools fixmate -m "$BAM" "$ALIGNED/${sample}.fixmate.bam"
  samtools sort -o "$ALIGNED/${sample}.fixmate.sorted.bam" "$ALIGNED/${sample}.fixmate.bam"

  # 2. Mark & remove duplicates (Picard)
  picard MarkDuplicates \
    I="$ALIGNED/${sample}.fixmate.sorted.bam" \
    O="$ALIGNED/${sample}.dedup.bam" \
    M="$ALIGNED/${sample}.dup_metrics.txt" \
    REMOVE_DUPLICATES=true

  samtools index "$ALIGNED/${sample}.dedup.bam"

  # 3. MAPQ filter >= 20
  samtools view -b -q 20 \
    "$ALIGNED/${sample}.dedup.bam" \
    > "$ALIGNED/${sample}.filtered.bam"
  samtools index "$ALIGNED/${sample}.filtered.bam"

  # 4. Cleanup intermediates
  rm "$ALIGNED/${sample}.fixmate.bam" \
     "$ALIGNED/${sample}.fixmate.sorted.bam" \
     "$ALIGNED/${sample}.dedup.bam"* \

  echo "Done: $sample → ${sample}.filtered.bam"
done
```

### Expected Output
```
aligned/
├── *.filtered.bam           # final analysis-ready BAMs (keep)
├── *.filtered.bam.bai       # indexes (keep)
└── *.dup_metrics.txt        # duplication stats (keep — fill table)
```

### Assignment Questions to Answer

| # | Question | Answer |
|---|----------|--------|
| 1 | What are PCR duplicates and how do they arise? | Identical read pairs from the same original DNA fragment amplified multiple times during PCR library prep; identified by identical start coordinates |
| 2 | Why do duplicates cause false signals in ChIP-seq? | They artificially inflate read depth at specific positions → appear as sharp "peaks" that are PCR artefacts, not true enrichment |
| 3 | Purpose of fixing mate information before dedup? | Paired-end dedup requires mate TLEN (insert size) info in the BAM; `fixmate` populates these fields so Picard can correctly identify duplicate pairs |

---

## Step 5 — Filtering Low-Quality Alignments

> ⚠️ Covered in Step 4 script (`samtools view -q 20`). Questions still need answers:

### Assignment Questions to Answer

| # | Question | Answer |
|---|----------|--------|
| 1 | What does MAPQ score represent? | Phred-scaled probability that the read is mapped to the wrong location; MAPQ=20 → 99% confidence in alignment position |
| 2 | Why is MAPQ ≥ 20 commonly used? | Removes multi-mapping reads (MAPQ=0 or 1) and low-confidence alignments; standard threshold balancing sensitivity vs specificity |
| 3 | What genomic regions produce low MAPQ? | Repetitive elements (SINEs, LINEs, satellite repeats), segmental duplications, centromeric regions — reads can't be uniquely placed |

---

## Step 6 — Peak Calling with MACS2

**Tool:** MACS2  
**Script to create:** `step6_peaks.sh`

### Commands
```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/IIITH/S26/mb/project"
ALIGNED="$PROJECT/aligned"
PEAKS="$PROJECT/peaks"

# === H3K4me3 peak calling (narrow peaks — sharp promoter mark) ===
declare -A H3K4ME3_SAMPLES=(
  ["IL7-H_H3k4me3_S294"]="IL7-H_igG_S286"
  ["Il7_H400_H3k4me3_S296"]="IL7-H400_IgG_S287"
  ["IL7-Low_H3k4me3_S298"]="IL7-Low_IgG_S288"
  ["IL7-L-400_H3k4me3_S300"]="IL7-L-400IgG_S289"
)

for chip in "${!H3K4ME3_SAMPLES[@]}"; do
  ctrl="${H3K4ME3_SAMPLES[$chip]}"
  macs2 callpeak \
    -t "$ALIGNED/${chip}.filtered.bam" \
    -c "$ALIGNED/${ctrl}.filtered.bam" \
    -f BAMPE \
    -g mm \
    -n "$chip" \
    --outdir "$PEAKS/H3K4me3" \
    -q 0.05 \
    2> "$PEAKS/H3K4me3/${chip}_macs2.log"
done

# === H3K27ac peak calling (broad peaks — diffuse enhancer mark) ===
declare -A H3K27AC_SAMPLES=(
  ["IL7-H_H3k27ac_S295"]="IL7-H_igG_S286"
  ["IL7-H400_H3k27ac_S297"]="IL7-H400_IgG_S287"
  ["IL7-_Low_H3k27ac_S299"]="IL7-Low_IgG_S288"
  ["IL7-L-400_H3k27ac_S301"]="IL7-L-400IgG_S289"
)

for chip in "${!H3K27AC_SAMPLES[@]}"; do
  ctrl="${H3K27AC_SAMPLES[$chip]}"
  macs2 callpeak \
    -t "$ALIGNED/${chip}.filtered.bam" \
    -c "$ALIGNED/${ctrl}.filtered.bam" \
    -f BAMPE \
    -g mm \
    -n "$chip" \
    --outdir "$PEAKS/H3K27ac" \
    --broad \
    -q 0.05 \
    2> "$PEAKS/H3K27ac/${chip}_macs2.log"
done

# === IgG self-peak calling (should give very few/no peaks) ===
for igG in IL7-H_igG_S286 IL7-H400_IgG_S287 IL7-Low_IgG_S288 IL7-L-400IgG_S289; do
  macs2 callpeak \
    -t "$ALIGNED/${igG}.filtered.bam" \
    -f BAMPE \
    -g mm \
    -n "$igG" \
    --outdir "$PEAKS/IgG" \
    -q 0.05 \
    2> "$PEAKS/IgG/${igG}_macs2.log"
done
```

### Expected Output
```
peaks/
├── H3K4me3/
│   ├── *_peaks.narrowPeak    # peak coordinates (keep ⭐)
│   ├── *_summits.bed         # peak summit positions (keep ⭐)
│   └── *_peaks.xls           # peak stats table (keep ⭐)
├── H3K27ac/
│   ├── *_peaks.broadPeak     # broad peak coordinates (keep ⭐)
│   └── *_peaks.xls
└── IgG/
    └── *_peaks.narrowPeak    # should be very few
```

### Assignment Questions to Answer

| # | Question | Answer |
|---|----------|--------|
| 1 | What does a "peak" represent in ChIP-seq? | A genomic region where ChIP reads are significantly enriched above background — indicates where the histone mark or protein is bound |
| 2 | Why is a control sample required? | IgG background defines local noise levels; MACS2 uses it to calculate fold-enrichment and p-values, preventing PCR/sequencing biases from being called as peaks |
| 3 | How does MACS2 distinguish real enrichment from noise? | Models fragment size, estimates local background from control, computes Poisson p-values, and applies q-value (FDR) threshold |

---

## Step 7 — Comparison: ChIP vs IgG Peak Counts

### Peak Count Summary Table *(fill after Step 6)*

| Sample | Condition | Mark | # Peaks | # IgG Peaks | Enrichment |
|--------|-----------|------|---------|-------------|------------|
| S294 | IL7-High | H3K4me3 | **[TODO]** | **[TODO]** | **[TODO]** |
| S296 | IL7-H400 | H3K4me3 | **[TODO]** | **[TODO]** | **[TODO]** |
| S298 | IL7-Low | H3K4me3 | **[TODO]** | **[TODO]** | **[TODO]** |
| S300 | IL7-L-400 | H3K4me3 | **[TODO]** | **[TODO]** | **[TODO]** |
| S295 | IL7-High | H3K27ac | **[TODO]** | **[TODO]** | **[TODO]** |
| S297 | IL7-H400 | H3K27ac | **[TODO]** | **[TODO]** | **[TODO]** |
| S299 | IL7-Low | H3K27ac | **[TODO]** | **[TODO]** | **[TODO]** |
| S301 | IL7-L-400 | H3K27ac | **[TODO]** | **[TODO]** | **[TODO]** |

> Quick count command: `wc -l peaks/H3K4me3/*_peaks.narrowPeak`

### Assignment Questions to Answer

| # | Question | Answer (fill after Step 6) |
|---|----------|---------------------------|
| 1 | How many peaks in your histone mark dataset? | **[TODO]** |
| 2 | How many peaks in IgG control? | **[TODO]** — expect very few (<100) |
| 3 | Why should IgG have significantly fewer peaks? | IgG is non-specific antibody — it samples genomic background uniformly without enriching any specific locus; any IgG "peaks" are sequencing/PCR artefacts |

---

## Step 8 — Genomic Distribution of Peaks (ChIPseeker)

**Tool:** ChIPseeker (R package)  
**Script to create:** `step8_chipseeker.R`

### Install (once)
```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                        "org.Mm.eg.db", "clusterProfiler"))
```

### Commands
```r
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Load peaks
h3k4me3_peaks <- list(
  High    = readPeakFile("peaks/H3K4me3/IL7-H_H3k4me3_S294_peaks.narrowPeak"),
  H400    = readPeakFile("peaks/H3K4me3/Il7_H400_H3k4me3_S296_peaks.narrowPeak"),
  Low     = readPeakFile("peaks/H3K4me3/IL7-Low_H3k4me3_S298_peaks.narrowPeak"),
  L400    = readPeakFile("peaks/H3K4me3/IL7-L-400_H3k4me3_S300_peaks.narrowPeak")
)

h3k27ac_peaks <- list(
  High    = readPeakFile("peaks/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak"),
  H400    = readPeakFile("peaks/H3K27ac/IL7-H400_H3k27ac_S297_peaks.broadPeak"),
  Low     = readPeakFile("peaks/H3K27ac/IL7-_Low_H3k27ac_S299_peaks.broadPeak"),
  L400    = readPeakFile("peaks/H3K27ac/IL7-L-400_H3k27ac_S301_peaks.broadPeak")
)

# Annotate
anno_k4 <- lapply(h3k4me3_peaks, annotatePeak, TxDb=txdb, 
                   tssRegion=c(-3000,3000), annoDb="org.Mm.eg.db")
anno_k27 <- lapply(h3k27ac_peaks, annotatePeak, TxDb=txdb,
                    tssRegion=c(-3000,3000), annoDb="org.Mm.eg.db")

# Plots
pdf("qc/chipseeker_annotation.pdf", width=12, height=8)
plotAnnoBar(anno_k4, title="H3K4me3 — Genomic Distribution")
plotAnnoBar(anno_k27, title="H3K27ac — Genomic Distribution")
plotDistToTSS(anno_k4, title="H3K4me3 — Distance to TSS")
plotDistToTSS(anno_k27, title="H3K27ac — Distance to TSS")
dev.off()
```

### Assignment Questions to Answer *(fill after running ChIPseeker)*

| # | Question | Answer |
|---|----------|--------|
| 1 | % peaks in promoter regions? | **[TODO]** — H3K4me3 expect ~60–80%; H3K27ac expect ~20–40% |
| 2 | % peaks in intergenic regions? | **[TODO]** — H3K27ac expect higher % (enhancers are distal) |
| 3 | Why do enhancer marks appear farther from promoters? | Enhancers are cis-regulatory elements that can be 10–100s of kb from their target gene; H3K27ac marks active enhancers which are predominantly in intronic/intergenic space |

---

## Step 9 — Biological Interpretation

### Assignment Questions to Answer

| # | Question | Answer |
|---|----------|--------|
| 1 | How does H3K4me3 distribution differ from H3K27ac? | H3K4me3 is tightly focused at TSS (narrow peaks, >60% promoter); H3K27ac is broader and more distal, marking both promoters and active enhancers |
| 2 | Why is H3K4me3 strongly associated with TSS? | H3K4me3 is deposited by the COMPASS/SET1 complex recruited by the Pol II pre-initiation complex; it recruits TFIID and promotes transcription initiation |
| 3 | What regulatory role do H3K27ac enhancers play? | H3K27ac enhancers loop to target promoters via cohesin/CTCF, recruit co-activators (p300/CBP), and drive cell-type-specific gene expression programs |

---

## Step 10 — Visualization in IGV

**Tool:** IGV (Integrative Genomics Viewer) + deepTools (bigWig generation)

### Generate bigWig files
```bash
# Add to step6 or create step10_bigwig.sh
for sample in ...; do
  bamCoverage \
    -b aligned/${sample}.filtered.bam \
    -o qc/${sample}.bw \
    --normalizeUsing RPKM \
    --binSize 10 \
    --extendReads \
    -p 8
done
```

### In IGV
1. Load genome: **mm10**
2. Load `.bw` files for ChIP + IgG
3. Navigate to a known active gene (e.g., *Actb*, *Gapdh*)
4. Screenshot a representative peak for submission

### Assignment Questions to Answer

| # | Question | Answer |
|---|----------|--------|
| 1 | What does peak height represent in genome browser? | Read depth / coverage at that position — higher = more ChIP-enriched reads, stronger histone modification signal |
| 2 | Why are H3K4me3 peaks sharper than H3K27ac? | H3K4me3 is precisely deposited at TSS (active promoter, ~200bp nucleosome), while H3K27ac marks broad enhancer domains spanning several kb |
| 3 | How does IgG signal differ from ChIP in IGV? | IgG shows flat, uniform low-level signal (background); ChIP shows distinct sharp spikes above background at biologically meaningful loci |

---

## Step 11 — Biological Conclusions

### Assignment Questions to Answer

| # | Question | Answer |
|---|----------|--------|
| 1 | What conclusions from promoter-associated histone marks? | H3K4me3 at promoters indicates actively transcribed genes; differences between IL7-High and IL7-Low conditions suggest IL7 signaling alters the transcriptional landscape at promoters |
| 2 | How do enhancer marks contribute to gene regulation? | H3K27ac-marked enhancers coordinate distal regulatory inputs, amplify promoter activity, and enable cell-type/stimulus-specific gene expression patterns |
| 3 | Why are chromatin modifications considered epigenetic regulators? | They are heritable through cell division (written by specific enzymes, read by effector proteins, erased by demethylases/deacetylases) without changing DNA sequence — they encode cellular memory of gene expression state |

---

## Submission Checklist

- [ ] **Peak count summary table** (Step 7 table filled in)
- [ ] **IGV genome browser screenshot** showing one representative peak with ChIP + IgG tracks
- [ ] **1–2 page written explanation** covering:
  - QC summary (Step 1)
  - Trimming rationale (Step 2)
  - Alignment stats (Step 3)
  - Peak counts ChIP vs IgG (Step 7)
  - Genomic distribution (Step 8)
  - Biological interpretation (Steps 9, 11)

---

## File Cleanup Guide

| After step | Delete | Keep |
|------------|--------|------|
| Step 2 | — | `*_val_1/2.fq.gz`, `*_trimming_report.txt` |
| Step 3 | Raw `.bam` | `*.sorted.bam`, `*.bai`, `*_bowtie2.log` |
| Step 4 | `*.fixmate.bam`, `*.dedup.bam` | `*.filtered.bam`, `*.bai`, `*.dup_metrics.txt` |
| Step 6 | `*_model.r`, `*_model.pdf` | `*.narrowPeak`, `*.broadPeak`, `*_summits.bed`, `*.xls` |
| Step 10 | — | `*.bw` bigWig files |
| After trimming | Trimmed FASTQs | Final filtered BAMs + peaks |
