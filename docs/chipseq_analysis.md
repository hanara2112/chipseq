# ChIP-seq Analysis Pipeline

**Project:** IL7 High/Low — H3K4me3 & H3K27ac ChIP-seq  
**Conda env:** `chipseq-env`  
**Project root:** `~/IIITH/S26/mb/project/`

---

## How to Run

Each step has its own bash script. Run them inside `chipseq-env`:
```bash
conda activate chipseq-env
bash step1_fastqc.sh
# bash step2_trim.sh   (coming next)
# bash step3_align.sh
# bash step4_peaks.sh
```

---

## Directory Structure

```
~/IIITH/S26/mb/project/
├── student_dataset/          # Raw FASTQ files (Only processing assigned S295 & S286)
├── qc/
│   └── fastqc_raw/
│       ├── multiqc_report.html      # ⭐ Main QC report (open this)
│       └── multiqc_report_data/     # TSV/JSON stats
├── trimmed/                  # Trimmed reads (Step 2)
├── aligned/                  # BAM files (Step 3)
├── peaks/                    # Called peaks (Step 4)
├── step1_fastqc.sh           # Bash script — Step 1
└── chipseq_analysis.md       # This file
```

---

## Samples

| Sample ID | Condition | Mark | R1 | R2 |
|-----------|-----------|------|----|----|
| S295 | IL7-High | H3K27ac | IL7-H_H3k27ac_S295_R1_1M.fastq.gz | IL7-H_H3k27ac_S295_R2_1M.fastq.gz |
| S286 | IL7-High | IgG ctrl | IL7-H_igG_S286_R1_1M.fastq.gz | IL7-H_igG_S286_R2_1M.fastq.gz |

---

## Step 1: Raw Read Quality Assessment

### Script

```bash
bash step1_fastqc.sh
```

See [`step1_fastqc.sh`](step1_fastqc.sh) for the full annotated script.

### What Was Kept (after cleanup)

```
qc/fastqc_raw/
├── multiqc_report.html      # ⭐ Single combined interactive report
└── multiqc_report_data/     # TSV/JSON stats for scripting
```

**Deleted** (redundant once MultiQC is run):
- 4× `*_fastqc.html` — individual per-sample reports, subsumed by MultiQC
- 4× `*_fastqc.zip` — intermediate data consumed by MultiQC
- `fastqc_run.log` — transient run log

### Results — General Statistics

> Tool: FastQC v0.12.1 | MultiQC v1.34 | All files subsampled to 1M reads, read length 151 bp

| Sample | Reads | %Dups | %GC | Seq Quality | Per-base Content | GC Content | Adapter Content |
|--------|-------|-------|-----|-------------|-----------------|------------|----------------|
| IL7-H_H3k27ac_S295_R1 | 1M | 5.8% | 51% | ✅ pass | ✅ pass | ⚠️ warn | ❌ fail |
| IL7-H_H3k27ac_S295_R2 | 1M | 5.4% | 51% | ✅ pass | ✅ pass | ⚠️ warn | ❌ fail |
| IL7-H_igG_S286_R1 | 1M | 9.0% | 43% | ✅ pass | ✅ pass | ❌ fail | ❌ fail |
| IL7-H_igG_S286_R2 | 1M | 8.6% | 42% | ✅ pass | ⚠️ warn | ❌ fail | ❌ fail |


### Interpretation

#### ✅ Per Base Sequence Quality
All assigned samples **pass**. Phred scores are consistently high across all 151 bp positions. No 3'-end quality drop.

#### ⚠️ Per Sequence GC Content — flagged but EXPECTED
- Most samples show a `warn` or `fail` on GC content.
- This is **normal for ChIP-seq** — the immunoprecipitation enriches for specific genomic regions that are inherently GC-biased (especially H3K4me3 at CpG-rich promoters, GC ~53%).
- IgG controls have lower GC (41–43%) as expected — they sample the genome more uniformly.
- **Not a problem. No action needed.**

#### ❌ Adapter Content — flagged, REQUIRES trimming
- All samples fail adapter content.
- Illumina TruSeq adapters detected (standard for paired-end libraries).
- **Action: Run Trim Galore in Step 2** (already installed in `chipseq-env`).

#### Sequence Duplication
- Low duplication: 5–10% for ChIP samples (**excellent library complexity**).
- Slightly higher for IgG controls (~9–14%) — normal.
- Well within acceptable range. No samples need to be excluded.

### Summary / Decision
> All assigned samples pass base quality QC. Adapter contamination is present in all samples (expected for Illumina paired-end). Proceeding to **Step 2: Trim Galore** to remove adapters before alignment.

---

## Step 2: Adapter Trimming (Trim Galore)

### Script

```bash
bash step2_trim.sh
```

### Expected Output
```
trimmed/
├── *_R1_1M_val_1.fq.gz     # 2 trimmed R1 files (discard after alignment)
├── *_R2_1M_val_2.fq.gz     # 2 trimmed R2 files (discard after alignment)
└── *_trimming_report.txt   # per-sample trim stats (keep these)
```

### Assignment Questions & Answers

| # | Question | Answer |
|---|----------|--------|
| 1 | Why must adapter sequences be removed before alignment? | Adapter bases are synthetic sequences not present in the reference genome. If left attached to the reads, they cause misalignments, reduce the overall mapping rate, and can introduce artefactual signals near read ends. |
| 2 | What could happen if adapters are not removed? | Reads with adapters will either fail to align or align incorrectly. This leads to reduced coverage, false-negative peak calls, and an inflated background noise level. |
| 3 | What quality threshold was used for trimming? | Q20 (Phred score of 20, which equals 99% base call accuracy). This is the standard threshold used by Trim Galore and ensures low-quality ends are removed along with adapters. |

---

## Step 3: Alignment (Bowtie2)

### Script

```bash
bash step3_align.sh
```

### Expected Output
```
aligned/
├── *.bam                    # sorted BAM files (keep)
├── *.bam.bai                # BAM indexes (keep)
└── *_bowtie2.log            # alignment stats (keep — fill table below)
```

### Assignment Questions & Answers

| # | Question | Answer |
|---|----------|--------|
| 1 | What percentage of reads aligned successfully to the genome? | *(Check your `*_bowtie2.log` files once alignment is done. Usually >70% for good ChIP, ~90% for IgG)* |
| 2 | What does "properly paired reads" mean? | It means both the forward (R1) and reverse (R2) reads aligned to the same chromosome, in the correct orientation facing each other, and within the expected fragment length. |
| 3 | Why is the maximum fragment length (`-X`) specified during alignment? | ChIP-seq libraries are size-selected during preparation (typically 200–500 bp). Restricting the maximum insert size prevents Bowtie2 from incorrectly pairing reads that are thousands of bases apart due to random genomic ligations or structural variations. |
| 4 | What reasons might cause reads to fail alignment? | Reads can fail to align due to low sequencing quality, adapter contamination that wasn't fully trimmed, sequencing errors, the reads coming from unmapped highly repetitive regions, or contamination from another organism. |

---

## Step 4: BAM Processing & PCR Duplicate Removal

### Script
```bash
bash 04_remove_pcr_duplicates.sh
```

Pipeline: name-sort → `samtools fixmate -m` → coord-sort → `samtools markdup -r` (dedup tool **B**, see note below).

> **Note on tool choice:** original plan used Picard `MarkDuplicates`, but it isn't installed in `chipseq-env`.
> Switched to `samtools markdup` — already in env, equivalent for paired-end PCR dedup, one less dependency.

### Output

```
aligned/
├── IL7-H_H3k27ac_S295_deduplicated.bam (+ .bai)   # 130 MB
└── IL7-H_igG_S286_deduplicated.bam (+ .bai)        # 143 MB
qc/
├── IL7-H_H3k27ac_S295_dup_metrics.txt
└── IL7-H_igG_S286_dup_metrics.txt
```

### Results — duplicate metrics

| Sample | Reads in | Reads out | Duplicate pairs | Dup rate | Library size (est.) |
|--------|----------|-----------|-----------------|----------|---------------------|
| S295 (H3K27ac) | 1,997,162 | 1,834,564 | 162,598 | **8.14%** | 5,304,792 |
| S286 (IgG)     | 1,996,996 | 1,802,666 | 194,330 | **9.74%** | 4,339,762 |

Both within healthy range (<15%). IgG has slightly higher duplication, which is expected — non-specific antibody enriches less, so PCR amplifies a smaller pool of unique fragments more times.

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | What are PCR duplicates and how do they arise? | Identical read pairs derived from the same original DNA fragment that was amplified multiple times during PCR library prep. They are identified by having identical mapped 5' coordinates for both mates. |
| 2 | Why do duplicates cause false signals in ChIP-seq? | They artificially inflate read depth at specific positions, producing sharp pileups that look like enriched peaks but actually reflect PCR amplification bias rather than true antibody binding. |
| 3 | Purpose of fixing mate information before dedup? | Paired-end deduplication needs each read to know its mate's position and the insert size. `samtools fixmate -m` fills in the MC (mate cigar) and ms (mate score) tags so `markdup` can correctly identify which read of a duplicate pair to keep. |

---

## Step 5: Filtering Low-Quality Alignments

### Script
```bash
bash 05_filter_low_mapq_alignments.sh
```

Filter applied: `samtools view -q 20 -f 2 -F 1804`
- `-q 20`: MAPQ ≥ 20 (drops multi-mappers and low-confidence alignments)
- `-f 2`: keep only properly-paired reads
- `-F 1804`: exclude unmapped, mate-unmapped, secondary, QC-fail, duplicate

### Output
```
aligned/
├── IL7-H_H3k27ac_S295_filtered_MAPQ20.bam (+ .bai)   # 116 MB — analysis-ready
└── IL7-H_igG_S286_filtered_MAPQ20.bam (+ .bai)        # 117 MB — analysis-ready
qc/
├── IL7-H_H3k27ac_S295_filter_stats.txt
└── IL7-H_igG_S286_filter_stats.txt
```

### Results — read retention

| Sample | After dedup | After MAPQ20 + proper-pair | Removed | % retained |
|--------|-------------|---------------------------|---------|------------|
| S295 (H3K27ac) | 1,834,564 | **1,660,742** | 173,822 | 90.5% |
| S286 (IgG)     | 1,802,666 | **1,461,448** | 341,218 | 81.1% |

IgG drops more reads at this step — expected. IgG samples the genome uniformly (including repetitive elements), so a larger fraction of its reads land in regions where unique placement is impossible (MAPQ 0–1) and get filtered out.

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | What does MAPQ score represent? | A Phred-scaled estimate of the probability that the read is mapped to the wrong location. MAPQ = 20 corresponds to a 1% chance of mis-mapping (99% confidence). Higher MAPQ = more uniquely mappable. |
| 2 | Why is MAPQ ≥ 20 commonly used? | It is a standard, well-balanced threshold: removes all multi-mappers (MAPQ 0–1) and low-confidence reads while retaining the bulk of uniquely-aligned, high-quality data. Stricter thresholds (e.g. ≥30) start removing real signal; looser ones let noise through. |
| 3 | Which genomic regions produce low MAPQ? | Repetitive elements (SINEs, LINEs, simple repeats, satellite DNA), segmental duplications, centromeric and pericentromeric regions, and recently duplicated paralogous gene families — anywhere the reference contains nearly-identical sequence at multiple loci, so the aligner cannot uniquely place the read. |

---

## Step 6: Peak Calling (MACS2)

### Script
```bash
bash 06_call_peaks_macs2_h3k27ac_vs_igg.sh
```

Two MACS2 runs:
- **6a** — H3K27ac (S295) vs IgG (S286) control, `--broad` mode (q ≤ 0.05 strong / 0.10 weak), `-f BAMPE`, `-g mm`
- **6b** — IgG (S286) self-call (no control) — establishes whether the IgG library has any locus-specific signal of its own

### Run summary (from MACS2 logs)

| | H3K27ac (6a) | IgG (6b) |
|---|---|---|
| Fragments after PE collapse | 830,371 | 730,724 |
| Mean fragment size | 250.6 bp | 220.7 bp |
| Effective genome size | 1.87e9 (mm10) | 1.87e9 |
| Mode | `--broad` | narrow |
| Control | IgG (S286) | none |
| q-value cutoff | 0.05 / 0.10 | 0.05 |

### Output

```
peaks/H3K27ac/
├── IL7-H_H3k27ac_S295_peaks.broadPeak    ⭐ 18,681 peaks
├── IL7-H_H3k27ac_S295_peaks.gappedPeak       (bed12 with sub-peak structure)
├── IL7-H_H3k27ac_S295_peaks.xls              (peak stats table)
└── IL7-H_H3k27ac_S295_macs2.log

peaks/IgG/
├── IL7-H_igG_S286_peaks.narrowPeak       ⭐ 0 peaks (empty)
├── IL7-H_igG_S286_summits.bed                (empty)
├── IL7-H_igG_S286_peaks.xls
└── IL7-H_igG_S286_macs2.log
```

### H3K27ac peak QC

| Metric | Value |
|---|---|
| Number of peaks | **18,681** |
| Median peak length | 1,218 bp |
| Mean peak length | 1,686 bp |
| Min / max peak length | 250 / 27,455 bp |
| Mean fold-enrichment over IgG | **4.54×** |
| Max fold-enrichment | 17.12× |

Peak length distribution is consistent with H3K27ac biology — broad domains spanning multiple nucleosomes (1–2 kb median) at active enhancers, with occasional very-broad super-enhancer-like regions (max ~27 kb).

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | What does a "peak" represent in ChIP-seq? | A contiguous genomic interval where the density of immunoprecipitated reads is statistically higher than the local background — i.e. a region the antibody (here, anti-H3K27ac) repeatedly bound, indicating that the histone modification is enriched there. |
| 2 | Why is a control sample required? | The IgG (or input) control measures sample-specific background: chromatin accessibility bias, sonication bias, GC bias, sequencing artefacts, and non-specific antibody pulldown. Without it, MACS2 would call peaks anywhere reads are simply more abundant — many of which would be technical noise rather than true H3K27ac binding. |
| 3 | How does MACS2 distinguish real enrichment from noise? | MACS2 (a) estimates fragment size from read pileups, (b) models the local read distribution under a dynamic Poisson model whose lambda parameter is the *maximum* of several local estimates (1 kb / 5 kb / 10 kb windows from the control), (c) computes a p-value at every position, then (d) applies Benjamini–Hochberg FDR (q-value) to get an FDR-controlled peak list. In `--broad` mode it merges adjacent significant regions allowing weaker (q ≤ 0.10) bridges between strong (q ≤ 0.05) cores to capture diffuse marks. |

---

## Step 7: ChIP vs IgG Peak Comparison

### Script
```bash
bash 07_count_peaks_chip_vs_igg.sh
```

### Peak count summary

| Sample | Condition | Mark | Peak file | # Peaks |
|---|---|---|---|---|
| **S295** | IL7-High | H3K27ac (broad) | `IL7-H_H3k27ac_S295_peaks.broadPeak` | **18,681** |
| **S286** | IL7-High | IgG (control)   | `IL7-H_igG_S286_peaks.narrowPeak`     | **0** |

**ChIP / IgG enrichment ratio: ∞** — the IgG self-call returned zero significant peaks, which is the ideal outcome.

### Interpretation

The IgG library produced no peaks even when called with no control (i.e. against MACS2's own genome-wide lambda model). This is exactly what we want: it means the IgG signal is approximately uniformly distributed across the genome, so there is no locus-specific binding to worry about. By extension, **all 18,681 H3K27ac peaks called against IgG can be confidently attributed to genuine H3K27ac enrichment**, not antibody non-specificity, sonication bias, or PCR artefacts.

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | How many peaks in your histone-mark dataset? | **18,681** broad H3K27ac peaks (q ≤ 0.05, mm10, IgG-controlled). |
| 2 | How many peaks in IgG control? | **0** — the IgG sample has no statistically significant enriched regions even when called against MACS2's null model. |
| 3 | Why should IgG datasets contain significantly fewer peaks? | IgG is a non-specific isotype-control antibody. It does not target any chromatin feature, so the immunoprecipitated DNA reflects bulk genome composition (uniformly distributed reads) rather than enrichment at biologically meaningful loci. Any "peaks" that *do* appear in an IgG sample are usually technical: PCR duplicates that survived dedup, mappability hotspots, or repetitive elements with inflated read counts. A well-prepared IgG library should produce 0 — or at most a handful — of peaks, which is exactly what we observe here. |

---

## Step 8: Genomic Distribution (ChIPseeker)

### Script
```bash
bash 08_run_chipseeker.sh           # wrapper
# internally calls: Rscript 08_chipseeker_annotate.R
```

Annotation reference: `TxDb.Mmusculus.UCSC.mm10.knownGene` + `org.Mm.eg.db`
TSS region: −3,000 / +3,000 bp.

### Output

```
peaks/H3K27ac/H3K27ac_annotated.csv     # 18,681 rows × 23 cols (per-peak annotation incl. nearest gene SYMBOL)
qc/chipseeker_H3K27ac.pdf               # pie + bar + distance-to-TSS plots
qc/chipseeker_summary.txt               # feature-class % table
```

### Genomic feature distribution

| Feature | % of peaks |
|---|---|
| Promoter (≤1 kb)        | 50.47 |
| Promoter (1–2 kb)       | 2.28 |
| Promoter (2–3 kb)       | 2.59 |
| 5′ UTR                  | 0.44 |
| 3′ UTR                  | 1.85 |
| 1st Exon                | 2.30 |
| Other Exon              | 3.13 |
| 1st Intron              | 7.24 |
| Other Intron            | 10.84 |
| Downstream (≤300 bp)    | 0.05 |
| Distal Intergenic       | 18.82 |

**Aggregated:**

| Region | % |
|---|---|
| **Promoter** (all sub-bins) | **55.34** |
| **Intergenic** (distal)     | **18.82** |
| **Intron** (1st + other)    | **18.08** |
| **Exon** (1st + other + UTRs) | **8.05** |

### Interpretation

H3K27ac marks **active** chromatin — both active promoters *and* active enhancers — so a substantial promoter share is biologically expected. The 55% promoter / 19% distal-intergenic / 18% intron split here indicates that in IL7-High T-cells the H3K27ac signal is dominated by transcriptionally active promoters, with a sizeable secondary pool at distal/intronic enhancer elements (~37% combined).

Two technical contributors likely amplify the promoter fraction beyond a "typical" H3K27ac profile (~20–40% promoter):
1. Broad peaks straddle TSSs — when a 1–2 kb broadPeak overlaps a TSS even partially, ChIPseeker classifies the *whole* peak as "Promoter".
2. The 1M-read subsample preferentially recovers the strongest, deepest-coverage peaks; the strongest H3K27ac peaks are usually at active promoters (where H3K4me3 + H3K27ac co-occur).

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | % peaks in promoter regions? | **55.34 %** (combining ≤1 kb + 1–2 kb + 2–3 kb sub-bins). The largest sub-bin is Promoter ≤1 kb (50.5%), reflecting precise nucleosomal positioning of H3K27ac at active TSSs. |
| 2 | % peaks in intergenic regions? | **18.82 %** (distal intergenic). A further 18.08% lie in introns, which functionally also host many enhancers; combined non-genic-promoter signal (intron + intergenic) is ~37%. |
| 3 | Why might enhancer marks appear farther from gene promoters? | Enhancers are *cis*-regulatory elements that physically loop to their target promoter via cohesin/CTCF-mediated 3D chromatin architecture, so they are not constrained to lie next to the gene they regulate. They typically live tens of kb away — most often in introns of the regulated gene or in distal intergenic space — and a single enhancer can act on multiple genes (or *vice versa*). H3K27ac, which marks *active* enhancers, therefore necessarily produces a long tail of peaks far from any annotated TSS. |

---

## Step 9: Biological Interpretation of Histone Marks

(No code — written interpretation drawing on the Steps 1–8 results.)

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | How does the genomic distribution of H3K4me3 differ from H3K27ac? | H3K4me3 is sharply localised to active TSSs (typically >70% of peaks fall within ±1 kb of a TSS, with narrow peaks ~200–500 bp wide). H3K27ac is broader and more evenly spread: it co-marks active promoters (so a substantial promoter share — 55% in our data) **but also** marks active enhancers in introns and distal intergenic space (~37% of our peaks). The peak shape also differs: H3K4me3 produces narrow, "spiky" peaks; H3K27ac produces broad domains 1–2 kb wide that can extend to >20 kb at super-enhancers. |
| 2 | Why is H3K4me3 strongly associated with TSSs? | H3K4me3 is deposited by the COMPASS / SET1A/B / MLL methyltransferase complex, which is recruited to chromatin by the RNA Pol II pre-initiation complex. The mark therefore tracks the position of paused/initiating polymerase — i.e. precisely at active gene starts. It then recruits TFIID (via the PHD finger of TAF3) and several chromatin remodellers, reinforcing an open, transcription-permissive state at the TSS. |
| 3 | What regulatory role do H3K27ac-marked enhancers play? | Active enhancers loop to their target promoter through cohesin/CTCF-anchored chromatin contacts, deliver bound transcription factors and co-activators (especially p300/CBP, which itself writes the H3K27ac mark), and amplify the rate of transcription initiation. Cell-type-specific enhancer usage is the principal determinant of cell identity: even though most cells share the same housekeeping promoters, each cell type marks a *distinct* set of enhancers with H3K27ac — and switching this enhancer landscape (e.g. in response to IL-7 signaling) is how a cell changes its gene expression program. |

---

## Step 10: Visualisation in IGV

### Script
```bash
bash 10_make_bigwig.sh
```

Generates RPKM-normalised bigWig coverage tracks via `deepTools bamCoverage`
(`--normalizeUsing RPKM --binSize 10 --extendReads`).

### Output
```
qc/
├── IL7-H_H3k27ac_S295.bw   ⭐ ChIP coverage track
└── IL7-H_igG_S286.bw       ⭐ IgG coverage track
```

### IGV loading procedure
1. Launch IGV → **Genome → mm10**.
2. **File → Load from File** ×3:
   - `qc/IL7-H_H3k27ac_S295.bw`
   - `qc/IL7-H_igG_S286.bw`
   - `peaks/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak`
3. Group autoscale on the two bigWigs (right-click track → *Group Autoscale*).
4. Navigate to a known active T-cell gene (e.g. **Actb**, **Gapdh**, **Cd3e**, **Il7r**) or jump to one of the top peaks from `extra_enhancer_analysis.Rmd`.
5. Take a screenshot showing both ChIP and IgG tracks plus the peak-call BED.

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | What does peak height represent in the genome browser? | Read coverage depth at that position, RPKM-normalised. Higher peaks mean more sequenced fragments are stacked at that locus, i.e. the antibody pulled down more H3K27ac-modified chromatin from that region. |
| 2 | Why are H3K4me3 peaks sharper than H3K27ac peaks? | H3K4me3 is deposited on a small number (typically 1–3) of nucleosomes immediately flanking the TSS, so the chromatin region carrying the mark is physically narrow (~500–1000 bp). H3K27ac, by contrast, decorates the entire body of an active enhancer or promoter — often a chromatin region tens of nucleosomes long, especially at super-enhancers — which produces a broad, "smeary" signal in the browser. |
| 3 | How does IgG signal differ from ChIP in IGV? | IgG signal appears as a **flat, low-amplitude track** with no locus-specific peaks — coverage is approximately uniform across the genome (see how the IgG self-call returned 0 peaks in Step 7). The ChIP track on the same x-axis shows tall, region-specific spikes that line up exactly with the `.broadPeak` calls. The visual contrast between the two tracks is what makes a peak call credible. |

---

## Step 11: Interpretation of Experimental Results

(No code — final write-up.)

### Assignment Q&A

| # | Question | Answer |
|---|----------|--------|
| 1 | What biological conclusions can be drawn from promoter-associated histone marks? | The 55% of H3K27ac peaks at promoters identifies the set of *actively transcribed* genes in the IL7-High condition — a snapshot of the cell's current transcriptional program. By overlaying with the (un-assigned) H3K4me3 promoter set one would distinguish **bivalent** loci (H3K4me3 only — poised) from **fully active** loci (H3K4me3 + H3K27ac — transcribing). Because H3K27ac is reversibly deposited by p300/CBP and erased by HDACs, comparing this profile across IL7 conditions (High vs Low ± dose) would reveal which promoters are dynamically reactivated by IL-7 signaling in T-cells. |
| 2 | How do enhancer histone marks contribute to gene regulation? | H3K27ac at distal/intronic sites identifies the cell's active enhancers — the ~37% of our peaks that lie outside promoters. Enhancers contact their target promoters through 3D loops and integrate inputs from multiple lineage transcription factors (in T-cells: GATA3, TCF1, RUNX1, STAT5 downstream of IL-7R). They are the principal determinants of cell-type-specific expression: housekeeping promoters are largely shared across tissues, but enhancer usage is unique to each cell state. Loss of an enhancer's H3K27ac (e.g. via inhibition of p300 or BRD4) collapses transcription of the connected gene without altering the promoter itself. |
| 3 | Why are chromatin modifications considered epigenetic regulators? | They satisfy the operational definition of epigenetic information: (a) they are **heritable** through DNA replication and cell division (writers re-deposit them on newly assembled nucleosomes), (b) they are **reversible** (specific enzymes erase or replace them), and (c) they alter **gene expression without changing DNA sequence**. Together with DNA methylation and 3D chromatin architecture, they encode a layer of cellular memory that records which genes were active in a parent cell and biases the daughter cells to maintain that program — which is how a single genome supports hundreds of stable, distinct cell types. |

---

## Submission Checklist (per PDF)

- [x] **Peak count summary table** — `qc/peak_count_summary.tsv` (Step 7)
- [ ] **IGV genome-browser screenshot** of one representative peak (manual; after Step 10)
- [ ] **1–2 page written explanation** — assemble from Q&A in Steps 1, 3, 4, 7, 8, 9, 11

---

## Extra (beyond the assignment)

`extra_enhancer_analysis.Rmd` — RMarkdown notebook covering:
1. Top-ranked H3K27ac peaks and their nearest genes.
2. Super-enhancer-style hockey-stick ranking.
3. GO Biological Process over-representation analysis on peak-associated genes (clusterProfiler).

Render with:
```bash
Rscript -e 'rmarkdown::render("extra_enhancer_analysis.Rmd")'
```
