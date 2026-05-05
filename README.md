# H3K27ac ChIP-seq — IL7-High condition (mm10)

ChIP-seq analysis assignment, MolBio S26 (IIIT-H).
Histone modification ChIP-seq with IgG control, mouse genome (mm10).

📊 **[View bonus enhancer analysis (rendered notebook)](https://hanara2112.github.io/chipseq/extra_enhancer_analysis.html)** — top peaks, super-enhancer ranking, GO enrichment.

> *Note:* the link works once GitHub Pages is enabled for this repo. **Settings → Pages → Source: `main`, folder: `/` (root) → Save.** The page will then be served from `https://hanara2112.github.io/chipseq/extra_enhancer_analysis.html`.

## Assigned samples

| Role | Sample | Library | Condition | Mark |
|---|---|---|---|---|
| ChIP | **S295** | `IL7-H_H3k27ac_S295` | IL7-High | **H3K27ac** (active enhancer mark) |
| Control | **S286** | `IL7-H_igG_S286` | IL7-High | IgG (non-specific) |

Each library: 1 M paired-end reads, 151 bp.

## Agenda — what this project does

1. **QC the raw reads** (Step 1) → confirm library quality + flag adapter contamination.
2. **Trim adapters + low-quality bases** (Step 2).
3. **Align to mm10** with Bowtie2 (Step 3).
4. **Remove PCR duplicates** with `samtools markdup` (Step 4).
5. **Filter low-quality alignments** — MAPQ ≥ 20 + properly-paired only (Step 5).
6. **Call peaks** with MACS2 (`--broad` for H3K27ac) using IgG as control (Step 6).
7. **Compare ChIP vs IgG peak counts** to confirm enrichment is real (Step 7).
8. **Annotate peaks** to genomic features (promoter / intron / intergenic) with ChIPseeker (Step 8).
9. **Interpret biologically** (Steps 9, 11).
10. **Visualise in IGV** with bigWig coverage tracks (Step 10).

Full per-step results, parameters, and assignment Q&A: see [`chipseq_analysis.md`](chipseq_analysis.md).

## Pipeline

```
Raw FASTQ → FastQC → Trim Galore → Bowtie2 → samtools markdup → MAPQ20 filter
            (01)      (02)         (03)         (04)            (05)

         → MACS2 (broad) ──┬─→ peak count vs IgG (07)
            (06)           ├─→ ChIPseeker annotation (08)
                           └─→ deepTools bigWig → IGV (10)

Steps 9, 11: written biological interpretation (no code).
extra_enhancer_analysis.Rmd: super-enhancer ranking + GO enrichment (bonus).
```

## Directory structure

```
project/
├── README.md                      ← you are here
├── chipseq_analysis.md            ← full per-step report + assignment Q&A
├── ChIP Project.pdf               ← original assignment
├── implementation_plan.md.resolved← initial plan (pre-narrowing)
│
├── student_dataset/               ← raw FASTQs (1M paired-end reads each)
├── genome/                        ← mm10 Bowtie2 index
│
├── 01_fastqc.sh                   ← Step 1: FastQC + MultiQC
├── 02_trim_adapters.sh            ← Step 2: Trim Galore (Q20)
├── 03_align_to_genome.sh          ← Step 3: Bowtie2 → sorted BAM
├── 04_remove_pcr_duplicates.sh    ← Step 4: fixmate + markdup
├── 05_filter_low_mapq_alignments.sh ← Step 5: MAPQ≥20, properly-paired
├── 06_call_peaks_macs2_h3k27ac_vs_igg.sh ← Step 6: MACS2 broad
├── 07_count_peaks_chip_vs_igg.sh  ← Step 7: peak count summary
├── 08_chipseeker_annotate.R       ← Step 8: peak annotation (R)
├── 08_run_chipseeker.sh           ← Step 8 wrapper
├── 10_make_bigwig.sh              ← Step 10: deepTools bamCoverage
├── extra_enhancer_analysis.Rmd    ← bonus: super-enhancers + GO enrichment
│
├── qc/
│   ├── fastqc_raw/multiqc_report.html      ← Step 1 report
│   ├── *_dup_metrics.txt                   ← Step 4 dedup metrics
│   ├── *_filter_stats.txt                  ← Step 5 filter stats
│   ├── peak_count_summary.tsv              ← Step 7 — submission deliverable
│   ├── chipseeker_summary.txt              ← Step 8 — % promoter/intergenic
│   ├── chipseeker_H3K27ac.pdf              ← Step 8 plots
│   └── *.bw                                ← Step 10 IGV coverage tracks
│
├── trimmed/                       ← Step 2: trimmed FASTQs + reports
├── aligned/
│   ├── *.bam (+ .bai)             ← Step 3: raw BAMs + bowtie2 logs
│   ├── *_deduplicated.bam         ← Step 4: dedup BAMs
│   └── *_filtered_MAPQ20.bam      ← Step 5: analysis-ready BAMs ⭐
└── peaks/
    ├── H3K27ac/
    │   ├── *_peaks.broadPeak           ⭐ 18,681 peaks (Step 6)
    │   ├── *_peaks.gappedPeak
    │   ├── *_peaks.xls
    │   └── H3K27ac_annotated.csv       ← Step 8 per-peak gene annotation
    └── IgG/
        └── *_peaks.narrowPeak          ← 0 peaks (Step 6 self-call)
```

## How to run

```bash
conda activate chipseq-env

bash 01_fastqc.sh                                # ~3 min
bash 02_trim_adapters.sh                          # ~5 min
bash 03_align_to_genome.sh                        # ~10 min
bash 04_remove_pcr_duplicates.sh                  # ~3 min
bash 05_filter_low_mapq_alignments.sh             # ~1 min
bash 06_call_peaks_macs2_h3k27ac_vs_igg.sh        # ~30 s
bash 07_count_peaks_chip_vs_igg.sh                # instant
bash 08_run_chipseeker.sh                         # ~1 min (+ first-run install)
bash 10_make_bigwig.sh                            # ~2 min

# bonus
Rscript -e 'rmarkdown::render("extra_enhancer_analysis.Rmd")'
```

Each script is idempotent — re-running overwrites its outputs.

## Headline results

| Metric | Value |
|---|---|
| Reads aligned (S295 / S286) | 95.80% / 95.33% |
| Duplicate rate (S295 / S286) | 8.14% / 9.74% |
| Reads after MAPQ20 filter (S295 / S286) | 1.66 M / 1.46 M |
| **H3K27ac peaks** (broad, q≤0.05, IgG-controlled) | **18,681** |
| **IgG peaks** (self-call, no control) | **0** |
| Mean peak fold-enrichment over IgG | 4.54× |
| Median peak length | 1,218 bp |
| Peaks at promoters (≤3 kb from TSS) | 55.34% |
| Peaks at distal intergenic regions | 18.82% |

## Submission deliverables (per assignment PDF)

- [x] Peak count summary table → `qc/peak_count_summary.tsv`
- [ ] IGV screenshot of one representative peak → manual after Step 10
- [ ] 1–2 page written explanation → assemble from Q&A in [`chipseq_analysis.md`](chipseq_analysis.md)

## Environment

Conda env: `chipseq-env` — fastqc, multiqc, trim-galore, bowtie2, samtools,
macs2, deeptools, R + Bioconductor (ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene,
org.Mm.eg.db, clusterProfiler).
# chipseq
