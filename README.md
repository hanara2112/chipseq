# H3K27ac ChIP-seq in IL7-High Mouse T-cells

> **Assignment:** Histone-modification ChIP-seq analysis with IgG control (mm10).
> **Course:** Molecular Biology, Spring 2026 — IIIT Hyderabad.
> **Student:** Aryaman Bahl (2022113010).

This repository contains the complete analysis for one ChIP library (H3K27ac, sample **S295**) and one matched IgG control (sample **S286**) — both from IL7-stimulated mouse T-cells, 1 M paired-end reads each. The pipeline takes raw FASTQs through FastQC → Trim Galore → Bowtie2 → samtools markdup → MAPQ filter → MACS2 → ChIPseeker → IGV, plus a bonus enhancer-landscape notebook.

---

## 📄 Submission deliverables (start here)

| # | What | Where |
|---|---|---|
| 1 | **Written report (1–2 pages)** | [`report.pdf`](report.pdf) — sources in [`report.tex`](report.tex) |
| 2 | **Peak-count summary table** | [`results/peak_count_summary.tsv`](results/peak_count_summary.tsv) |
| 3 | **IGV genome-browser screenshot** | [`igv.png`](igv.png) |
| ★ | **Bonus: enhancer-analysis notebook (rendered)** | [`docs/extra_enhancer_analysis.html`](docs/extra_enhancer_analysis.html) — top peaks + super-enhancer ranking + GO BP enrichment. Live preview: <https://hanara2112.github.io/chipseq/docs/extra_enhancer_analysis.html> |
| ★ | **Detailed per-step report with full Q&A** | [`docs/chipseq_analysis.md`](docs/chipseq_analysis.md) |

---

## 📊 Headline results

| Metric | Value |
|---|---|
| Bowtie2 overall alignment rate (S295 / S286) | **95.80 % / 95.33 %** |
| PCR duplication rate (S295 / S286) | 8.14 % / 9.74 % |
| Reads after MAPQ ≥ 20 + properly-paired filter (S295 / S286) | 1.66 M / 1.46 M |
| **H3K27ac peaks** (MACS2 broad, q ≤ 0.05, IgG-controlled) | **18,681** |
| **IgG peaks** (MACS2 self-call, no control) | **0** ← ideal |
| Mean peak fold-enrichment over IgG | 4.54 × (max 17.12 ×) |
| Median peak length | 1,218 bp |
| Peak distribution: promoter / intergenic / intron / exon-UTR | 55.34 % / 18.82 % / 18.08 % / 8.05 % |

The IgG self-call returning **0** peaks at the same statistical threshold as the ChIP confirms that all 18,681 H3K27ac calls are genuine antibody-driven enrichment, not technical artefact.

---

## 🗂️ Directory tour

```
project/
├── README.md                  ← you are here
├── ChIP Project.pdf           ← original assignment brief
├── report.pdf  /  report.tex  ← submitted 1–2 page report (compiled + source)
├── igv.png                    ← submitted IGV screenshot (Figure 1 of report)
├── .gitignore
│
├── scripts/                   ← entire pipeline as numbered shell + R scripts
│   ├── 01_fastqc.sh                                ── FastQC + MultiQC
│   ├── 02_trim_adapters.sh                          ── Trim Galore, Q20
│   ├── 03_align_to_genome.sh                        ── Bowtie2 → mm10
│   ├── 04_remove_pcr_duplicates.sh                  ── samtools fixmate + markdup
│   ├── 05_filter_low_mapq_alignments.sh             ── MAPQ ≥ 20, properly-paired
│   ├── 06_call_peaks_macs2_h3k27ac_vs_igg.sh        ── MACS2 broad, IgG-controlled
│   ├── 07_count_peaks_chip_vs_igg.sh                ── peak-count summary
│   ├── 08_run_chipseeker.sh                         ── wrapper for the R script
│   ├── 08_chipseeker_annotate.R                     ── peak → gene annotation
│   ├── 10_make_bigwig.sh                            ── deepTools bamCoverage
│   └── extra_enhancer_analysis.Rmd                  ── bonus notebook (R Markdown)
│
├── docs/                      ← long-form documentation
│   ├── chipseq_analysis.md            ← per-step results + every assignment Q&A
│   ├── implementation_plan.md         ← initial planning doc
│   └── extra_enhancer_analysis.html   ← rendered bonus notebook (open in browser)
│
├── results/                   ← all outputs an evaluator wants to see
│   ├── peak_count_summary.tsv         ⭐ submission deliverable #2
│   ├── peaks/
│   │   ├── H3K27ac/
│   │   │   ├── *_peaks.broadPeak      ⭐ 18,681 peaks
│   │   │   ├── *_peaks.gappedPeak     ── bed12 with sub-peak structure
│   │   │   ├── *_peaks.xls            ── full peak table with stats
│   │   │   ├── *_macs2.log            ── MACS2 run log
│   │   │   └── H3K27ac_annotated.csv  ── per-peak nearest-gene annotation
│   │   └── IgG/
│   │       └── *_peaks.narrowPeak     ── 0 peaks (file exists, empty)
│   ├── qc/
│   │   ├── multiqc_report.html        ── aggregated FastQC report
│   │   ├── chipseeker_summary.txt     ── % promoter / intergenic / intron / exon
│   │   ├── chipseeker_H3K27ac.pdf     ── pie + bar + dist-to-TSS plots
│   │   ├── *_bowtie2.log              ── alignment stats
│   │   ├── *_dup_metrics.txt          ── samtools markdup output
│   │   └── *_filter_stats.txt         ── before/after MAPQ filter
│   └── coverage/                      (gitignored — bigWig tracks for IGV, regenerable)
│
└── data/                      ← raw + intermediate data (gitignored, regenerable)
    ├── student_dataset/       ── raw FASTQs (24 files, ~1.5 GB)
    ├── genome/                ── mm10 Bowtie2 index (~3.5 GB)
    │   └── make_mm10.sh        ← committed: index-build script
    ├── trimmed/               ── post-trim FASTQs + Trim Galore reports
    └── aligned/               ── BAMs (raw, dedup, filtered) + bowtie2 logs
```

The principle: **`results/` and `docs/` are everything an evaluator needs to read; `data/` is regenerable bulk and is gitignored.**

---

## ▶️ How to reproduce

Prerequisites: a conda environment named `chipseq-env` with `fastqc`, `multiqc`, `trim-galore`, `bowtie2`, `samtools`, `macs2`, `deeptools`, and R + Bioconductor (`ChIPseeker`, `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`, `clusterProfiler`).

```bash
conda activate chipseq-env

bash scripts/01_fastqc.sh                                 # ~3 min
bash scripts/02_trim_adapters.sh                           # ~5 min
bash scripts/03_align_to_genome.sh                         # ~10 min
bash scripts/04_remove_pcr_duplicates.sh                   # ~3 min
bash scripts/05_filter_low_mapq_alignments.sh              # ~1 min
bash scripts/06_call_peaks_macs2_h3k27ac_vs_igg.sh         # ~30 s
bash scripts/07_count_peaks_chip_vs_igg.sh                 # instant
bash scripts/08_run_chipseeker.sh                          # ~1 min (+ first-run package install)
bash scripts/10_make_bigwig.sh                             # ~2 min

# bonus enhancer / GO analysis
Rscript -e 'rmarkdown::render("scripts/extra_enhancer_analysis.Rmd",
                              output_dir = "docs")'
```

Each script is idempotent — re-running overwrites its outputs. Raw data must be placed in `data/student_dataset/` and the mm10 Bowtie2 index in `data/genome/` (run `bash data/genome/make_mm10.sh` to fetch+build the index from scratch).

---

## 🔬 What the analysis shows (one-paragraph summary)

H3K27ac is the canonical molecular signature of active regulatory chromatin — deposited by p300/CBP at both active promoters and active enhancers. In IL7-stimulated mouse T-cells we recovered **18,681 H3K27ac peaks** at q ≤ 0.05 against **zero** peaks in the matched IgG self-call: an ideal specificity outcome. The peaks distribute as **55 % promoter / 19 % distal intergenic / 18 % intron / 8 % exon-UTR** — sitting between the textbook profile of a pure promoter mark (H3K4me3) and a pure enhancer mark, and faithfully reflecting H3K27ac's dual role at both classes of regulatory element. The promoter-proximal majority identifies actively transcribed genes; the ~37 % intronic + intergenic pool identifies active enhancers — the chromatin regions where lineage transcription factors and IL-7-driven STAT5 converge to specify T-cell state. Full per-step methodology, every assignment question, and the biological discussion are in [`docs/chipseq_analysis.md`](docs/chipseq_analysis.md) and [`report.pdf`](report.pdf).

---

## 🛠 Tool versions

FastQC 0.12.1 · MultiQC 1.34 · Trim Galore (Cutadapt 5.2) · Bowtie2 (mm10 pre-built index) · samtools 1.23.1 · MACS2 2.2.9.1 · deepTools 3.5.6 · ChIPseeker (Bioconductor release 3.22) · clusterProfiler (Bioconductor)
