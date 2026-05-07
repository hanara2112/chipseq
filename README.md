# H3K27ac ChIP-seq in IL7-stimulated Mouse T-cells

End-to-end ChIP-seq analysis of one H3K27ac library (S295) and a matched IgG control (S286) from IL7-stimulated mouse T-cells. Raw FASTQs through FastQC → Trim Galore → Bowtie2 (mm10) → samtools markdup → MAPQ filter → MACS2 → ChIPseeker → IGV, plus a bonus enhancer-landscape notebook.

> Molecular Biology, Spring 2026 — IIIT Hyderabad · Aryaman Bahl (2022113010)

---

## What the data say

H3K27ac is the canonical mark of active regulatory chromatin — deposited by p300/CBP at both active promoters and active enhancers. In IL7-stimulated mouse T-cells we recovered **18,681 H3K27ac peaks** at q ≤ 0.05 against **zero** peaks in the matched IgG self-call: an ideal specificity outcome. Peaks distribute as **55 % promoter / 19 % distal intergenic / 18 % intron / 8 % exon-UTR** — between the textbook profile of a pure promoter mark (H3K4me3) and a pure enhancer mark, faithfully reflecting H3K27ac's dual role. The promoter-proximal majority marks actively transcribed genes; the ~37 % intronic + intergenic pool marks active enhancers — chromatin regions where lineage transcription factors and IL-7-driven STAT5 converge to specify T-cell state.

---

## Headline numbers

| Metric | Value |
|---|---|
| Bowtie2 alignment rate (S295 / S286) | 95.80 % / 95.33 % |
| PCR duplication rate (S295 / S286) | 8.14 % / 9.74 % |
| Reads after MAPQ ≥ 20 + properly-paired (S295 / S286) | 1.66 M / 1.46 M |
| **H3K27ac peaks** (MACS2 broad, q ≤ 0.05, IgG-controlled) | **18,681** |
| **IgG peaks** (self-call, no control) | **0** |
| Mean / max peak fold-enrichment over IgG | 4.54 × / 17.12 × |
| Median peak length | 1,218 bp |
| Promoter / intergenic / intron / exon-UTR | 55.34 / 18.82 / 18.08 / 8.05 % |

The IgG self-call returning zero peaks at the same threshold confirms all 18,681 H3K27ac calls are genuine antibody-driven enrichment, not technical artefact.

---

## Deliverables

- **Report (1–2 pages)** — [`report.pdf`](report.pdf) · sources [`report.tex`](report.tex)
- **Peak-count summary** — [`results/peak_count_summary.tsv`](results/peak_count_summary.tsv)
- **IGV screenshot** — [`igv.png`](igv.png)
- **Bonus: enhancer notebook** — top peaks, super-enhancer ranking, GO BP enrichment. [`docs/extra_enhancer_analysis.html`](docs/extra_enhancer_analysis.html) · [live preview](https://hanara2112.github.io/chipseq/docs/extra_enhancer_analysis.html)
- **Per-step report with full Q&A** — [`docs/chipseq_analysis.md`](docs/chipseq_analysis.md)

---

## Stack

FastQC · MultiQC · Trim Galore · Bowtie2 · samtools · MACS2 · deepTools · ChIPseeker · clusterProfiler · IGV. Pipeline lives in [`scripts/`](scripts/) as numbered shell + R steps; reproduction instructions are in [`docs/chipseq_analysis.md`](docs/chipseq_analysis.md).
