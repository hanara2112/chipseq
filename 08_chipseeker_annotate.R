# Step 8 — Genomic distribution of H3K27ac peaks (ChIPseeker)
# Inputs:  peaks/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak
# Outputs: peaks/H3K27ac/H3K27ac_annotated.csv
#          qc/chipseeker_H3K27ac.pdf
#          qc/chipseeker_summary.txt

suppressPackageStartupMessages({

  # auto-install missing Bioconductor packages
  needed <- c("ChIPseeker",
              "TxDb.Mmusculus.UCSC.mm10.knownGene",
              "org.Mm.eg.db")
  missing <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }

  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

PROJECT  <- path.expand("~/IIITH/S26/mb/project")
PEAK_FN  <- file.path(PROJECT, "peaks/H3K27ac/IL7-H_H3k27ac_S295_peaks.broadPeak")
OUT_CSV  <- file.path(PROJECT, "peaks/H3K27ac/H3K27ac_annotated.csv")
OUT_PDF  <- file.path(PROJECT, "qc/chipseeker_H3K27ac.pdf")
OUT_TXT  <- file.path(PROJECT, "qc/chipseeker_summary.txt")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

cat(">>> reading peaks:", PEAK_FN, "\n")
peaks <- readPeakFile(PEAK_FN)
cat("    n =", length(peaks), "peaks\n")

cat(">>> annotating against mm10 TxDb (TSS region: -3000 / +3000)\n")
anno <- annotatePeak(peaks,
                     TxDb     = txdb,
                     tssRegion = c(-3000, 3000),
                     annoDb   = "org.Mm.eg.db",
                     verbose  = FALSE)

# --- write annotated peaks table -----------------------------------------
df <- as.data.frame(anno)
write.csv(df, OUT_CSV, row.names = FALSE)
cat(">>> wrote:", OUT_CSV, "\n")

# --- summary stats (for assignment Q&A) ----------------------------------
feat <- anno@annoStat                 # data.frame: Feature, Frequency
feat$Frequency <- round(feat$Frequency, 2)

promoter_pct  <- sum(feat$Frequency[grepl("Promoter", feat$Feature)])
intergenic_pct <- sum(feat$Frequency[grepl("Intergenic", feat$Feature)])
intron_pct    <- sum(feat$Frequency[grepl("Intron", feat$Feature)])
exon_pct      <- sum(feat$Frequency[grepl("Exon", feat$Feature)])

sink(OUT_TXT)
cat("=== ChIPseeker summary — H3K27ac (S295) ===\n")
cat("Peaks annotated:           ", length(peaks), "\n\n")
cat("--- Genomic feature distribution (%) ---\n")
print(feat, row.names = FALSE)
cat("\n--- aggregated ---\n")
cat(sprintf("Promoter  : %.2f %%\n", promoter_pct))
cat(sprintf("Intergenic: %.2f %%\n", intergenic_pct))
cat(sprintf("Intron    : %.2f %%\n", intron_pct))
cat(sprintf("Exon      : %.2f %%\n", exon_pct))
sink()
cat(">>> wrote:", OUT_TXT, "\n")

# --- plots ---------------------------------------------------------------
pdf(OUT_PDF, width = 10, height = 7)
plotAnnoPie(anno)
plotAnnoBar(anno)
plotDistToTSS(anno, title = "H3K27ac peaks — distance to TSS")
dev.off()
cat(">>> wrote:", OUT_PDF, "\n")

cat("=== Step 8 complete ===\n")
