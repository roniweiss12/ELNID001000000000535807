library(rtracklayer)
library(GenomicRanges)
library(vizzy)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)

samples <- list(
  high = c(
    "03_0KC9_02SLTytgat_High-fatigue-1_ATAC_hs_i29-48_REP1",
    "08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1",
    "09_High-fatigue-3"
  ),
  medium = c(
    "01_0KC8_02SLTytgat_Medium-fatigue-1_ATAC_hs_i28-47_REP1",
    "05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1",
    "06_Medium-fatigue-3"
  )
)
marks <- c("H3K27ac","H3K4me1","H3K4me2","H3K4me3",
           "H3K122ac","H3K27me3","H3K36me3",
           "H3K9ac","H3K9me3","H4K20me1")
base_dir <- "~/ELNID001000000000486355/analysis/dHICA/output"

bw_files <- lapply(marks, function(mark) {
  list(
    high = file.path(base_dir, samples$high, "predicted_bw",
                     paste0("out-R-", mark, ".bw")),
    medium = file.path(base_dir, samples$medium, "predicted_bw",
                       paste0("out-R-", mark, ".bw"))
  )
})
names(bw_files) <- marks

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# get TSS GRanges (one TSS per gene)
genes <- genes(txdb)
tss <- resize(promoters(genes, upstream=0, downstream=1), width=1, fix="start") # TSS point
tss_win <- promoters(tss, upstream=4000, downstream=4000) # ±2kb

# attach gene symbols
eg <- mapIds(org.Hs.eg.db, keys=names(tss_win), column="SYMBOL", keytype="ENTREZID", multiVals="first")
mcols(tss_win)$symbol <- eg
# Ensure tss_win is GRanges, not GRangesList
tss_gr <- as(tss_win, "GRanges")

i = 1
for(hist in marks){
  plot_profiles(
    bigwigs = unlist(bw_files)[i:(i + 5)],
    ranges  = tss_gr, 
    names = names(unlist(bw_files))[i:(i + 5)], 
    title = paste0(hist, " TSS enrichemnt"), 
    colors = c("#bae4b3", "#74c476", "#238b45", "#fdcc8a", "#fc8d59", "#d7301f")
  ) +
    guides(color = guide_legend(nrow = 3, ncol = 2))
  ggsave(paste0("~/ELNID001000000000535807/analysis/predicted_histone_modifications/", hist, "_TSSenrich_4k.png"), plot = last_plot(), width = 8, height = 6, units = "in")
  i = i + 6
}


library(IRanges)

# Load your bigwig and TSS annotations
bw <- import(unlist(bw_files)[1], as = "RleList")
tss <- resize(promoters(genes, upstream=0, downstream=1), width=1, fix="start") # TSS point

# Resize TSS to a window, e.g. ±2kb
window <- 2000
tss_windows <- resize(tss, width = window * 2 + 1, fix = "center")

# Extract signal over TSS windows using deepTools-style scoring
extract_matrix <- function(bw_rle, regions, bins = 400) {
  mat <- matrix(NA, nrow = length(regions), ncol = bins)
  for (i in seq_along(regions)) {
    region <- regions[i]
    chr <- as.character(seqnames(region))
    if (!chr %in% names(bw_rle)) next
    scores <- as.numeric(bw_rle[[chr]][start(region):end(region)])
    if (length(scores) < 2) next
    # Bin the signal
    breaks <- cut(seq_along(scores), bins, labels = FALSE)
    mat[i, ] <- tapply(scores, breaks, mean, na.rm = TRUE)
  }
  mat
}

mat <- extract_matrix(bw, tss_windows)

# TSS enrichment score: center signal / background signal
# Background = average of first and last 10% of bins
n_bins <- ncol(mat)
flank_bins <- round(n_bins * 0.1)

center_bin  <- round(n_bins / 2)
center_signal <- mean(colMeans(mat[, (center_bin - 5):(center_bin + 5)], na.rm = TRUE))
bg_signal     <- mean(colMeans(mat[, c(1:flank_bins, (n_bins - flank_bins):n_bins)], na.rm = TRUE))

tss_enrichment_score <- center_signal / bg_signal
cat("TSS Enrichment Score:", tss_enrichment_score)


scores <- sapply(unlist(bw_files), function(f) {
  bw <- import(f, as = "RleList")
  mat <- extract_matrix(bw, tss_windows)
  
  n <- ncol(mat)
  flank <- round(n * 0.1)
  center <- round(n / 2)
  
  center_sig <- mean(colMeans(mat[, (center - 5):(center + 5)], na.rm = TRUE))
  bg_sig     <- mean(colMeans(mat[, c(1:flank, (n - flank):n)], na.rm = TRUE))
  
  center_sig / bg_sig
})

scores_df <- data.frame(
  sample = basename(unlist(bw_files)),
  TSS_enrichment = scores
)
print(scores_df)