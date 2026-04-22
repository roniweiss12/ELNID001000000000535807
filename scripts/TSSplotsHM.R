library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)
library(circlize)

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
base_dir <- "~/ELNID001000000000535807/analysis/dHICA/output"

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
tss_win <- promoters(tss, upstream=10000, downstream=10000) # ±2kb

# attach gene symbols
eg <- mapIds(org.Hs.eg.db, keys=names(tss_win), column="SYMBOL", keytype="ENTREZID", multiVals="first")
mcols(tss_win)$symbol <- eg
# Ensure tss_win is GRanges, not GRangesList
tss_gr <- as(tss_win, "GRanges")

for(histname in marks){
  print(paste("starting", histname))
  ## H3K27ac hf signal extraction
  
  cov_list_hf <- lapply(bw_files[[histname]][["high"]], function(bw) {
    import(bw, format = "BigWig", which = tss_gr)
  })
  
  # convert imported RleList to numeric matrices with EnrichedHeatmap utilities:
  mat_list_hf <- lapply(cov_list_hf, function(gr) {
    normalizeToMatrix(
      signal = gr,
      target = tss,
      value_col = "score",
      extend = 7000,
      w = 10,
      mean_mode = "w0"
    )
  })
  
  # Compute average signal per gene across all samples
  gene_mean_signal <- sapply(mat_list_hf, function(m) rowMeans(m, na.rm = TRUE))
  gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)
  
  # Rank and pick top 1000 genes
  top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:min(1000, length(gene_mean_signal))]
  mat_list_top_hf <- lapply(mat_list_hf, function(m) m[top_idx, ])
  
  # Average matrices by condition
  avg_mat_high <- (mat_list_top_hf[[1]] + mat_list_top_hf[[2]] + mat_list_top_hf[[3]]) / 3
  
  # Compute a comparable enrichment range manually
  mean_profile_high <- colMeans(avg_mat_high, na.rm = TRUE)
  
  ## hist lf signal extraction
  cov_list_lf <- lapply(bw_files[[histname]][["medium"]], function(bw) {
    import(bw, format = "BigWig", which = tss_gr)
  })
  
  # convert imported RleList to numeric matrices with EnrichedHeatmap utilities:
  mat_list_lf <- lapply(cov_list_lf, function(gr) {
    normalizeToMatrix(
      signal = gr,
      target = tss,
      value_col = "score",
      extend = 7000,
      w = 10,
      mean_mode = "w0"
    )
  })
  
  # Compute average signal per gene across all samples
  gene_mean_signal <- sapply(mat_list_lf, function(m) rowMeans(m, na.rm = TRUE))
  gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)
  
  # Rank and pick top 1000 genes
  top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:min(1000, length(gene_mean_signal))]
  mat_list_top_lf <- lapply(mat_list_lf, function(m) m[top_idx, ])
  
  # Average matrices by condition
  avg_mat_medium <- (mat_list_top_lf[[1]] + mat_list_top_lf[[2]] + mat_list_top_lf[[3]]) / 3
  
  # Compute a comparable enrichment range manually
  mean_profile_medium <- colMeans(avg_mat_medium, na.rm = TRUE)
  
  shared_rng <- range(c(mat_list_top_hf, mat_list_top_lf), na.rm = TRUE)
  col_fun <- colorRamp2(c(shared_rng[1], shared_rng[2]), c("white", "red"))
  
  max_y2 <- max(mean_profile_high, mean_profile_medium, na.rm = TRUE)
  ylim_shared_hist <- c(0, max_y2)

  # Plot averaged heatmaps
  png(paste0("~/ELNID001000000000535807/analysis/predicted_histone_modifications/", histname, "_TSSheatmap_7k.png"), height = 7, width = 6, units = "in", res = 300)
  EnrichedHeatmap(avg_mat_high, name = paste0(histname, " - HF"), col = col_fun,
                  top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                            ylim = ylim_shared_hist)),
                  column_title = "High Fatigue") +
    EnrichedHeatmap(avg_mat_medium, name = paste0(histname, " - LF"), col = col_fun,
                    top_annotation = HeatmapAnnotation(
                      enrich = anno_enriched(
                        ylim = ylim_shared_hist
                      )
                    ),
                    column_title = "Low Fatigue")
  dev.off()
  
  # Calculate average signals across genes for each bin
  profile_hist_high <- colMeans(avg_mat_high, na.rm = TRUE)
  profile_hist_low <- colMeans(avg_mat_medium, na.rm = TRUE)
  
  # Create a vector of positions relative to TSS, assuming 10 bp bins over ± 2000 bp
  positions <- seq(-7000, 7000, length.out = 800) # Midpoint of bins
  ylim_range <- c(0, max(c(profile_hist_high, profile_hist_low), na.rm = TRUE))
  png(paste0("~/ELNID001000000000535807/analysis/predicted_histone_modifications/TSSarea_predicted", histname, ".png"),
      height = 4, width = 6, units = "in", res = 300)
  plot(positions, profile_hist_high, type = "l", col = "#ca0020", lwd = 2,
         xlab = "Position relative to TSS (bp)", ylab = "Average signal",
         main = paste0(histname, " signals around TSS"), ylim = ylim_range)
  lines(positions, profile_hist_low, col = "#f4a582", lwd = 2)
  # Plot average signals overlaid
  legend("topright", legend = c("HF", "LF"),
         col = c("#ca0020", "#f4a582"), lwd = 2)
  dev.off()
  print(paste("finished", histname))
  
  rm(list = ls(pattern = "mat|profile_|gene_mean|top_idx|avg|list|max_y"))
  gc()
}


library(GenomicRanges)
library(rtracklayer)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(circlize)



# ── helper: build normalised matrix for ONE sample, streaming ───────────────
# Crucially: import → normalise → discard raw GRanges immediately
make_norm_mat <- function(bw_path, tss_point, tss_windows, extend, w) {
  sig <- import(bw_path, format="BigWig", which=tss_windows)   # raw import
  mat <- normalizeToMatrix(sig, tss_point,
                           value_col = "score",
                           extend    = extend,
                           w         = w,
                           mean_mode = "w0")
  rm(sig); gc(verbose=FALSE)   # <-- drop the big GRanges right away
  mat
}

# ── ranking pass: use only ONE sample per condition to pick top genes ────────
# This avoids holding all 6 matrices in memory at the same time.
# (Using rep-1 of each condition as the ranking sample is a reasonable proxy.)
cat("Ranking genes using rep-1 of each condition...\n")

rank_mat_high <- make_norm_mat(bw_files$high[1],   tss, tss_gr, EXTEND, BIN_W)
rank_mat_med  <- make_norm_mat(bw_files$medium[1], tss, tss_gr, EXTEND, BIN_W)

gene_score <- rowMeans(rank_mat_high, na.rm=TRUE) +
  rowMeans(rank_mat_med,  na.rm=TRUE)
top_idx    <- order(gene_score, decreasing=TRUE)[1:min(N_TOP, length(gene_score))]

rm(rank_mat_high, rank_mat_med, gene_score); gc()

# ── accumulate averaged matrix ONE sample at a time ─────────────────────────
# Never hold more than 1 full matrix in RAM beyond the running sum.

cat("Accumulating high-fatigue matrices...\n")
avg_mat_high <- NULL
for (bw in bw_files$high) {
  m <- make_norm_mat(bw, tss, tss_gr, EXTEND, BIN_W)[top_idx, ]
  if (is.null(avg_mat_high)) avg_mat_high <- m
  else                        avg_mat_high <- avg_mat_high + m
  rm(m); gc(verbose=FALSE)
}
avg_mat_high <- avg_mat_high / length(bw_files$high)

cat("Accumulating medium-fatigue matrices...\n")
avg_mat_medium <- NULL
for (bw in bw_files$medium) {
  m <- make_norm_mat(bw, tss, tss_gr, EXTEND, BIN_W)[top_idx, ]
  if (is.null(avg_mat_medium)) avg_mat_medium <- m
  else                          avg_mat_medium <- avg_mat_medium + m
  rm(m); gc(verbose=FALSE)
}
avg_mat_medium <- avg_mat_medium / length(bw_files$medium)

# ── plotting (unchanged logic) ───────────────────────────────────────────────
shared_rng <- range(c(avg_mat_high, avg_mat_medium), na.rm=TRUE)
col_fun    <- colorRamp2(c(shared_rng[1], shared_rng[2]), c("white","red"))

mean_profile_high   <- colMeans(avg_mat_high,   na.rm=TRUE)
mean_profile_medium <- colMeans(avg_mat_medium, na.rm=TRUE)
ylim_shared         <- c(0, max(mean_profile_high, mean_profile_medium, na.rm=TRUE))

out_dir <- "~/ELNID001000000000535807/analysis/predicted_histone_modifications"

png(file.path(out_dir, paste0(histname, "_TSSheatmap_7k.png")),
    height=7, width=6, units="in", res=300)
draw(
  EnrichedHeatmap(avg_mat_high, name=paste0(histname," - HF"), col=col_fun,
                  top_annotation=HeatmapAnnotation(
                    enrich=anno_enriched(axis_param=list(side="left"), ylim=ylim_shared)),
                  column_title="High Fatigue") +
    EnrichedHeatmap(avg_mat_medium, name=paste0(histname," - LF"), col=col_fun,
                    top_annotation=HeatmapAnnotation(
                      enrich=anno_enriched(ylim=ylim_shared)),
                    column_title="Low Fatigue")
)
dev.off()

positions  <- seq(-EXTEND, EXTEND, length.out=ncol(avg_mat_high))
ylim_range <- c(0, max(mean_profile_high, mean_profile_medium, na.rm=TRUE))

png(file.path(out_dir, paste0("TSSarea_predicted", histname, ".png")),
    height=4, width=6, units="in", res=300)
plot(positions, mean_profile_high,  type="l", col="#ca0020", lwd=2,
     xlab="Position relative to TSS (bp)", ylab="Average signal",
     main=paste0(histname, " signals around TSS"), ylim=ylim_range)
lines(positions, mean_profile_medium, col="#f4a582", lwd=2)
legend("topright", legend=c("HF","LF"), col=c("#ca0020","#f4a582"), lwd=2)
dev.off()

cat("Done:", histname, "\n")