library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # change as needed
library(org.Hs.eg.db)

atac_bw_files <- c(High.Fatigue1="~/ELNID001000000000486355/ATACseq/BigWig/03_0KC9_02SLTytgat_High-fatigue-1_ATAC_hs_i29-48_REP1.mLb.clN.bigWig", 
              High.Fatigue2="~/ELNID001000000000486355/Q52867_ATACseq/BigWig/08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1.mLb.clN.bigWig", 
              High.Fatigue3="~/ELNID001000000000486355/Q52867_ATACseq/BigWig/09_0KXR_02UCAUMC_High-fatigue-3_ATAC_i7-15_REP1.mLb.clN.bigWig", 
              Medium.Fatigue1="~/ELNID001000000000486355/ATACseq/BigWig/01_0KC8_02SLTytgat_Medium-fatigue-1_ATAC_hs_i28-47_REP1.mLb.clN.bigWig", 
              Medium.Fatigue2="~/ELNID001000000000486355/Q52867_ATACseq/BigWig/05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1.mLb.clN.bigWig",
              Medium.Fatigue3="~/ELNID001000000000486355/Q52867_ATACseq/BigWig/06_0KXP_02UCAUMC_Medium-fatigue-3_ATAC_i5-521_REP1.mLb.clN.bigWig")
ct_bw_files <- c(High.Fatigue1="~/ELNID001000000000486355/CUT_TAG/BigWig/03_0KCB_02SMTytgat_High-fatigue-1_H3Q5ser_hs_i702-508_R1.bigWig", 
                 Medium.Fatigue1="~/ELNID001000000000486355/CUT_TAG/BigWig/01_0KCA_02SMTytgat_Medium-fatigue-1_H3Q5ser_hs_i701-508_R1.bigWig")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 1) get TSS GRanges (one TSS per gene)
genes <- genes(txdb)
tss <- resize(promoters(genes, upstream=0, downstream=1), width=1, fix="start") # TSS point
tss_win <- promoters(tss, upstream=2000, downstream=2000) # ±2kb

# attach gene symbols
eg <- mapIds(org.Hs.eg.db, keys=names(tss_win), column="SYMBOL", keytype="ENTREZID", multiVals="first")
mcols(tss_win)$symbol <- eg

##ATAC
# Ensure tss_win is GRanges, not GRangesList
tss_gr <- as(tss_win, "GRanges")

cov_list <- lapply(atac_bw_files, function(bw) {
  import(bw, format = "BigWig", which = tss_gr)
})
# convert imported RleList to numeric matrices with EnrichedHeatmap utilities:
mat_list <- lapply(cov_list, function(gr) {
  normalizeToMatrix(
    signal = gr,
    target = tss_win,
    value_col = "score",
    extend = 2000,
    w = 10,
    mean_mode = "w0"
  )
})

# Compute average signal per gene across all samples
gene_mean_signal <- sapply(mat_list, function(m) rowMeans(m, na.rm = TRUE))
gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)

# Rank and pick top 1000 genes
top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:1000]

# Subset matrices to these genes
mat_list_top <- lapply(mat_list, function(m) m[top_idx, ])

# Confirm
sapply(mat_list_top, nrow)

# Average matrices by condition
avg_mat_high <- (mat_list_top$High.Fatigue1 + mat_list_top$High.Fatigue2 + mat_list_top$High.Fatigue3) / 3
avg_mat_med  <- (mat_list_top$Medium.Fatigue1 + mat_list_top$Medium.Fatigue2 + mat_list_top$Medium.Fatigue3) / 3

# Compute a comparable enrichment range manually
mean_profile_high <- colMeans(avg_mat_high, na.rm = TRUE)
mean_profile_med  <- colMeans(avg_mat_med, na.rm = TRUE)

# Determine the shared maximum (and minimum, usually 0)
max_y <- max(mean_profile_high, mean_profile_med, na.rm = TRUE)
ylim_shared <- c(0, max_y)

# Plot averaged heatmaps
png("~/ELNID001000000000486355/DiffAnalysis/ATAC/TSSheatmap.png", height = 7, width = 6, units = "in", res = 300)
EnrichedHeatmap(avg_mat_high, name = "ATAC-seq - HF", col = c("white", "red"),
                top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                          ylim = ylim_shared)),
                column_title = "High Fatigue") +
  EnrichedHeatmap(avg_mat_med, name = "ATAC-seq - LF", col = c("white", "red"),
                  top_annotation = HeatmapAnnotation(
                    enrich = anno_enriched(
                      ylim = ylim_shared
                    )
                  ),
                  column_title = "Low Fatigue")
dev.off()


##high fatigue oriented:

# Compute average signal per gene across all samples
gene_mean_signal <- sapply(mat_list[1:3], function(m) rowMeans(m, na.rm = TRUE))
gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)

# Rank and pick top 1000 genes
top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:1000]

# Subset matrices to these genes
mat_list_top <- lapply(mat_list, function(m) m[top_idx, ])

# Confirm
sapply(mat_list_top, nrow)

# Average matrices by condition
avg_mat_high <- (mat_list_top$High.Fatigue1 + mat_list_top$High.Fatigue2 + mat_list_top$High.Fatigue3) / 3
avg_mat_med  <- (mat_list_top$Medium.Fatigue1 + mat_list_top$Medium.Fatigue2 + mat_list_top$Medium.Fatigue3) / 3

# Compute a comparable enrichment range manually
mean_profile_high <- colMeans(avg_mat_high, na.rm = TRUE)
mean_profile_med  <- colMeans(avg_mat_med, na.rm = TRUE)

# Determine the shared maximum (and minimum, usually 0)
max_y <- max(mean_profile_high, mean_profile_med, na.rm = TRUE)
ylim_shared <- c(0, max_y)

# Plot averaged heatmaps
png("~/ELNID001000000000486355/DiffAnalysis/ATAC/TSSheatmap_HF.png", height = 7, width = 6, units = "in", res = 300)
EnrichedHeatmap(avg_mat_high, name = "ATAC-seq - HF", col = c("white", "red"),
                top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                          ylim = ylim_shared)),
                column_title = "High Fatigue") +
  EnrichedHeatmap(avg_mat_med, name = "ATAC-seq - LF", col = c("white", "red"),
                  top_annotation = HeatmapAnnotation(
                    enrich = anno_enriched(
                      ylim = ylim_shared
                    )
                  ),
                  column_title = "Low Fatigue")
dev.off()


##low fatigue oriented:
# Compute average signal per gene across all samples
gene_mean_signal <- sapply(mat_list[4:6], function(m) rowMeans(m, na.rm = TRUE))
gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)

# Rank and pick top 1000 genes
top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:1000]

# Subset matrices to these genes
mat_list_top <- lapply(mat_list, function(m) m[top_idx, ])

# Confirm
sapply(mat_list_top, nrow)

# Average matrices by condition
avg_mat_high <- (mat_list_top$High.Fatigue1 + mat_list_top$High.Fatigue2 + mat_list_top$High.Fatigue3) / 3
avg_mat_med  <- (mat_list_top$Medium.Fatigue1 + mat_list_top$Medium.Fatigue2 + mat_list_top$Medium.Fatigue3) / 3

# Compute a comparable enrichment range manually
mean_profile_high <- colMeans(avg_mat_high, na.rm = TRUE)
mean_profile_med  <- colMeans(avg_mat_med, na.rm = TRUE)

# Determine the shared maximum (and minimum, usually 0)
max_y <- max(mean_profile_high, mean_profile_med, na.rm = TRUE)
ylim_shared <- c(0, max_y)

# Plot averaged heatmaps
png("~/ELNID001000000000486355/DiffAnalysis/ATAC/TSSheatmap_LF.png", height = 7, width = 6, units = "in", res = 300)
EnrichedHeatmap(avg_mat_high, name = "ATAC-seq - HF", col = c("white", "red"),
                top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                          ylim = ylim_shared)),
                column_title = "High Fatigue") +
  EnrichedHeatmap(avg_mat_med, name = "ATAC-seq - LF", col = c("white", "red"),
                  top_annotation = HeatmapAnnotation(
                    enrich = anno_enriched(
                      ylim = ylim_shared
                    )
                  ),
                  column_title = "Low Fatigue")
dev.off()


# ---- Identify top accessible TSS regions ----

# 1) Compute the mean signal per gene (already done above)
gene_mean_signal <- sapply(mat_list, function(m) rowMeans(m, na.rm = TRUE))
gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)

# 2) Get indices of top 50 most accessible genes/TSSs
top_idx_50 <- order(gene_mean_signal, decreasing = TRUE)[1:50]

# 3) Extract the corresponding GRanges (TSS ±2kb window)
top_tss_regions <- tss_win[top_idx_50]

# 4) Add accessibility values
top_tss_df <- data.frame(
  seqnames = as.character(seqnames(top_tss_regions)),
  start = start(top_tss_regions),
  end = end(top_tss_regions),
  symbol = mcols(top_tss_regions)$symbol,
  entrez_id = names(top_tss_regions),
  mean_accessibility = gene_mean_signal[top_idx_50]
)

# 5) Inspect or export
head(top_tss_df, 10)
write.csv(top_tss_df, "~/ELNID001000000000486355/DiffAnalysis/ATAC/top50_accessible_TSS.csv", row.names = FALSE)


##CUT&TAG

# Ensure tss_win is GRanges, not GRangesList
tss_gr <- as(tss_win, "GRanges")

cov_list_ct <- lapply(ct_bw_files, function(bw) {
  import(bw, format = "BigWig", which = tss_gr)
})
# convert imported RleList to numeric matrices with EnrichedHeatmap utilities:
mat_list_ct <- lapply(cov_list_ct, function(gr) {
  normalizeToMatrix(
    signal = gr,
    target = tss_win,
    value_col = "score",
    extend = 2000,
    w = 10,
    mean_mode = "w0"
  )
})

# Compute average signal per gene across all samples
gene_mean_signal <- sapply(mat_list_ct, function(m) rowMeans(m, na.rm = TRUE))
gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)

# Rank and pick top 1000 genes
top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:50]

# Subset matrices to these genes
mat_list_top_ct <- lapply(mat_list_ct, function(m) m[top_idx, ])

# Confirm
sapply(mat_list_top_ct, nrow)

# Compute a comparable enrichment range manually
mean_profile_high <- colMeans(mat_list_top_ct$High.Fatigue1, na.rm = TRUE)
mean_profile_med  <- colMeans(mat_list_top_ct$Medium.Fatigue1, na.rm = TRUE)

# Determine the shared maximum (and minimum, usually 0)
max_y <- max(mean_profile_high, mean_profile_med, na.rm = TRUE)
ylim_shared <- c(0, max_y)

# Plot averaged heatmaps
png("~/ELNID001000000000486355/DiffAnalysis/CUT_TAG/TSSheatmap_top50.png", height = 7, width = 6, units = "in", res = 300)
ht <- EnrichedHeatmap(mat_list_top_ct$High.Fatigue1, name = "CUT&TAG - HF", col = c("white", "red"),
                      top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                                ylim = ylim_shared)),
                      column_title = "High Fatigue") +
  EnrichedHeatmap(mat_list_top_ct$Medium.Fatigue1, name = "CUT&TAG - LF", col = c("white", "red"),
                  top_annotation = HeatmapAnnotation(
                    enrich = anno_enriched(
                      ylim = ylim_shared
                    )
                  ), 
                  column_title = "Low Fatigue")
# Draw with extra left padding (top, right, bottom, left)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"))  # <-- increase left padding (last value)
dev.off()

#high oriented
# Compute average signal per gene across all samples
gene_mean_signal <- rowMeans(mat_list_ct$High.Fatigue1, na.rm = TRUE)

# Rank and pick top 1000 genes
top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:1000]

# Subset matrices to these genes
mat_list_top_ct <- lapply(mat_list_ct, function(m) m[top_idx, ])

# Confirm
sapply(mat_list_top_ct, nrow)

# Compute a comparable enrichment range manually
mean_profile_high <- colMeans(mat_list_top_ct$High.Fatigue1, na.rm = TRUE)
mean_profile_med  <- colMeans(mat_list_top_ct$Medium.Fatigue1, na.rm = TRUE)

# Determine the shared maximum (and minimum, usually 0)
max_y <- max(mean_profile_high, mean_profile_med, na.rm = TRUE)
ylim_shared <- c(0, max_y)

# Plot averaged heatmaps
png("~/ELNID001000000000486355/DiffAnalysis/CUT_TAG/TSSheatmap_HF.png", height = 7, width = 6, units = "in", res = 300)
ht <- EnrichedHeatmap(mat_list_top_ct$High.Fatigue1, name = "CUT&TAG - HF", col = c("white", "red"),
                      top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                                ylim = ylim_shared)),
                      column_title = "High Fatigue") +
  EnrichedHeatmap(mat_list_top_ct$Medium.Fatigue1, name = "CUT&TAG - LF", col = c("white", "red"),
                  top_annotation = HeatmapAnnotation(
                    enrich = anno_enriched(
                      ylim = ylim_shared
                    )
                  ), 
                  column_title = "Low Fatigue")
# Draw with extra left padding (top, right, bottom, left)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"))  # <-- increase left padding (last value)
dev.off()

#low oriented
# Compute average signal per gene across all samples
gene_mean_signal <- rowMeans(mat_list_ct$Medium.Fatigue1, na.rm = TRUE)

# Rank and pick top 1000 genes
top_idx <- order(gene_mean_signal, decreasing = TRUE)[1:1000]

# Subset matrices to these genes
mat_list_top_ct <- lapply(mat_list_ct, function(m) m[top_idx, ])

# Confirm
sapply(mat_list_top_ct, nrow)

# Compute a comparable enrichment range manually
mean_profile_high <- colMeans(mat_list_top_ct$High.Fatigue1, na.rm = TRUE)
mean_profile_med  <- colMeans(mat_list_top_ct$Medium.Fatigue1, na.rm = TRUE)

# Determine the shared maximum (and minimum, usually 0)
max_y <- max(mean_profile_high, mean_profile_med, na.rm = TRUE)
ylim_shared <- c(0, max_y)

# Plot averaged heatmaps
png("~/ELNID001000000000486355/DiffAnalysis/CUT_TAG/TSSheatmap_LF.png", height = 7, width = 6, units = "in", res = 300)
ht <- EnrichedHeatmap(mat_list_top_ct$High.Fatigue1, name = "CUT&TAG - HF", col = c("white", "red"),
                      top_annotation = HeatmapAnnotation(enrich = anno_enriched(axis_param = list(side = "left"),
                                                                                ylim = ylim_shared)),
                      column_title = "High Fatigue") +
  EnrichedHeatmap(mat_list_top_ct$Medium.Fatigue1, name = "CUT&TAG - LF", col = c("white", "red"),
                  top_annotation = HeatmapAnnotation(
                    enrich = anno_enriched(
                      ylim = ylim_shared
                    )
                  ), 
                  column_title = "Low Fatigue")
# Draw with extra left padding (top, right, bottom, left)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"))  # <-- increase left padding (last value)
dev.off()

# ---- Identify top accessible TSS regions ----

# 1) Compute the mean signal per gene (already done above)
gene_mean_signal <- sapply(mat_list, function(m) rowMeans(m, na.rm = TRUE))
gene_mean_signal <- rowMeans(gene_mean_signal, na.rm = TRUE)

# 2) Get indices of top 50 most accessible genes/TSSs
top_idx_50 <- order(gene_mean_signal, decreasing = TRUE)[1:50]

# 3) Extract the corresponding GRanges (TSS ±2kb window)
top_tss_regions <- tss_win[top_idx_50]

# 4) Add accessibility values
top_tss_df <- data.frame(
  seqnames = as.character(seqnames(top_tss_regions)),
  start = start(top_tss_regions),
  end = end(top_tss_regions),
  symbol = mcols(top_tss_regions)$symbol,
  entrez_id = names(top_tss_regions),
  mean_accessibility = gene_mean_signal[top_idx_50]
)

# 5) Inspect or export
head(top_tss_df, 10)
write.csv(top_tss_df, "~/ELNID001000000000486355/DiffAnalysis/CUT_TAG/top50_accessible_TSS.csv", row.names = FALSE)

# Calculate average signals across genes for each bin
profile_atac_high <- colMeans(avg_mat_high, na.rm = TRUE)
profile_cuttag_high <- colMeans(mat_list_top_ct$High.Fatigue1, na.rm = TRUE)
profile_atac_low <- colMeans(avg_mat_med, na.rm = TRUE)
profile_cuttag_low <- colMeans(mat_list_top_ct$Medium.Fatigue1, na.rm = TRUE)

# Create a vector of positions relative to TSS, assuming 10 bp bins over ± 2000 bp
positions <- seq(-2000, 2000, length.out = 667)  # Midpoint of bins

# Plot average signals overlaid
png("~/ELNID001000000000486355/DiffAnalysis/TSSarea_ATAC_CUTTAG_overlay.png",
    height = 4, width = 6, units = "in", res = 300)
plot(positions, profile_atac_high, type = "l", col = "#ca0020", lwd = 2,
     xlab = "Position relative to TSS (bp)", ylab = "Average signal",
     main = "Average ATAC and CUT&TAG signals around TSS")
lines(positions, profile_cuttag_high, col = "#0571b0", lwd = 2)
lines(positions, profile_atac_low, col = "#f4a582", lwd = 2)
lines(positions, profile_cuttag_low, col = "#92c5de", lwd = 2)
legend("topright", legend = c("ATAC - High Fatigue", "ATAC - Low Fatigue", "CUT&TAG - High Fatigue", "CUT&TAG - Low Fatigue"),
       col = c("#ca0020", "#f4a582", "#0571b0", "#92c5de"), lwd = 2)
dev.off()

# Plot ATAC-seq average signals 
png("~/ELNID001000000000486355/DiffAnalysis/TSSarea_ATAConly_overlay.png",
    height = 4, width = 6, units = "in", res = 300)
plot(positions, profile_atac_high, type = "l", col = "#ca0020", lwd = 2,
     xlab = "Position relative to TSS (bp)", ylab = "Average signal",
     main = "Average ATAC signal around TSS")
lines(positions, profile_atac_low, col = "#f4a582", lwd = 2)
legend("topright", legend = c("ATAC - High Fatigue", "ATAC - Low Fatigue"),
       col = c("#ca0020", "#f4a582"), lwd = 2)
dev.off()

# Plot CUT&TAG average signals
png("~/ELNID001000000000486355/DiffAnalysis/TSSarea_CUTTAGonly_overlay.png",
    height = 4, width = 9, units = "in", res = 300)
ylim_shared <- range(c(profile_cuttag_high, profile_cuttag_low), na.rm = TRUE)
plot(positions, profile_cuttag_high, type = "l", col = "#0571b0", lwd = 2,
     xlab = "Position relative to TSS (bp)", ylab = "Average signal",
     main = "Average CUT&TAG signals around TSS",
     ylim = ylim_shared)
lines(positions, profile_cuttag_low, col = "#92c5de", lwd = 2)
legend("topright", legend = c("CUT&TAG - High Fatigue", "CUT&TAG - Low Fatigue"),
       col = c("#0571b0", "#92c5de"), lwd = 2)
dev.off()
