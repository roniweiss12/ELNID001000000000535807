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

# Flatten all bw files into a single vector for multiBigwigSummary
all_bw_files <- unlist(lapply(marks, function(mark) {
  c(
    file.path(base_dir, samples$high, "predicted_bw", paste0("out-R-", mark, ".bw")),
    file.path(base_dir, samples$medium, "predicted_bw", paste0("out-R-", mark, ".bw"))
  )
}))

# Verify files exist before running
missing <- all_bw_files[!file.exists(all_bw_files)]
if (length(missing) > 0) {
  cat("Missing files:\n", paste(missing, collapse = "\n"))
} else {
  cat("All", length(all_bw_files), "files found\n")
}

# Build and run the command for each mark separately (easier to manage)
lapply(marks, function(mark) {
  bw_high   <- file.path(base_dir, samples$high,   "predicted_bw", paste0("out-R-", mark, ".bw"))
  bw_medium <- file.path(base_dir, samples$medium, "predicted_bw", paste0("out-R-", mark, ".bw"))
  
  bw_string <- paste(c(bw_high, bw_medium), collapse = " ")
  
  labels <- paste(
    c(paste0("high_", 1:3), paste0("medium_", 1:3)),
    collapse = " "
  )
  conda_path <- "/appdata/users/P094415/miniconda3/bin/conda"
  
  cmd <- paste(
    conda_path, "run -n deept multiBigwigSummary BED-file",
    "--bwfiles", bw_string,
    "--BED ~/ELNID001000000000535807/analysis/tss_windows.bed",
    "--outFileName", paste0("~/ELNID001000000000535807/analysis/predicted_histone_modifications/", mark, "_summary.npz"),
    "--outRawCounts", paste0("~/ELNID001000000000535807/analysis/predicted_histone_modifications/", mark, "_counts.tab"),
    "--labels", labels,
    "-p 6"
  )
  
  cat("Running:", mark, "\n")
  system(cmd)
})


# Directory where your output files are
out_dir <- "~/ELNID001000000000535807/analysis/predicted_histone_modifications"  # adjust if different

# Read all count tables into a named list
counts_list <- lapply(marks, function(mark) {
  f <- file.path(out_dir, paste0(mark, "_counts.tab"))
  df <- read.table(f, header = FALSE, sep = "\t", comment.char = "#")
  colnames(df) <- c("chr", "start", "end", "high_1", "high_2","high_3","medium_1","medium_2","medium_3")
  df
})
names(counts_list) <- marks

# Check one looks right
head(counts_list[["H3K27ac"]])
dim(counts_list[["H3K27ac"]])



library(DESeq2)

# Run differential binding analysis for each mark
results_list <- lapply(marks, function(mark) {
  df <- counts_list[[mark]]
  
  # Extract count matrix (drop chr/start/end)
  mat <- df[, c("high_1", "high_2", "high_3", "medium_1", "medium_2", "medium_3")]
  
  # DESeq2 expects integers - scale and round
  mat_int <- round(mat * 1e6)
  
  # Remove rows with all zeros
  mat_int <- mat_int[rowSums(mat_int) > 0, ]
  
  # Sample metadata
  coldata <- data.frame(
    condition = factor(c("high", "high", "high", "medium", "medium", "medium"),
                       levels = c("medium", "high")),
    row.names = colnames(mat_int)
  )
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = mat_int,
    colData   = coldata,
    design    = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results: positive LFC = higher in high fatigue
  res <- results(dds, contrast = c("condition", "high", "medium"))
  res_df <- as.data.frame(res)
  
  # Add coordinates back
  coords <- df[rowSums(round(df[, 4:9] * 1e6)) > 0, c("chr", "start", "end")]
  res_df <- cbind(coords, res_df)
  
  res_df
})
names(results_list) <- marks

# Quick summary of significant regions per mark
lapply(marks, function(mark) {
  res <- results_list[[mark]]
  sig <- sum(res$padj < 0.05 & !is.na(res$padj))
  cat(mark, "- significant regions:", sig, "\n")
})



library(ggplot2)
library(dplyr)

# 1. Volcano plot for each mark
lapply(marks, function(mark) {
  res <- results_list[[mark]] %>%
    mutate(
      significant = !is.na(padj) & padj < 0.05,
      direction   = case_when(
        significant & log2FoldChange > 0 ~ "Higher in High",
        significant & log2FoldChange < 0 ~ "Higher in Medium",
        TRUE ~ "NS"
      )
    )
  
  ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = direction)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c(
      "Higher in High"   = "#E41A1C",
      "Higher in Medium" = "#377EB8",
      "NS"               = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(title = paste(mark, "- High vs Medium Fatigue"),
         x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_bw()
  
  ggsave(paste0(out_dir, "/volcano_", mark, ".pdf"),
         width = 6, height = 5)
})

# 2. Summary table across all marks
summary_df <- do.call(rbind, lapply(marks, function(mark) {
  res <- results_list[[mark]]
  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  data.frame(
    mark          = mark,
    n_significant = nrow(sig),
    n_higher_high   = sum(sig$log2FoldChange > 0),
    n_higher_medium = sum(sig$log2FoldChange < 0),
    median_lfc    = median(res$log2FoldChange, na.rm = TRUE)
  )
}))
print(summary_df)

# 3. Look at the significant regions in detail for top marks
for (mark in c("H3K4me2", "H3K4me3", "H3K9me3", "H4K20me1")) {
  cat("\n---", mark, "---\n")
  res <- results_list[[mark]]
  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  sig <- sig[order(sig$padj), ]
  print(head(sig[, c("chr", "start", "end", "log2FoldChange", "padj")], 10))
}


# Re-run summary with relaxed thresholds
summary_relaxed <- do.call(rbind, lapply(marks, function(mark) {
  res <- results_list[[mark]]
  
  # Three thresholds
  sig_strict  <- res[!is.na(res$padj)   & res$padj   < 0.05, ]
  sig_relaxed <- res[!is.na(res$padj)   & res$padj   < 0.10, ]
  sig_pval    <- res[!is.na(res$pvalue) & res$pvalue < 0.05, ]
  
  data.frame(
    mark              = mark,
    padj_0.05         = nrow(sig_strict),
    padj_0.10         = nrow(sig_relaxed),
    pval_0.05         = nrow(sig_pval),
    n_higher_high     = sum(sig_pval$log2FoldChange > 0),
    n_higher_medium   = sum(sig_pval$log2FoldChange < 0),
    median_lfc        = round(median(res$log2FoldChange, na.rm = TRUE), 3)
  )
}))

print(summary_relaxed)
