library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(plyranges)
library(apeglm)
library(pheatmap)
library(matrixStats)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(readr)
library(tidyverse)
library(biomaRt)



combine_count_tables <- function(first_path, second_path){
  ct1 <- read.table(first_path, sep = "\t", header = TRUE)
  ct2 <- read.table(second_path, sep = "\t", header = TRUE)
  
  # Convert each to GRanges
  gr1 <- makeGRangesFromDataFrame(ct1, keep.extra.columns = TRUE)
  gr2 <- makeGRangesFromDataFrame(ct2, keep.extra.columns = TRUE)
  
  # Combine all
  gr_all <- c(gr1, gr2)
  
  # Reduce with revmap to keep track of original intervals merged into each reduced interval
  gr_union <- GenomicRanges::reduce(gr_all, with.revmap=TRUE)
  
  # Aggregate counts over the merged intervals by summing counts of the original intervals
  # Access metadata columns using mcols()
  all_mcols <- as.data.frame(mcols(gr_all))
  agg_counts <- sapply(count_columns, function(col) {
    col_counts <- all_mcols[[col]]
    vapply(gr_union$revmap, function(ix) sum(col_counts[ix], na.rm=TRUE), numeric(1))
  })
  
  # Add aggregated counts as metadata to the reduced GRanges
  mcols(gr_union)[count_columns] <- agg_counts
  mcols(gr_union)$revmap <- NULL  # Remove revmap to clean up
  
  # gr_union now has merged intervals with summed counts from overlapping original intervals
  
  #rename sample columns, remove low fatigue and order columns.
  colnames(mcols(gr_union)) <- new_count_columns
  mcols(gr_union) <- mcols(gr_union)[ , -which(colnames(mcols(gr_union)) %in% cols_to_remove) ]
  mcols(gr_union) <- mcols(gr_union)[order(colnames(mcols(gr_union)))]
  
  full_count_table <- as.data.frame(gr_union)
  row.names(full_count_table) <- paste0(full_count_table$seqnames, "_", full_count_table$start, "_", full_count_table$end)
  full_count_table <- full_count_table %>% dplyr::select(-seqnames, -start, -end, -width, -strand)
  #rename medium fatigue as low
  colnames(full_count_table) <- c("High.fatigue.1",  "High.fatigue.2", "High.fatigue.3", "Low.fatigue.1", "Low.fatigue.2", "Low.fatigue.3")
  names(mcols(gr_union)) <- c("High.fatigue.1",  "High.fatigue.2", "High.fatigue.3", "Low.fatigue.1", "Low.fatigue.2", "Low.fatigue.3")
  return(list(ct = full_count_table, gr = gr_union))
}
filter_random_chr <- function(de_results){
  rn <- rownames(de_results)
  
  # Create a logical vector for rows with "chrUn" or "random"
  to_remove <- grepl("chrUn", rn) | grepl("random", rn)
  
  # Subset the data to exclude those rows
  filtered_res <- de_results[!to_remove, ]
  return(filtered_res)
}
run_fgsea <- function(anno_df, res_path){
  geneRank <- anno_df %>%
    filter(!is.na(SYMBOL)) %>%
    filter(grepl("Promoter", annotation)) %>%
    group_by(SYMBOL) %>%
    # Select the peak with the smallest p-value per gene
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    # Optionally, calculate a combined score
    mutate(score = -log10(pvalue) * sign(log2FoldChange)) %>%
    ungroup() %>%
    # Rank genes by the combined score decreasingly
    arrange(desc(score))

  rankings <-  geneRank$score
  names(rankings) <- geneRank$SYMBOL
  rankings <- sort(rankings, decreasing = TRUE)

  msigdbr_df <- msigdbr(species = "Homo sapiens")
  msigdbr_gobp <- msigdbr_df[msigdbr_df$gs_subcollection == "GO:BP", ]
  msigdbr_gobp_list <- split(x = msigdbr_gobp$gene_symbol, f = msigdbr_gobp$gs_name)
  
  fgsea_res <- fgsea(pathways = msigdbr_gobp_list, 
                     stats    = rankings,
                     minSize  = 15,
                     maxSize  = 500)
  fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) paste(unlist(x), collapse = ", "))
  write_csv(fgsea_res, res_path)

  sig_fgsea_res <- fgsea_res %>%
    dplyr::arrange(NES) %>%
    dplyr::filter(padj < 0.05 & abs(NES) > 1) %>%
    mutate(pathway = (gsub("GOBP_", "", pathway)))
  # Reorder the pathway factor based on ES
  sig_fgsea_res$pathway <- factor(sig_fgsea_res$pathway, levels = sig_fgsea_res$pathway)
  
  return(sig_fgsea_res)
  
}
plot_pca <- function(vsd_counts, color_by, plot_path){
  mat <- assay(vsd_counts)
  
  # Transpose: samples as rows, peaks as columns
  pca <- prcomp(t(mat), scale. = TRUE)
  
  # Create a data frame with metadata
  pcaData <- data.frame(pca$x, colData(vsd_counts))
  
  # Plot
  pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = !!sym(color_by))) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = rownames(pcaData)), size = 3)+
    scale_color_brewer(palette = "Set2") +
    labs(x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
         y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)"),
         title = "PCA - ATAC-seq") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggsave(plot_path, pca_plot, width = 6, height = 5, units = "in")
}
add_geneHancer_annotation <- function(gr, geneHancer_file){
  geneHancer <- read.table(geneHancer_file, sep = "\t")
  
  # Parse the 'V9' column into a long format
  geneHancer_long <- geneHancer %>%
    mutate(V9 = str_replace_all(V9, "^\\s+|\\s+$", "")) %>%  # trim
    mutate(genehancer_id = str_extract(V9, "(?<=genehancer_id=)[^;]+")) %>% # extract ID
    mutate(connected_info = str_extract_all(V9, "connected_gene=[^;]+;score=[^;]+")) %>% # extract gene-score pairs
    dplyr::select(-V9) %>%
    unnest(connected_info) %>%  # expand one row per connected gene
    separate(connected_info, into = c("connected_gene", "score"), sep = ";") %>%
    mutate(
      connected_gene = str_remove(connected_gene, "connected_gene="),
      score = as.numeric(str_remove(score, "score="))
    ) %>%
    dplyr::filter(!grepl("^(piR-|LOC|ENSG0|MIR|lnc-|HSALNG)", connected_gene))
  
  enhancers_gr <- GRanges(
    seqnames = geneHancer_long$V1,
    ranges = IRanges(start = geneHancer_long$V4, end = geneHancer_long$V5),
    ID = geneHancer_long$genehancer_id,
    connected_gene = geneHancer_long$connected_gene,
    score = geneHancer_long$score
  )
  seqlevelsStyle(gr) <- seqlevelsStyle(enhancers_gr)
  peaks_enh_hits <- findOverlaps(gr, enhancers_gr)
  
  # 7. Combine enh info with mapped genes
  peaks_with_genes <- data.frame(
    chr = seqnames(gr)[queryHits(peaks_enh_hits)],
    start = start(gr)[queryHits(peaks_enh_hits)],
    end = end(gr)[queryHits(peaks_enh_hits)],
    GHenhancer = enhancers_gr[subjectHits(peaks_enh_hits)]
  )
  # Collapse multiple genes per enh into one row
  peaks_genes_top <- peaks_with_genes %>%
    group_by(GHenhancer.ID) %>%
    slice_max(order_by = GHenhancer.score, n = 1, with_ties = FALSE) %>%  # keep highest-scoring gene per enhancer
    ungroup()
  
  # Join back to original DMPs
  peaks_annotated <- filtered_res %>%
    tibble::rownames_to_column("region") %>%
    separate(region, into = c("chr", "start", "end"), sep = "_", convert = TRUE) %>%
    left_join(peaks_genes_top, by = c("chr", "start", "end"))
  return(peaks_annotated)
  
}
add_pchic_annotation <- function(pchic_file, gr_counts){
  pchic <- read_csv(pchic_file)
  enhancers_gr <- GRanges(
    seqnames = pchic$EnhChr,
    ranges = IRanges(start = pchic$EnhStart, end = pchic$EnhEnd),
    ID = pchic$EnhID
  )
  seqlevelsStyle(gr_counts) <- seqlevelsStyle(enhancers_gr)
  atac_enh_hits <- findOverlaps(gr_counts, enhancers_gr)

  # 7. Combine DMP info with mapped genes
  atac_with_genes <- data.frame(
    chr = seqnames(gr_counts)[queryHits(atac_enh_hits)],
    start = start(gr_counts)[queryHits(atac_enh_hits)],
    end = end(gr_counts)[queryHits(atac_enh_hits)],
    enhancer = enhancers_gr[subjectHits(atac_enh_hits)]
  ) %>%
    left_join(pchic, by = join_by("enhancer.ID" == "EnhID"))

  # Collapse multiple genes per DMP into one row
  atac_genes_collapsed <- atac_with_genes %>%
    group_by(chr, start, end) %>%
    summarise(potentialEnhancers = paste(unique(HGNC), collapse = ";"),
              .groups = "drop")

  # Join back to original DMPs
  atac_annotated <- filtered_res %>%
    tibble::rownames_to_column("region") %>%
    separate(region, into = c("chr", "start", "end"), sep = "_", convert = TRUE) %>%
    left_join(atac_genes_collapsed, by = c("chr", "start", "end"))
  return(atac_annotated)
}
volcano_plot <- function(df, plot_path){
  df$StatSig <- "NS"
  # if log2Foldchange > 0 and pvalue < 0.05, set as "UP"
  df$StatSig[df$log2FoldChange > 1 & df$pvalue < 0.05] <- "Upregulated"
  # if log2Foldchange < 0 and pvalue < 0.05, set as "DOWN"
  df$StatSig[df$log2FoldChange < -1 & df$pvalue < 0.05] <- "Downregulated"
  df$StatSig <- factor(df$StatSig, levels = c("Upregulated", "Downregulated", "NS"))
  #df$delabel <- ifelse((df$StatSig != "NS") & (df$padj < 0.01), df$SYMBOL, NA)
  df$shape <- ifelse(!is.na(df$padj) & df$padj < 0.01 , TRUE, FALSE)
  # Order df by padj ascending, keep row indices of top 10
  top_idx <- order(df$padj)[1:10]
  # Set all labels to NA first
  df$delabel <- NA
  
  # Assign SYMBOL as label only for top 10 with lowest padj and significant StatSig (not "NS")
  df$delabel[top10_idx] <- ifelse(df$StatSig[top10_idx] != "NS", df$SYMBOL[top_idx], NA)
  
  ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), color = StatSig, label = delabel, shape = shape)) +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = 1, col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = -1, col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("Upregulated" = "#00AFBB", "Downregulated" = "#bb0c00", "NS" = "grey"), # to set the colours of our variable
                       labels = c("Upregulated", "Downregulated", "NS")) + # to set the labels 
    theme_bw() +
    labs(color = 'Expression Change', shape = "padj < 0.01",
         x = "log2FoldChange", y = expression("-log"[10]*"P.Value")) +
    ggtitle("Volcano plot - RNA-seq") + # Plot title
    geom_text_repel(max.overlaps = Inf) #+
  xlim(-12,12)
  
  ggsave(filename = plot_path, plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)
  
}
pathway_dotplot <- function(df, title, plot_path){
  p <- ggplot(df, aes(x=NES, y=pathway, color = pval, size = size)) +
    geom_point() +
    theme_bw() +
    ggtitle(title)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 14, hjust = 0.5))
  ggsave(plot_path, plot = p, dpi = 400, width = 12, height = 8)
}
count_table_to_heatmap <- function(count_table, grouping_vector, plot_path){
  # Assume `mat` is a peaks x samples matrix
  # Calculate variance across samples
  peak_var <- rowVars(as.matrix(count_table))
  
  # Select top 1000–5000 most variable peaks
  top_peaks <- order(peak_var, decreasing = TRUE)[1:1000]
  mat_top <- count_table[top_peaks, ]
  # Create and store the heatmap object
  colSide <- brewer.pal(3, "Set2")[grouping_vector]
  coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
  png(plot_path, width = 6, height = 8, units = "in", res = 300)
  heatmap <- heatmap(as.matrix(mat_top), col = coul, ColSideColors = colSide, labRow = NA, Colv = NA)
  dev.off()
  # png(paste0(RESULTS_PATH_RNA, "heatmap_legend.png"), width = 4, height = 4, units = "in", res = 300)
  # plot.new() 
  # legend(
  #   "center",
  #   legend = c("Low", "Medium", "High"),
  #   fill = coul[c(1, round(length(coul)/2), length(coul))],
  #   title = "Expression", bty = "n"
  # )
  # dev.off()
}
count_table1_file <- "~/ELNID001000000000486355/data/Q52867_ATACseq/PeakCalling/consensus/consensus_peaks.mLb.clN.featureCounts.txt"
count_table2_file <- "~/ELNID001000000000535807/data/ATACseq/PeakCalling/consensus/consensus_peaks.mLb.clN.featureCounts.txt"
geneHancer_file <- "~/ELNID001000000000535807/analysis/geneHancer_data/GeneHancer_v5.25.gff"
pchic_file <- "~/ELNID001000000000466898/output/epic/resources/pchic_annotated_grch38.csv"
RESULTS_PATH_ATAC <- "~/ELNID001000000000535807/analysis/ATAC/"
count_columns <- c("X03_0KXN_02UCAUMC_Low.fatigue.3_ATAC_i3.3_REP1.mLb.clN.sorted.bam", 
                   "X08_0KXQ_02UCAUMC_High.fatigue.2_ATAC_i6.522_REP1.mLb.clN.sorted.bam",
                   "X02_0KXM_02UCAUMC_Low.fatigue.2_ATAC_i2.2_REP1.mLb.clN.sorted.bam",
                   "X01_0KXL_02UCAUMC_Low.fatigue.1_ATAC_i1.1_REP1.mLb.clN.sorted.bam",      
                   "X09_0KXR_02UCAUMC_High.fatigue.3_ATAC_i7.15_REP1.mLb.clN.sorted.bam",
                   "X06_0KXP_02UCAUMC_Medium.fatigue.3_ATAC_i5.521_REP1.mLb.clN.sorted.bam", 
                   "X05_0KXO_02UCAUMC_Medium.fatigue.2_ATAC_hs_i4.4_REP1.mLb.clN.sorted.bam",
                   "X03_0KC9_02SLTytgat_High.fatigue.1_ATAC_hs_i29.48_REP1.mLb.clN.sorted.bam",
                   "X01_0KC8_02SLTytgat_Medium.fatigue.1_ATAC_hs_i28.47_REP1.mLb.clN.sorted.bam"
)
new_count_columns <- c("Low.fatigue.3", 
                       "High.fatigue.2",
                       "Low.fatigue.2",
                       "Low.fatigue.1",      
                       "High.fatigue.3",
                       "Medium.fatigue.3", 
                       "Medium.fatigue.2",
                       "High.fatigue.1",
                       "Medium.fatigue.1")
cols_to_remove <- c("Low.fatigue.3", 
                    "Low.fatigue.2",
                    "Low.fatigue.1")

metadata <- data.frame(patientID = c("CD28", "CD20", "CD3", "CD26", "CD14", "CD6"),
                       sex = c("female", "female", "female", "male", "female", "male"), 
                       MFI = c(80, 80, 86, 40, 43, 46), 
                       group = c("High_fatigue", "High_fatigue", "High_fatigue", "Low_fatigue", "Low_fatigue", "Low_fatigue"), 
                       batch = c(1,2,2,1,2,2),
                       row.names = c("High.fatigue.1", "High.fatigue.2", "High.fatigue.3", "Low.fatigue.1", "Low.fatigue.2", "Low.fatigue.3"))
metadata$sex <- as.factor(metadata$sex)
metadata$batch <- as.factor(metadata$batch)
metadata$group <- factor(metadata$group, levels = c("Low_fatigue", "High_fatigue"))
metadata$MFI_scaled <- scale(metadata$MFI)



set.seed(123)
#combine data from 2 experiments
counts_data <- combine_count_tables(count_table1_file, count_table2_file)

full_count_table <- counts_data[["ct"]]
gr_counts <- counts_data[["gr"]]

# differentiall analysis
dds <- DESeqDataSetFromMatrix(countData = full_count_table,
                              colData = metadata,
                              design = ~ group + batch)

# Increase filtering threshold to keep peaks with sum counts >= 20
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep, ]
# remove sex chromosomes
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "chromosome_name"),values=genes,mart= mart)
autosomal_genes <- G_list$ensembl_gene_id[!G_list$chromosome_name %in% c("X", "Y", "chrX", "chrY")]
dds <- dds[!grepl(c("chrX|chrY"), rownames(dds)), ]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
vsd <- vst(dds, blind = TRUE)

# Apply LFC shrinkage for stable effect size estimates
res <- lfcShrink(dds, coef = "group_High_fatigue_vs_Low_fatigue", type = "apeglm")
res <- as.data.frame(res) %>% arrange(padj)
# Subset the data to exclude those rows
filtered_res <- filter_random_chr(res)

#write.csv(filtered_res, paste0(RESULTS_PATH_ATAC, "differentiallAnalysis_ATACseq_HighVSLow_noSexChr.csv"), row.names = TRUE, quote = FALSE)
filtered_res <- read.table(paste0(RESULTS_PATH_ATAC, "differentiallAnalysis_ATACseq_HighVSLow.csv"), sep = ",", row.names = 1, header = TRUE)
# clustering and heatmap
plot_pca(vsd, "group", paste0(RESULTS_PATH_ATAC, "pca.png"))

count_table <- data.frame(assay(dds))
#write.csv(count_table, paste0(RESULTS_PATH_ATAC, "countTable.csv"), row.names = TRUE, quote = FALSE)
count_table <- read.table(paste0(RESULTS_PATH_ATAC, "countTable.csv"), sep = ",", row.names = 1, header = TRUE)
count_table_to_heatmap(count_table, metadata$group, paste0(RESULTS_PATH_ATAC, "heatmap.png"))

#obtain gene annotations from PCHiC/geneHancer
########################################
#peaks_pchic_annotated <- add_pchic_annotation(pchic_file, gr_counts)
#peaks_gh_annotated <- add_geneHancer_annotation(gr_counts, geneHancer_file)

#peak annotation - nearest TSS
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Annotate peaks
peak_annot <- annotatePeak(gr_counts, TxDb=txdb, annoDb="org.Hs.eg.db")
anno_df <- as.data.frame(peak_annot)
#write.table(anno_df, paste0(RESULTS_PATH_ATAC, "annotatedDifferentialPeaks.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
#merge with differential accessability data
filtered_res$rowname <- row.names(filtered_res)
anno_df <- anno_df %>%
  mutate(rowname = paste0(seqnames, "_", start, "_", end)) %>%
  right_join(filtered_res, by = "rowname")

#write_csv(anno_df, paste0(RESULTS_PATH_ATAC, "annotatedDifferentialPeaks_noSexChr.csv"))
anno_df <- read_csv(paste0(RESULTS_PATH_ATAC, "annotatedDifferentialPeaks.csv"))

geneRank <- anno_df %>%
  filter(!is.na(SYMBOL)) %>%
  filter(grepl("Promoter", annotation)) %>%
  group_by(SYMBOL) %>%
  # Select the peak with the smallest p-value per gene
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  # Optionally, calculate a combined score
  mutate(score = -log10(pvalue) * sign(log2FoldChange)) %>%
  ungroup() %>%
  # Rank genes by the combined score decreasingly
  arrange(desc(score))

volcano_plot(geneRank, paste0(RESULTS_PATH_ATAC, "volcano_plot_noSexChr.png"))

sig_fgsea_res <- run_fgsea(anno_df, paste0(RESULTS_PATH_ATAC, "enrichedPathwaysHighFatigue_noSexChr.csv"))
#show top unbiased pathways
sig_fgsea_res$pathway <- factor(sig_fgsea_res$pathway, levels = sig_fgsea_res$pathway)
pathway_dotplot(sig_fgsea_res %>%
                  arrange(padj) %>%
                  slice_tail(n = 40), "Pathways Enriched in High Fatigue: ATAC-seq", 
                paste0(RESULTS_PATH_ATAC, "pathwayEnrichment_top40unbiased.png"))
#show relevant pathways
pattern <- "OX|MITOC|METAB|GLYCOL|CIRC|HYDROGEN_PEROXIDE"
specific_fgsea_res <- sig_fgsea_res %>%
  dplyr::filter(grepl(pattern, pathway, ignore.case = TRUE))
# Reorder the pathway factor based on ES
specific_fgsea_res$pathway <- factor(specific_fgsea_res$pathway, levels = sig_fgsea_res$pathway)
pathway_dotplot(specific_fgsea_res, "Pathways Enriched in High Fatigue: ATAC-seq", 
                paste0(RESULTS_PATH_ATAC, "pathwayEnrichment_selectedPathways_noSexChr.png"))
