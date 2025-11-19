library(dplyr)
library(DESeq2)
library(readr)
library(biomaRt)
library(ggrepel)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(enrichplot)

volcano_plot <- function(df){
  df$StatSig <- "NS"
  # if log2Foldchange > 0 and pvalue < 0.05, set as "UP"
  df$StatSig[df$log2FoldChange > 1 & df$pvalue < 0.05] <- "Upregulated"
  # if log2Foldchange < 0 and pvalue < 0.05, set as "DOWN"
  df$StatSig[df$log2FoldChange < -1 & df$pvalue < 0.05] <- "Downregulated"
  df$StatSig <- factor(df$StatSig, levels = c("Upregulated", "Downregulated", "NS"))
  df$delabel <- ifelse(df$StatSig != "NS", df$hgnc_symbol, NA)
  df$shape <- ifelse(!is.na(df$padj) & df$padj < 0.05 , TRUE, FALSE)
  
  ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), color = StatSig, label = delabel, shape = shape)) +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = 1, col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = -1, col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("Upregulated" = "#00AFBB", "Downregulated" = "#bb0c00", "NS" = "grey"), # to set the colours of our variable
                       labels = c("Upregulated", "Downregulated", "NS")) + # to set the labels 
    theme_bw() +
    labs(color = 'Expression Change', shape = "padj < 0.05",
         x = "log2FoldChange", y = expression("-log"[10]*"P.Value")) +
    ggtitle("Volcano plot - RNA-seq") + # Plot title
    geom_text_repel(max.overlaps = Inf) +
    xlim(-7.5,10)
  
  ggsave(filename = plot_path, plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)
  
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
count_table_to_heatmap <- function(count_table, G_list, grouping_vector, plot_path){
  # Assume `mat` is a peaks x samples matrix
  # Calculate variance across samples
  peak_var <- rowVars(as.matrix(count_table))
  
  # Select top 1000–5000 most variable peaks
  top_peaks <- order(peak_var, decreasing = TRUE)[1:1000]
  mat_top <- count_table[top_peaks, ]
  mat_top <- mat_top %>% 
    tibble::rownames_to_column(var = "ensembl_gene_id") %>%
    left_join(G_list, by = "ensembl_gene_id") %>% 
    tibble::column_to_rownames(var = "hgnc_symbol") %>%
    dplyr::select(-chromosome_name, -ensembl_gene_id)
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

RESULTS_PATH_RNA <- "~/ELNID001000000000486355/DiffAnalysis/RNA/"
counts_table_file <- "~/ELNID001000000000466898/bulkRNAseq_pipeline/output/rnaseq/counts/counts.txt"
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
metadata$patientID <- factor(metadata$patientID)

ct <- read.table(counts_table_file, sep = "\t", header = TRUE)

# Vector of CD samples
oredered_samples <- c("CD14", "CD20", "CD26", "CD28", "CD3", "CD6")
# Build regex with word boundaries
samples_regex <- paste0("\\b(", paste(oredered_samples, collapse="|"), ")\\b")
# Filter columns
ct_filtered <- ct[, !grepl("output", colnames(ct)) | grepl(samples_regex, colnames(ct))]
colnames(ct_filtered)[7:12] <- oredered_samples
row.names(ct_filtered) <- ct_filtered$Geneid
ct_filtered[1:6] <- NULL
ct_filtered <- ct_filtered[,metadata$patientID]
colnames(ct_filtered) <- row.names(metadata)

set.seed(123)
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = ct_filtered,
                              colData = metadata,
                              design = ~ group)
# Filter lowly expressed genes (sum counts >= 20)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
# remove sex chromosomes
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
rownames(dds) <- sub("\\..*", "", rownames(dds))
genes <- rownames(dds)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "chromosome_name"),values=genes,mart= mart)
autosomal_genes <- G_list$ensembl_gene_id[!G_list$chromosome_name %in% c("X", "Y", "chrX", "chrY")]
dds <- dds[rownames(dds) %in% autosomal_genes, ]
# Run DESeq2
dds <- DESeq(dds)  # default fitType = "parametric"
# Variance stabilization for visualization
vsd <- vst(dds, blind = TRUE)
#res <- results(dds, name="group_High_fatigue_vs_Low_fatigue")
# LFC shrinkage
res <- lfcShrink(dds, coef = "group_High_fatigue_vs_Low_fatigue", type = "apeglm")

res <- as.data.frame(res) %>% arrange(padj)
res$ensembl_gene_id <- sub("\\..*", "", row.names(res))
res <- res %>% left_join(G_list, by = "ensembl_gene_id")

#write_csv(res, paste0(RESULTS_PATH_RNA, "differentiallAnalysis_RNAseq_HighVSLow.csv"))
res <- read_csv(paste0(RESULTS_PATH_RNA, "differentiallAnalysis_RNAseq_HighVSLow.csv"))
res <- res %>% dplyr::filter(!is.na(hgnc_symbol))

# clustering and heatmap
plot_pca(vsd, "group", paste0(RESULTS_PATH_RNA, "pca.png"))

count_table <- data.frame(assay(dds))
count_table_to_heatmap(count_table, G_list, metadata$group, paste0(RESULTS_PATH_RNA, "heatmap.png"))

count_table$ensembl_gene_id <- sub("\\..*", "", row.names(count_table))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(count_table)
genes_clean <- sub("\\..*", "", genes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "entrezgene_id", "hgnc_symbol", "description"),values=genes_clean,mart= mart)
count_table_annotated <- count_table %>% left_join(G_list, by = "ensembl_gene_id")
#write_csv(count_table_annotated, paste0(RESULTS_PATH_RNA, "countTable.csv"))
count_table <- read_csv(paste0(RESULTS_PATH_RNA, "countTable.csv"))

volcano_plot(res, paste0(RESULTS_PATH_RNA, "volcano_plot.png"))

geneRank <- res %>%
  filter(!is.na(hgnc_symbol)) %>%
  group_by(hgnc_symbol) %>%
  # Select the peak with the smallest p-value per gene
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  # calculate a combined score
  mutate(score = -log10(pvalue) * sign(log2FoldChange)) %>%
  ungroup() %>%
  # Rank genes by the combined score decreasingly
  arrange(desc(score))

rankings <-  geneRank$score
names(rankings) <- geneRank$hgnc_symbol
rankings <- sort(rankings, decreasing = TRUE)

msigdbr_df <- msigdbr(species = "Homo sapiens")
msigdbr_gobp <- msigdbr_df[msigdbr_df$gs_subcollection == "GO:BP", ]
msigdbr_gobp_list <- split(x = msigdbr_gobp$gene_symbol, f = msigdbr_gobp$gs_name)

fgsea_res <- fgsea(pathways = msigdbr_gobp_list, 
                   stats    = rankings,
                   minSize  = 15,
                   maxSize  = 500)
fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) paste(unlist(x), collapse = ", "))
write.csv(fgsea_res, paste0(RESULTS_PATH_RNA, "enrichedPathwaysHighFatigue.csv"), row.names = FALSE, quote = FALSE)

##top pathways anbiased
sig_fgsea_res <- fgsea_res %>%
  dplyr::arrange(NES) %>%
  dplyr::filter(padj < 0.05 & abs(NES) > 1) %>%
  mutate(pathway = (gsub("GOBP_", "", pathway)))
# Reorder the pathway factor based on ES
sig_fgsea_res$pathway <- factor(sig_fgsea_res$pathway, levels = sig_fgsea_res$pathway)
pathway_dotplot(sig_fgsea_res, "Pathways Enriched in High Fatigue: RNA-seq", paste0(RESULTS_PATH_RNA, "pathwayEnrichment_topPathways.png"))
##pathways related to specific topics
pattern <- "OX|MITOC|METAB|GLYCOL|CIRC|HYDROGEN_PEROXIDE"
specific_fgsea_res <- sig_fgsea_res %>%
  dplyr::filter(grepl(pattern, pathway, ignore.case = TRUE))
# Reorder the pathway factor based on ES
specific_fgsea_res$pathway <- factor(specific_fgsea_res$pathway, levels = specific_fgsea_res$pathway)
pathway_dotplot(specific_fgsea_res, "Pathways Enriched in High Fatigue: RNA-seq", paste0(RESULTS_PATH_RNA, "pathwayEnrichment_selectedPathways.png"))
