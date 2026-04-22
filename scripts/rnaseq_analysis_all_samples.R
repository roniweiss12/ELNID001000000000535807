library(dplyr)
library(DESeq2)
library(readr)
library(biomaRt)
library(ggrepel)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(readxl)
library(vegan)

counts_table_file <- "~/ELNID001000000000466898/bulkRNAseq_pipeline/output/rnaseq/counts/counts.txt"
ct <- read.table(counts_table_file, sep = "\t", header = TRUE)

# Vector of CD samples

full_metadata <- read_excel("~/ELNID001000000000466898/metadata/20250220_allmetadata_FIBD22_final-labeled_DNA.xlsx")
cd_fatigued_samples <- full_metadata %>%
  dplyr::filter(Type == "CD", MFI_totalscore >= 40) %>%
  dplyr::filter((MFI_totalscore >= 80) | (MFI_totalscore < 60)) 
  #dplyr::filter(Patient_ID != "CD7")
#oredered_samples <- c("CD14", "CD20", "CD26", "CD28", "CD3", "CD6")

# Build regex with word boundaries
samples_regex <- paste0("\\b(", paste(cd_fatigued_samples$Patient_ID, collapse="|"), ")\\b")

# Filter columns
ct_filtered <- ct[, !grepl("output", colnames(ct)) | grepl(samples_regex, colnames(ct))]

colnames(ct_filtered)[7:24] <- cd_fatigued_samples$Patient_ID
row.names(ct_filtered) <- ct_filtered$Geneid
ct_filtered[1:6] <- NULL

cd_fatigued_samples$Seks <- as.factor(cd_fatigued_samples$Seks)
cd_fatigued_samples$group <- ifelse(cd_fatigued_samples$MFI_totalscore < 80, "Low_fatigue", "High_fatigue")
cd_fatigued_samples$group <- factor(cd_fatigued_samples$group, levels = c("Low_fatigue", "High_fatigue"))
cd_fatigued_samples$MFI_scaled <- scale(cd_fatigued_samples$MFI_totalscore)
cd_fatigued_samples$Patient_ID <- factor(cd_fatigued_samples$Patient_ID)
row.names(cd_fatigued_samples) <- cd_fatigued_samples$Patient_ID
ct_filtered <- ct_filtered[,cd_fatigued_samples$Patient_ID]
colnames(ct_filtered) <- row.names(cd_fatigued_samples)

ct_filtered <- ct_filtered[rowSums(ct_filtered) != 0, ]

set.seed(123)
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = ct_filtered,
                              colData = cd_fatigued_samples,
                              design = ~ group + Seks)

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
res <- results(dds, name="group_High_fatigue_vs_Low_fatigue")
# LFC shrinkage
#res <- lfcShrink(dds, coef = "group_High_fatigue_vs_Low_fatigue", type = "apeglm")


res <- as.data.frame(res) %>% arrange(padj)
res$ensembl_gene_id <- sub("\\..*", "", row.names(res))
res <- res %>% left_join(G_list, by = "ensembl_gene_id")

#write_csv(res, "~/ELNID001000000000486355/analysis/RNA_allSamples/differentiallAnalysis_RNAseq_HighVSLow.csv")
res <- read_csv("~/ELNID001000000000535807/analysis/RNA_allSamples/differentiallAnalysis_RNAseq_HighVSLow.csv")
res <- res %>% dplyr::filter(!is.na(hgnc_symbol))


# clustering and heatmap
# Extract the normalized matrix
mat <- assay(vsd)

# Transpose: samples as rows, peaks as columns
pca <- prcomp(t(mat), scale. = TRUE)

# Create a data frame with metadata
pcaData <- data.frame(pca$x, colData(vsd))
# PERMANOVA test
dist_mat <- dist(t(mat))
res_permanova <- adonis2(dist_mat ~ group, data = as.data.frame(colData(vsd)), permutations = 999)

# Plot
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  #stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.5) +
  #geom_text_repel(aes(label = rownames(pcaData)), size = 3)+
  scale_color_brewer(palette = "Set2") +
  labs(x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)"),
       title = "PCA - RNA-seq") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("~/ELNID001000000000535807/DiffAnalysis/RNA_allSamples/pca.png", pca_plot, width = 6, height = 5, units = "in")

count_table <- data.frame(assay(dds))
# Assume `mat` is a peaks x samples matrix
# Calculate variance across samples
peak_var <- rowVars(as.matrix(count_table))

# Select top 1000–5000 most variable peaks
top_peaks <- order(peak_var, decreasing = TRUE)[1:1000]
mat_top <- count_table[top_peaks, ]
mat_top <- mat_top[, c(cd_fatigued_samples$Patient_ID[cd_fatigued_samples$group == "High_fatigue"], cd_fatigued_samples$Patient_ID[cd_fatigued_samples$group == "Low_fatigue"])]
mat_top <- mat_top %>% 
  tibble::rownames_to_column(var = "ensembl_gene_id") %>%
  left_join(G_list, by = "ensembl_gene_id") %>% 
  tibble::column_to_rownames(var = "hgnc_symbol") %>%
  dplyr::select(-chromosome_name, -ensembl_gene_id)
# Create and store the heatmap object
colSide <- brewer.pal(3, "Set2")[unlist(cd_fatigued_samples[colnames(mat_top), "group"])]
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
png("~/ELNID001000000000535807/DiffAnalysis/RNA_allSamples/heatmap.png", width = 6, height = 8, units = "in", res = 300)
heatmap <- heatmap(as.matrix(mat_top), col = coul, ColSideColors = colSide, labRow = NA, Colv = NA)
dev.off()
# png("~/ELNID001000000000535807/DiffAnalysis/RNA_allSamples/heatmap_legend.png", width = 4, height = 4, units = "in", res = 300)
# plot.new() 
# legend(
#   "center",
#   legend = c("Low", "Medium", "High"),
#   fill = coul[c(1, round(length(coul)/2), length(coul))],
#   title = "Expression", bty = "n"
# )
# dev.off()

count_table$ensembl_gene_id <- sub("\\..*", "", row.names(count_table))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(count_table)
genes_clean <- sub("\\..*", "", genes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "entrezgene_id", "hgnc_symbol", "description"),values=genes_clean,mart= mart)
count_table_annotated <- count_table %>% left_join(G_list, by = "ensembl_gene_id")
#write_csv(count_table_annotated, "~/ELNID001000000000535807/analysis/RNA_allSamples/countTable.csv")
count_table <- read_csv("~/ELNID001000000000535807/analysis/RNA_allSamples/countTable.csv")

volcano_plot <- function(df){
  df$StatSig <- "NS"
  # if log2Foldchange > 0 and pvalue < 0.05, set as "UP"
  df$StatSig[df$log2FoldChange > 1 & df$pvalue < 0.05] <- "Upregulated"
  # if log2Foldchange < 0 and pvalue < 0.05, set as "DOWN"
  df$StatSig[df$log2FoldChange < -1 & df$pvalue < 0.05] <- "Downregulated"
  df$StatSig <- factor(df$StatSig, levels = c("Upregulated", "Downregulated", "NS"))
  # Order df by padj ascending, keep row indices of top 10
  top10_idx <- order(df$padj)[1:30]
  
  # Set all labels to NA first
  df$delabel <- NA
  
  # Assign SYMBOL as label only for top 10 with lowest padj and significant StatSig (not "NS")
  df$delabel[top10_idx] <- ifelse(df$StatSig[top10_idx] != "NS", df$hgnc_symbol[top10_idx], NA)
  ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), color = StatSig, label = delabel)) +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = 1, col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = -1, col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("Upregulated" = "#00AFBB", "Downregulated" = "#bb0c00", "NS" = "grey"), # to set the colours of our variable
                       labels = c("Upregulated", "Downregulated", "NS")) + # to set the labels 
    theme_bw() +
    labs(color = 'Expression Change',
         x = "log2FoldChange", y = expression("-log"[10]*"P.Value")) +
    ggtitle("Volcano plot - RNA-seq") + # Plot title
    geom_text_repel(max.overlaps = Inf) +
    xlim(-7.5,7.5)
  
  ggsave(filename = "~/ELNID001000000000535807/analysis/RNA_allSamples/volcano_plot.png", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)
  
}

volcano_plot(res)


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
rankings <- rankings[!is.na(names(rankings)) & names(rankings) != ""]
rankings <- rankings[!duplicated(names(rankings))]


msigdbr_df <- msigdbr(species = "Homo sapiens")
msigdbr_gobp <- msigdbr_df[msigdbr_df$gs_subcollection == "GO:BP", ]
msigdbr_gobp_list <- split(x = msigdbr_gobp$gene_symbol, f = msigdbr_gobp$gs_name)

fgsea_res <- fgsea(pathways = msigdbr_gobp_list, 
                   stats    = rankings,
                   minSize  = 15,
                   maxSize  = 500)
fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) paste(unlist(x), collapse = ", "))
write.csv(fgsea_res, "~/ELNID001000000000535807/analysis/RNA_allSamples/enrichedPathwaysHighFatigue.csv", row.names = FALSE, quote = FALSE)
pattern <- "OX|MITOC|METAB|GLYCOL|CIRC|HYDROGEN_PEROXIDE"

sig_fgsea_res <- fgsea_res %>%
  dplyr::arrange(NES) %>%
  dplyr::filter(padj < 0.01 & abs(NES) > 1) %>%
  mutate(pathway = (gsub("GOBP_", "", pathway)))#%>%
  #dplyr::filter(grepl(pattern, pathway, ignore.case = TRUE))
# Reorder the pathway factor based on ES
sig_fgsea_res$pathway <- factor(sig_fgsea_res$pathway, levels = sig_fgsea_res$pathway)

p <- ggplot(sig_fgsea_res, aes(x=NES, y=pathway, color = pval, size = size)) +
  geom_point() +
  #scale_color_gradient(low = "white", high = "blue") +
  theme_bw() +
  ggtitle("Pathways Enriched in High Fatigue: RNA-seq ")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 14, hjust = 0.5))
ggsave("~/ELNID001000000000535807/analysis/RNA_allSamples/pathwayEnrichment.png", plot = p, dpi = 400, width = 10, height = 10, units = "in")
