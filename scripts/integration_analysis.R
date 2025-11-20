library(dplyr)
library(GenomicRanges)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(rlang)
library(tidyverse)
library(clusterProfiler)
library("org.Hs.eg.db")


create_volcano_df <- function(first, second, first_name, second_name){
  merged_df <- merge(first, second, by.x = first_name, by.y = second_name)
  if ("pvalue.x" %in% colnames(merged_df)) {
    merged_df$Significance <- ifelse(merged_df$pvalue.x < 0.05,
                                     ifelse(merged_df$pvalue.y < 0.05, "ATAC & RNA", "ATAC"),
                                     ifelse(merged_df$pvalue.y < 0.05, "RNA", "NS" ))
  } else {
    merged_df$Significance <- ifelse(merged_df$pval.x < 0.05,
                                     ifelse(merged_df$pval.y < 0.05, "ATAC & RNA", "ATAC"),
                                     ifelse(merged_df$pval.y < 0.05, "RNA", "NS" ))
  }
  merged_df$Significance <- ifelse(is.na(merged_df$Significance), "NS", merged_df$Significance)
  merged_df$Significance <- factor(merged_df$Significance, levels = c("ATAC & RNA",
                                                                      "ATAC",  "RNA",
                                                                      "NS"))
  return(merged_df)
}
volcano_plot <- function(df1, df2, df1_name, df2_name, axis_col, ttl){
  df <- create_volcano_df(df1, df2, df1_name, df2_name)
  df <- df %>%
    dplyr::mutate(Significance = factor(Significance, 
                                        levels = c("NS", "RNA", "ATAC", "ATAC & RNA"))) %>%
    dplyr::arrange(Significance)  # ensures "ATAC & RNA" is drawn last
  df[[df1_name]] <- gsub("GOBP_", "",df[[df1_name]])

  ggplot(df, aes(x=!!sym(paste0(axis_col, ".x")), y=!!sym(paste0(axis_col, ".y")), color = Significance)) + 
    geom_point()+ 
    geom_label(
      data=df %>% dplyr::filter(Significance != "NS"),
      aes(label=!!sym(df1_name)),size = 3, show.legend = F) +
    geom_vline(xintercept = 0, color="black", linetype="dashed", alpha = 0.75) + geom_hline(yintercept = 0, color="black", linetype="dashed", alpha = 0.75) +
    labs(title = ttl, x = "ATAC-seq", y = "RNA-seq") +
    guides(color=guide_legend(title="P.Value < 0.05"))+ 
    xlim(-2, 2)+
    theme_bw() +
    scale_color_manual(values = c("ATAC & RNA" = "#e31a1c", 
                                  "ATAC" = "#fb9a99", 
                                  "RNA" = "#fdbf6f", 
                                  "NS" = "#e0e0e0"))
}

RNA_RES_PATH <- "~/ELNID001000000000486355/DiffAnalysis/RNA/"
ATAC_RES_PATH <- "~/ELNID001000000000486355/DiffAnalysis/ATAC/"


#combine ATAC-seq with RNA-seq data
#atac data = geneRnak
#rna data = DE results
rna_res <- read_csv(paste0(RNA_RES_PATH, "differentiallAnalysis_RNAseq_HighVSLow.csv"))
rna_res <- rna_res %>% dplyr::filter(!is.na(hgnc_symbol)) %>% dplyr::filter(!is.na(pvalue))

atac_res <- read_csv(paste0(ATAC_RES_PATH, "annotatedDifferentialPeaks_noSexChr.csv"))
atac_res <- atac_res %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::filter(!is.na(pvalue)) %>% filter(grepl("Promoter", annotation)) %>%
  group_by(SYMBOL) %>%
  # Select the peak with the smallest p-value per gene
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() 

volcano_plot(atac_res, rna_res, "SYMBOL", "hgnc_symbol", "log2FoldChange", ttl = "ATAC-seq Vs RNA-seq DEGs")
ggsave(paste0(RNA_RES_PATH, "DEGscatterATAC_RNA.png"), plot = last_plot(), width = 8, height = 6, units = "in")

atac_fgsea <- read_csv(paste0(ATAC_RES_PATH, "enrichedPathwaysHighFatigue_noSexChr.csv"))
rna_fgsea <- read_csv(paste0(RNA_RES_PATH, "enrichedPathwaysHighFatigue.csv"))
volcano_plot(atac_fgsea, rna_fgsea, "pathway", "pathway", "NES", ttl = "ATAC-seq Vs RNA-seq pathways")
ggsave(paste0(RNA_RES_PATH, "pathwaysScatterATAC_RNA.png"), plot = last_plot(), width = 18, height = 10, units = "in")

merged_fgsea <- merge(atac_fgsea, rna_fgsea, by = "pathway")
merged_long <- merged_fgsea  %>%
  dplyr::filter(padj.x < 0.05 & padj.y < 0.05)%>%
  dplyr::select(pathway,
         starts_with("pval"), starts_with("padj"),
         starts_with("log2err"), starts_with("ES"),
         starts_with("NES"), starts_with("size"), starts_with("leadingEdge")) %>%
  pivot_longer(
    cols = -pathway,
    names_to = c(".value", "omics_layer"),
    names_pattern = "(.*)\\.(x|y)"
  ) %>%
  mutate(omics_layer = gsub("x", "ATAC-seq", omics_layer))%>%
  mutate(omics_layer = gsub("y", "RNA-seq", omics_layer)) %>%
  mutate(pathway = gsub("GOBP_", "",pathway)) %>%
  mutate(omics_layer = factor(omics_layer))

p <- ggplot(merged_long , aes(x=omics_layer, y=pathway, color = NES, size = size)) +
  geom_point() +
  scale_color_gradient2(low = "red", mid = "white", high = "blue") +
  theme_bw() +
  ggtitle("Pathways Enriched in High Fatigue: ATAC-seq Vs. RNA-seq")
  
ggsave(paste0(RNA_RES_PATH, "pathwayEnrichment_ATACvsRNA_padj.png"), plot = p, dpi = 400, width = 8, height = 6)

# first integrate genes and then run fgsea:
# upregulated in both ATAC and RNA seq
atac_up_genes <- atac_res %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::filter(pvalue < 0.05)
rna_up_genes <- rna_res %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::filter(pvalue < 0.05)

double_up_genes <- unique(atac_up_genes$SYMBOL[atac_up_genes$SYMBOL %in% rna_up_genes$hgnc_symbol])
# downregulated in RNA and upregulated in ATAC
rna_down_genes <- rna_res %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::filter(pvalue < 0.05)

upATACdownRNA_genes <- unique(atac_up_genes$SYMBOL[atac_up_genes$SYMBOL %in% rna_down_genes$hgnc_symbol])
# downregulated in both ATAC and RNA seq
atac_down_genes <- atac_res %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::filter(pvalue < 0.05)

double_down_genes <- unique(atac_down_genes$SYMBOL[atac_down_genes$SYMBOL %in% rna_down_genes$hgnc_symbol])
# upregulated in RNA and downregulated in ATAC
downATACupRNA_genes <- unique(atac_down_genes$SYMBOL[atac_down_genes$SYMBOL %in% rna_up_genes$hgnc_symbol])


enrich_go <- enrichGO(
  gene = double_up_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",        # or "ENTREZID"
  ont = "BP",                # BP, CC, or MF
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

dotplot(enrich_go, title = "Gene Enrichment - upregulated in both ATAC-seq and RNA-seq")

enrich_go <- enrichGO(
  gene = double_down_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",        # or "ENTREZID"
  ont = "BP",                # BP, CC, or MF
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

dotplot(enrich_go, title = "Gene Enrichment - downregulated in both ATAC-seq and RNA-seq")

enrich_go <- enrichGO(
  gene = upATACdownRNA_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",        # or "ENTREZID"
  ont = "BP",                # BP, CC, or MF
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

dotplot(enrich_go, title = "Gene Enrichment - upregulated in ATAC-seq and downregulated in RNA-seq")

enrich_go <- enrichGO(
  gene = downATACupRNA_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",        # or "ENTREZID"
  ont = "BP",                # BP, CC, or MF
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

dotplot(enrich_go, title = "Gene Enrichment - downregulated in ATAC-seq and upregulated in RNA-seq")
