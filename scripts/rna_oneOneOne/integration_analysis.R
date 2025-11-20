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
volcano_plot <- function(df1, df2, df1_name, df2_name, axis_col, ttl, ct_genes){
  df <- create_volcano_df(df1, df2, df1_name, df2_name)
  df <- df %>%
    dplyr::mutate(Significance = factor(Significance, 
                                        levels = c("NS", "RNA", "ATAC", "ATAC & RNA"))) %>%
    dplyr::arrange(Significance)  # ensures "ATAC & RNA" is drawn last
  df[[df1_name]] <- gsub("GOBP_", "",df[[df1_name]])

  ggplot(df, aes(x=!!sym(paste0(axis_col, ".x")), y=!!sym(paste0(axis_col, ".y")), color = Significance)) + 
    geom_point()+ 
    geom_label(
      data=df %>% dplyr::filter(SYMBOL %in% ct_genes),
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
rnaseq_diff_results_file <- "~/ELNID001000000000486355/RNAseq/Diff_analysis/diff.Medium.fatigue_vs_High.fatigue/edgeR_results_allgenes.csv"
atacseq_diff_results_file <- "~/ELNID001000000000486355/DiffAnalysis/ATAC/annotatedDifferentialPeaks_noSexChr.csv"
ct_diff_results_file <- "/net/beegfs/groups/tytgat/ELNID001000000000486355/DiffAnalysis/CUT_TAG/annotatedDifferentialPeaks.csv"
rnaseq_fgsea_results_file <- "~/ELNID001000000000486355/DiffAnalysis/RNA_oneOnOne/enrichedPathwaysHighFatigue.csv"
atacseq_fgsea_results_file <- "~/ELNID001000000000486355/DiffAnalysis/ATAC/enrichedPathwaysHighFatigue_noSexChr.csv"
gene_scatter_file <- "~/ELNID001000000000486355/DiffAnalysis/RNA_oneOnOne//DEGscatterATAC_RNA_CUT_TAGgenes.png"
pathway_scatter_file <- "~/ELNID001000000000486355/DiffAnalysis/RNA_oneOnOne/pathwaysScatterATAC_RNA.png"
pathway_enrichment_point_plot <- "~/ELNID001000000000486355/DiffAnalysis/RNA_oneOnOne/pathwayEnrichment_ATACvsRNA_padj.png"

#combine ATAC-seq with RNA-seq data
#differentially expressed genes
rna_res <- read_csv(rnaseq_diff_results_file) 
rna_res <- rna_res %>% dplyr::filter(!is.na(gene)) %>% dplyr::filter(!is.na(pvalue))

atac_res <- read_csv(atacseq_diff_results_file) 
atac_res <- atac_res %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::filter(!is.na(pvalue)) %>% filter(grepl("Promoter", annotation)) %>%
  group_by(SYMBOL) %>%
  # Select the peak with the smallest p-value per gene
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() 

ct_res <- read_csv(ct_diff_results_file)
ct_genes <- ct_res$annotated_gene[ct_res$annotation != "Distal Intergenic"]

volcano_plot(atac_res, rna_res, "SYMBOL", "gene", "log2FoldChange", ttl = "ATAC-seq Vs RNA-seq DEGs", ct_genes)
ggsave(gene_scatter_file, plot = last_plot(), width = 8, height = 6, units = "in")

#differentially enriched pathways
atac_fgsea <- read_csv(atacseq_fgsea_results_file)
rna_fgsea <- read_csv(rnaseq_fgsea_results_file)
volcano_plot(atac_fgsea, rna_fgsea, "pathway", "pathway", "NES", ttl = "ATAC-seq Vs RNA-seq pathways")
ggsave(pathway_scatter_file, plot = last_plot(), width = 18, height = 10, units = "in")

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
  
ggsave(pathway_enrichment_point_plot, plot = p, dpi = 400, width = 8, height = 6)

# first integrate genes and then run fgsea:
# upregulated in both ATAC and RNA seq
atac_up_genes <- atac_res %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::filter(pvalue < 0.05)
rna_up_genes <- rna_res %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::filter(pvalue < 0.05)

double_up_genes <- unique(atac_up_genes$SYMBOL[atac_up_genes$SYMBOL %in% rna_up_genes$gene])
# downregulated in RNA and upregulated in ATAC
rna_down_genes <- rna_res %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::filter(pvalue < 0.05)

upATACdownRNA_genes <- unique(atac_up_genes$SYMBOL[atac_up_genes$SYMBOL %in% rna_down_genes$gene])
# downregulated in both ATAC and RNA seq
atac_down_genes <- atac_res %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::filter(pvalue < 0.05)

double_down_genes <- unique(atac_down_genes$SYMBOL[atac_down_genes$SYMBOL %in% rna_down_genes$gene])
# upregulated in RNA and downregulated in ATAC
downATACupRNA_genes <- unique(atac_down_genes$SYMBOL[atac_down_genes$SYMBOL %in% rna_up_genes$gene])


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


#show downregulated genes in ATAC
rna_res_downregulated <- rna_res %>% 
  dplyr::arrange(log2FoldChange) %>%
  dplyr::select(gene)
top_downregulated <- rna_res_downregulated$gene[1:10]
source("/net/beegfs/groups/tytgat/ELNID001000000000486355/peak_and_expression_correlation/peaksExpressionAnalysis.R")  

for(gene in top_downregulated){
  plot_gene_track(
    gene,
    atac_hf1_bw, atac_hf2_bw, atac_hf3_bw,
    atac_lf1_bw, atac_lf2_bw, atac_lf3_bw,
    atac_hf1_bed, atac_hf2_bed, atac_hf3_bed,
    atac_lf1_bed, atac_lf2_bed, atac_lf3_bed,
    ser_hf_bigwig_file, ser_mf_bigwig_file,
    ser_hf_bed_file, ser_mf_bed_file,
    rna_hf1_bw, rna_hf2_bw, rna_hf3_bw,
    rna_lf1_bw, rna_lf2_bw, rna_lf3_bw,
    "upregulated", -1.308151e-06,
    genome = "hg38",
    edb = EnsDb.Hsapiens.v86,
    window = 50000, 
    ct_y_top = 0.03,
    atac_y_top = 1.5,
    rna_y_top = 6
  )
  
}

plot_gene_track(
  "BCORP1",
  atac_hf1_bw, atac_hf2_bw, atac_hf3_bw,
  atac_lf1_bw, atac_lf2_bw, atac_lf3_bw,
  atac_hf1_bed, atac_hf2_bed, atac_hf3_bed,
  atac_lf1_bed, atac_lf2_bed, atac_lf3_bed,
  ser_hf_bigwig_file, ser_mf_bigwig_file,
  ser_hf_bed_file, ser_mf_bed_file,
  rna_hf1_bw, 
  rna_lf1_bw, 
  "upregulated", -1.308151e-06,
  genome = "hg38",
  edb = EnsDb.Hsapiens.v86,
  window = 5000, 
  ct_y_top = 0.015,
  atac_y_top = 0.08,
  rna_y_top = 20
)

