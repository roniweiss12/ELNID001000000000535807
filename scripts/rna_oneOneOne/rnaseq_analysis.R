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
  ggsave(plot_pathplot_path, plot = p, dpi = 400, width = 12, height = 8)
}
volcano_plot <- function(df, plot_path){
  df$StatSig <- "NS"
  # if log2Foldchange > 0 and pvalue < 0.05, set as "UP"
  df$StatSig[df$log2FoldChange > 1 & df$pvalue < 0.05] <- "Upregulated"
  # if log2Foldchange < 0 and pvalue < 0.05, set as "DOWN"
  df$StatSig[df$log2FoldChange < -1 & df$pvalue < 0.05] <- "Downregulated"
  df$StatSig <- factor(df$StatSig, levels = c("Upregulated", "Downregulated", "NS"))
  df$delabel <- ifelse((df$StatSig != "NS") & (df$FDR < 0.05), df$gene, NA)
  df$shape <- ifelse(!is.na(df$FDR) & df$FDR < 0.05 , TRUE, FALSE)
  
  ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), color = StatSig, label = delabel, shape = shape)) +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = 1, col = "gray", linetype = 'dashed') +
    geom_vline(xintercept = -1, col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("Upregulated" = "#00AFBB", "Downregulated" = "#bb0c00", "NS" = "grey"), # to set the colours of our variable
                       labels = c("Upregulated", "Downregulated", "NS")) + # to set the labels 
    theme_bw() +
    labs(color = 'Expression Change', shape = "FDR < 0.05",
         x = "log2FoldChange", y = expression("-log"[10]*"P.Value")) +
    ggtitle("Volcano plot - RNA-seq") + # Plot title
    geom_text_repel(max.overlaps = Inf) +
    xlim(-7.5,10)
  ggsave(filename = plot_path, plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)
  
}

RESULTS_PATH_RNA <- "~/ELNID001000000000486355/DiffAnalysis/RNA_oneOnOne/"
res <- read_csv("~/ELNID001000000000486355/RNAseq/Diff_analysis/diff.Medium.fatigue_vs_High.fatigue/edgeR_results_allgenes.csv")

volcano_plot(res)
volcano_plot(res, paste0(RESULTS_PATH_RNA, "volcano_plot.png"))
##rank differentially expressed genes
geneRank <- res %>%
  # calculate a combined score
  mutate(score = -log10(pvalue) * sign(log2FoldChange)) %>%
  # Rank genes by the combined score decreasingly
  arrange(desc(score))

rankings <-  geneRank$score
names(rankings) <- geneRank$gene
rankings <- sort(rankings, decreasing = TRUE)
rankings <- rankings[!is.na(names(rankings)) & names(rankings) != ""]
rankings <- rankings[!duplicated(names(rankings))]

##obtain pathways database
msigdbr_df <- msigdbr(species = "Homo sapiens")
msigdbr_gobp <- msigdbr_df[msigdbr_df$gs_subcollection == "GO:BP", ]
msigdbr_gobp_list <- split(x = msigdbr_gobp$gene_symbol, f = msigdbr_gobp$gs_name)
##perform pathway enrichment analysis
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
