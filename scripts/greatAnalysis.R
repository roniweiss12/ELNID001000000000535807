library("rGREAT")
library(GenomicRanges)
library(dplyr)
library(ggplot2)

hf_bed_file <- "~/ELNID001000000000535807/analysis/ATAC/peaks_files/consensus_peaks/peaks_highFatigue.consensus.bed"
lf_bed_file <- "~/ELNID001000000000535807/analysis/ATAC/peaks_files/consensus_peaks/peaks_mediumFatigue.consensus.bed"

peaks_hf <- read.table(hf_bed_file, sep = "\t", header = FALSE)
peaks_lf <- read.table(lf_bed_file, sep = "\t", header = FALSE)

gr_hf <- makeGRangesFromDataFrame(peaks_hf, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE)
gr_lf <- makeGRangesFromDataFrame(peaks_lf, seqnames.field = "V1", start.field = "V2", end.field = "V3", ignore.strand = TRUE)

all_peaks <- reduce(c(gr_hf, gr_lf))  

gr_list <- list(gr_hf, gr_lf)

enrich_tables <- lapply(gr_list, function(gr) {
  res <- great(gr, background = all_peaks, "c5.go.bp", "TxDb.Hsapiens.UCSC.hg38.knownGene")
  tb = getEnrichmentTable(res)
  sig_tb <- tb %>% 
    dplyr::filter(p_adjust <= 0.05) %>%
    dplyr::arrange(fold_enrichment) %>%
    mutate(id = gsub("GOBP_", "", id)) %>%
    mutate(id = factor(id, levels = id))
  
})
# write.csv(enrich_tables[[1]], "~/ELNID001000000000486355/analysis/ATAC/pathwayEnrichment_highFatigueGREAT.csv", quote = FALSE, row.names = FALSE)
# write.csv(enrich_tables[[2]], "~/ELNID001000000000486355/analysis/ATAC/pathwayEnrichment_LOWFatigueGREAT.csv", quote = FALSE, row.names = FALSE)
     
p <- ggplot(enrich_tables[[1]] %>% dplyr::filter(p_adjust_hyper <= 0.05), aes(x=fold_enrichment, y=id, color = p_value, size = observed_gene_hits)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 14, hjust = 0.5))
ggsave("~/ELNID001000000000535807/analysis/ATAC/pathwayEnrichment_highFatigueGREATfull.png", plot = p, dpi = 400, width = 12, height = 8)

p <- ggplot(enrich_tables[[2]], aes(x=fold_enrichment, y=id, color = p_value, size = observed_gene_hits)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 14, hjust = 0.5))
ggsave("~/ELNID001000000000535807/analysis/ATAC/pathwayEnrichment_lowFatigueGREAT.png", plot = p, dpi = 400, width = 12, height = 8)