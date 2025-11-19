library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(dplyr)

enriched_motifs_hf_file <- "/net/beegfs/groups/tytgat/ELNID001000000000486355/DiffAnalysis/ATAC/motif_analysis/diff_highFatigue_motifs/homerResults.csv"
enriched_motifs_lf_file <- "/net/beegfs/groups/tytgat/ELNID001000000000486355/DiffAnalysis/ATAC/motif_analysis/diff_lowFatigue_motifs/homerResults.csv"

#diffExpGenes_file <- "/net/beegfs/groups/tytgat/ELNID001000000000486355/DiffAnalysis/RNA_remission/differentiallAnalysis_RNAseq_HighVSLow.csv"

enriched_motifs_hf <- read.table(enriched_motifs_hf_file, sep = ",", header = TRUE)
enriched_motifs_lf <- read.table(enriched_motifs_lf_file, sep = "\t", header = TRUE)

enrich_go <- enrichGO(
  gene = enriched_motifs_hf$TF_name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",        # or "ENTREZID"
  ont = "BP",                # BP, CC, or MF
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)
enr_sum <- summary(enrich_go)

dotplot(enrich_go, title = "Gene Enrichment - Motifs Enriched in High Fatigue")

enrich_go <- enrichGO(
  gene = enriched_motifs_lf$TF_name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",        # or "ENTREZID"
  ont = "BP",                # BP, CC, or MF
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

dotplot(enrich_go, title = "Gene Enrichment - Motifs Enriched in Low Fatigue")
