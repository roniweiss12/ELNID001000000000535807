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
set.seed(123)

filter_random_chr <- function(de_results){
  rn <- rownames(de_results)
  
  # Create a logical vector for rows with "chrUn" or "random"
  to_remove <- grepl("chrUn", rn) | grepl("random", rn)
  
  # Subset the data to exclude those rows
  filtered_res <- de_results[!to_remove, ]
  return(filtered_res)
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

count_table_file <- "~/ELNID001000000000535807/CUT_TAG/PeakCalling/featureCounts.txt"
diffAnalysisRes_file <- "~/ELNID001000000000535807/DiffAnalysis/CUT_TAG/differentiallAnalysis_CUT_TAG_HighVSLow.csv"
geneHancer_file <- "~/ELNID001000000000535807/DiffAnalysis/geneHancer_data/GeneHancer_v5.25.gff"
annotated_peaks_file <- "~/ELNID001000000000535807/DiffAnalysis/CUT_TAG/annotatedDifferentialPeaks.txt"
samples <- c("Low_fatigue", "High_fatigue")

##read in and process count table and metadata
##############################################
ct <- read.table(count_table_file, sep = "\t", header = TRUE)

# Convert each to GRanges
gr <- makeGRangesFromDataFrame(ct, keep.extra.columns = TRUE)

row.names(ct) <- paste0(ct$Chr, "_", ct$Start, "_", ct$End)
ct <- ct %>% dplyr::select(-Chr, -Start, -End, -Length, -Strand, -Geneid)
#rename medium fatigue as low
colnames(ct) <- samples
names(mcols(gr)) <- samples

metadata <- data.frame(group = samples, 
                       row.names = samples)
metadata$group <- factor(metadata$group, levels = samples)

# perform differentiall analysis
################################

# Construct DESeq2 dataset (even though DE cannot be done properly due to low sample size)
dds <- DESeqDataSetFromMatrix(countData = ct,
                              colData = metadata,
                              design = ~ group)

# Estimate size factors (library size normalization)
dds <- estimateSizeFactors(dds)

# Skip full DESeq fitting — not possible with 1 replicate per condition
# Instead, use normalized counts + pseudo-log transform
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Compute simple log2 fold change manually
group_levels <- levels(metadata$group)
lfc <- log2((norm_counts[, metadata$group == group_levels[2]] + 1) /
              (norm_counts[, metadata$group == group_levels[1]] + 1))

# Create a result-like data frame
res <- data.frame(
  gene = rownames(norm_counts),
  log2FoldChange = as.numeric(lfc),
  baseMean = rowMeans(norm_counts),
  pvalue = NA,
  padj = NA
)

# Sort by fold change magnitude (optional)
res <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]

# Subset the data to exclude those rows
filtered_res <- filter_random_chr(res)

#write.csv(filtered_res, diffAnalysisRes_file, row.names = TRUE, quote = FALSE)
filtered_res <- read.table(diffAnalysisRes_file, sep = ",", row.names = 1, header = TRUE)

count_table <- data.frame(assay(dds))

#save a processed version of count table
#write.csv(ct, "~/ELNID001000000000486355/DiffAnalysis/CUT_TAG/countTable.csv", row.names = TRUE, quote = FALSE)

#make clustering heatmap
########################

# Assume `mat` is a peaks x samples matrix
# Calculate variance across samples
peak_var <- rowVars(as.matrix(ct))

# Select top 1000–5000 most variable peaks
top_peaks <- order(peak_var, decreasing = TRUE)
mat_top <- ct[top_peaks, ]
# Create and store the heatmap object
colSide <- brewer.pal(3, "Set2")[metadata$group]
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap <- heatmap(as.matrix(mat_top), col = coul, ColSideColors = colSide, labRow = NA)


#obtain gene annotations from nearest TSS
########################################
#peak annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Annotate peaks
peak_annot <- annotatePeak(gr, TxDb=txdb, annoDb="org.Hs.eg.db")
anno_df <- as.data.frame(peak_annot)
#write.table(anno_df, annotated_peaks_file, row.names = FALSE, quote = FALSE, sep = "\t")
#merge with differential accessability data
filtered_res$rowname <- row.names(filtered_res)
anno_df <- anno_df %>%
  mutate(rowname = paste0(seqnames, "_", start, "_", end)) %>%
  right_join(filtered_res, by = "rowname")


#obtain gene annotations from geneHancer
########################################
peaks_gh_annotated <- add_geneHancer_annotation(gr, geneHancer_file)
#merge with genehancer annotation, if exists
anno_df_gh <- anno_df %>% rename(chr = "seqnames") %>% left_join(peaks_gh_annotated, by = c("chr", "start", "end"))
anno_df_gh$annotated_gene <- ifelse(!is.na(anno_df_gh$GHenhancer.connected_gene), anno_df_gh$GHenhancer.connected_gene, anno_df_gh$SYMBOL)
#write.table(anno_df_gh, gsub(".txt", "_GH.txt",annotated_peaks_file), row.names = FALSE, quote = FALSE, sep = "\t")
