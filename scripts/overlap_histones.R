library(ChIPseeker)
library(clusterProfiler)
library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # change as needed
library(org.Hs.eg.db)
library(circlize)
library("ATACseqQC")
library(GenomicAlignments)
library(tidyr)

atac_peaks_file <- "~/ELNID001000000000486355/analysis/ATAC/differentiallAnalysis_ATACseq_HighVSLow_noSexChr.csv"

data_path <- "~/ELNID001000000000535807/data/public_data/"
histones_metadata_file <- paste0(data_path, "metadata.tsv")
histones_metadata <- read.table(histones_metadata_file, header = TRUE, sep = "\t")

histones_metadata <- histones_metadata %>%
  dplyr::select(-c(Biosample.treatments, Biosample.genetic.modifications.categories, Biosample.genetic.modifications.gene.targets, 
                   Biosample.treatments.duration, Biosample.treatments.amount, Biosample.genetic.modifications.methods, 
                   Biosample.genetic.modifications.targets, Biosample.genetic.modifications.site.coordinates, Biosample.genetic.modifications.zygosity,
                   Library.made.from, Library.depleted.in, Library.extraction.method, Library.lysis.method, Library.crosslinking.method, 
                   Library.strand.specific, Experiment.date.released, Project, RBNS.protein.concentration, Library.fragmentation.method, Library.size.range, 
                   Biological.replicate.s., Technical.replicate.s., Read.length, Mapped.read.length, Run.type, Paired.end, Paired.with, Index.of, 
                   Derived.from, Genome.annotation, Platform, Controlled.by)) %>%
  dplyr::mutate(Experiment.target = gsub("-human", "", Experiment.target)) %>%
  dplyr::filter(File.format == "bigBed narrowPeak", File.analysis.status == "released", File.assembly == "GRCh38", Audit.ERROR == "") %>%
  dplyr::arrange(Experiment.target)# %>%
#  dplyr::filter(!grepl("insufficient", Audit.NOT_COMPLIANT))
row.names(histones_metadata) <- make.unique(histones_metadata$Experiment.target)

gr_list <- lapply(histones_metadata$File.accession, function(acc){
  bigbed_file <- paste0(data_path, acc, ".bigBed")
  gr <- import(bigbed_file)
})
gr_list <- GRangesList(gr_list)
names(gr_list) <- row.names(histones_metadata)

# #vennplot(gr_list[15:18], by = "gplots")
# get_setcounts <- function(sub_gr_list, labels){
#   # Compute overlap matrix for all regions across all sets
#   all_regions <- Reduce(c, sub_gr_list)
#   n_sets <- length(sub_gr_list)
#   n_regions <- length(all_regions)
#   
#   # Matrix: rows=regions, cols=sets (TRUE if region overlaps set i)
#   overlap_matrix <- vapply(sub_gr_list, function(gr) {
#     overlapsAny(all_regions, gr)
#   }, logical(n_regions))
#   
#   # Create binary pattern string for each region (e.g. "101" for sets 1+3)
#   overlap_pattern <- apply(overlap_matrix, 1, paste0, collapse = "")
#   
#   # Count regions per overlap pattern
#   pattern_counts <- table(overlap_pattern)
#   
#   # Convert patterns to eulerr format (e.g. "101" -> "Samp1&Samp3")
#   set_names <- character(length(pattern_counts))
#   names(set_names) <- names(pattern_counts)
#   
#   for (pattern in names(pattern_counts)) {
#     pattern_bin <- gsub("TRUE", "1", gsub("FALSE", "0", pattern))
#     active_sets <- which(strsplit(pattern_bin, "")[[1]] == "1")
#     if (length(active_sets) == 0) next
#     set_labels <- labels[active_sets]
#     set_names[pattern] <- paste(set_labels, collapse = "&")
#   }
#   
#   
#   print(pattern_counts)
#   # Return named vector for eulerr (exclude "000...0" pattern)
#   set_counts <- as.numeric(pattern_counts[set_names != ""])
#   names(set_counts) <- set_names[set_names != ""]
#   print(set_counts)
#   return(set_counts)
# }
# sub_gr_list <- gr_list[15:18]
# labels <- c("H3K9me3_1", "H3K9me3_2", "H3K9me3_3", "H3K9me3_4")
# 
# set_counts <- get_setcounts(sub_gr_list, labels)
# # Prepare eulerr input
# fit <- euler(set_counts)
# 
# # Plot
# plot(fit,
#      fills = list(alpha = 0.5, fill = c("#fbb4ae", "#b3cde3", "#decbe4", "#ccebc5", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2")),
#      labels = labels,
#      quantities = TRUE,
#      main = "H3K9me3 Peak Overlaps")

histone_groups <- list("H3K27ac" = c(1,2),
                    "H3K27me3" = c(3,4,5), 
                    "H3K36me3" = c(6,7),
                    "H3K4me1" = c(8,9),
                    "H3K4me3" = c(10,11,12), 
                    "H3K9ac" = c(13,14), 
                    "H3K9me3" = c(15,16,17,18))


merge_histone_group <- function(hist, histone_groups, gr_list, outdir) { 
  idx <- histone_groups[[hist]] 
  reps <- gr_list[idx] 
  n <- length(reps) 
  if (n == 2) { 
    # Direct pairwise overlap 
    ol <- countOverlaps(reps[[1]], reps[[2]]) > 0 
    message("Overlap ", hist, ": ", sum(ol)) 
    merged <- reps[[1]][ol] 
  } else {
    # N ≥ 3: build union then require ≥2-replicate support 
    all_peaks <- GenomicRanges::reduce(unlist(reps))
    overlap_mat <- sapply(reps, function(gr) countOverlaps(all_peaks, gr) > 0) 
    consensus <- rowSums(overlap_mat) >= 2 
    merged <- GenomicRanges::reduce(all_peaks[consensus]) 
    message(hist, ": ", sum(consensus), " consensus peaks (≥2 replicates)") 
  } 
  strand(merged) <- "+"
  # outfile <- file.path(outdir, paste0(hist, "_merged_peaks.bed")) 
  # export(merged, outfile) 
  merged 
}

outdir <- "~/ELNID001000000000535807/data/public_data/merged_histones"

merged_hist <- lapply(
  names(histone_groups),
  merge_histone_group,
  histone_groups = histone_groups,
  gr_list = gr_list,
  outdir = outdir
)

merged_hist <- GRangesList(merged_hist)
names(merged_hist) <- names(histone_groups)
# Merge the ranges from all GRanges objects
merged_gr <-  reduce(unlist(merged_hist))

# Initialize metadata columns with logical values
metadata <- data.frame(matrix(NA, nrow = length(merged_gr), ncol = length(merged_hist)))
colnames(metadata) <- names(merged_hist)

# Loop through each GRanges object in the list and check for overlaps
for (i in seq_along(merged_hist)) {
  overlaps <- findOverlaps(merged_gr, merged_hist[[i]])  # Find overlaps between merged GRanges and current GRanges object
  metadata[, i] <- as.integer(1:length(merged_gr) %in% queryHits(overlaps))  # 1 if overlap exists, 0 otherwise
}

# Add metadata to the merged GRanges object
mcols(merged_gr) <- metadata

#overlap with top 100 differentially accessible atac peaks
atac_peaks <- read.csv(atac_peaks_file)
atac_peaks <- atac_peaks %>% separate(X, into = c("chr", "start", "end"), sep = "_", convert = TRUE)
top_idx <- order(atac_peaks$padj)[1:100]
top_atac_peaks <- atac_peaks[top_idx,]
atac_peaks_gr <- GRanges(
  seqnames = top_atac_peaks$chr,
  ranges   = IRanges(start = top_atac_peaks$start, end = top_atac_peaks$end),
  mcols    = top_atac_peaks[ , setdiff(names(top_atac_peaks), c("chr", "start", "end")), drop = FALSE]
)

hits <- findOverlaps(merged_gr, atac_peaks_gr)

# subset
mg <- merged_gr[queryHits(hits)]
ap <- atac_peaks_gr[subjectHits(hits)]

# combine metadata
mcols(mg) <- cbind(mcols(mg), mcols(ap))

#annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

atac_peak_annot <- annotatePeak(mg, TxDb=txdb, annoDb="org.Hs.eg.db")
anno_df <- as.data.frame(atac_peak_annot)

write.csv(anno_df, "~/ELNID001000000000535807/analysis/histone_overlap/top100_hf_atac_peaks_overlap_histones_annotated.csv")



#vennplot(gr_list[15:18], by = "gplots")
get_setcounts <- function(sub_gr_list, labels){
  # Compute overlap matrix for all regions across all sets
  all_regions <- Reduce(c, sub_gr_list)
  n_sets <- length(sub_gr_list)
  n_regions <- length(all_regions)

  # Matrix: rows=regions, cols=sets (TRUE if region overlaps set i)
  overlap_matrix <- vapply(sub_gr_list, function(gr) {
    overlapsAny(all_regions, gr)
  }, logical(n_regions))

  # Create binary pattern string for each region (e.g. "101" for sets 1+3)
  overlap_pattern <- apply(overlap_matrix, 1, paste0, collapse = "")

  # Count regions per overlap pattern
  pattern_counts <- table(overlap_pattern)

  # Convert patterns to eulerr format (e.g. "101" -> "Samp1&Samp3")
  set_names <- character(length(pattern_counts))
  names(set_names) <- names(pattern_counts)

  for (pattern in names(pattern_counts)) {
    pattern_bin <- gsub("TRUE", "1", gsub("FALSE", "0", pattern))
    active_sets <- which(strsplit(pattern_bin, "")[[1]] == "1")
    if (length(active_sets) == 0) next
    set_labels <- labels[active_sets]
    set_names[pattern] <- paste(set_labels, collapse = "&")
  }


  print(pattern_counts)
  # Return named vector for eulerr (exclude "000...0" pattern)
  set_counts <- as.numeric(pattern_counts[set_names != ""])
  names(set_counts) <- set_names[set_names != ""]
  print(set_counts)
  return(set_counts)
}
sub_gr_list <- gr_list[15:18]
labels <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3")

set_counts <- get_setcounts(merged_hist, labels)
# Prepare eulerr input
fit <- euler(set_counts)

# Plot
plot(fit,
     fills = list(alpha = 0.5, fill = c("#fbb4ae", "#b3cde3", "#decbe4", "#ccebc5", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2")),
     labels = labels,
     quantities = TRUE,
     main = "Histones Peak Overlaps")

