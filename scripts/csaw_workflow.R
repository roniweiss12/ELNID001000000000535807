library(GenomicRanges)
library(csaw)
library("edgeR")

pe.bams <- c("~/ELNID001000000000486355/data/Q52867_ATACseq/BAM/05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1.mLb.clN.sorted.bam", "~/ELNID001000000000535807/data/Q52867_ATACseq/BAM/06_0KXP_02UCAUMC_Medium-fatigue-3_ATAC_i5-521_REP1.mLb.clN.sorted.bam",
	     "~/ELNID001000000000535807/data/Q52867_ATACseq/BAM/08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1.mLb.clN.sorted.bam", "~/ELNID001000000000535807/data/Q52867_ATACseq/BAM/09_0KXR_02UCAUMC_High-fatigue-3_ATAC_i7-15_REP1.mLb.clN.sorted.bam")

##############################
# read hg38 blacklist
blacklist <- read.table("~/ELNID001000000000535807/data/public_data/hg38-blacklist.v2.bed.gz", sep="\t")
colnames(blacklist) <- c("chrom", "start", "end", "reason")
blacklist <- GRanges(blacklist)
start(blacklist) <- start(blacklist) + 1

# define read parameters
standard.chr <- paste0("chr", c(1:22)) # only use standard autosomal chromosomes
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)

# few or no peaks should be removed; modify as desired
##############################
counts <- windowCounts(
  pe.bams,
  width = 150,        # ~1 nucleosome
  spacing = 50,
  param = param
)

neighbor <- resize(rowRanges(counts), 2000, "center")
wider <- regionCounts(pe.bams, neighbor, param=param)

filt <- filterWindowsLocal(counts, wider)
counts <- counts[filt$filter > log2(2), ]

binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)
counts <- normFactors(counts, se.out=TRUE, background=binned)

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(counts)
colnames(y$counts) <- c("lf1", "lf2", "hf1", "hf2")
rownames(y$samples) <- c("lf1", "lf2", "hf1", "hf2")
y$samples$group <- c("lf", "lf", "hf", "hf")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("hf", "lf") # CONFIRM THAT THESE COLUMNS CORRECTLY ALIGN!!
stopifnot(all(colnames(design) == c("hf", "lf")))
# design
# IMPORTANT: the user should manually confirm that the design matrix is correctly labeled according to sample metadata!

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(hf-lf, levels=design))
# head(results$table)
rowData(counts) <- cbind(rowData(counts), results$table) # combine GRanges rowdata with differential statistics
 
tab <- as.data.frame(results$table)
tab$PValue <- results$table$PValue
tab$FDR <- p.adjust(tab$PValue, method = "BH")
rowData(counts)$FDR <- tab$FDR

# choose FDR threshold for significant windows
sig.windows <- tab$FDR <= 0.05
# merge nearby significant windows; tol controls max gap (in bp)
merged <- mergeWindows(
  rowRanges(counts)[sig.windows],
  tol = 100,
  max.width = 5000
)

# summarize statistics across windows in each merged region
combined <- combineTests(
  id       = merged$ids,
  tab = tab[sig.windows, ]
)

# construct final peak-like GRanges
peak.regions <- merged$region
mcols(peak.regions) <- combined

# optional: attach average logCPM per region
logCPM <- cpm(y, log = TRUE, prior.count = 2)
region.means <- rowsum(
  logCPM[sig.windows, , drop = FALSE],
  group = merged$ids
) / as.numeric(table(merged$ids))
peak.regions$avgLogCPM <- rowMeans(region.means)

##############################
# write out csaw peaks
##############################

peak.df <- as.data.frame(peak.regions)
write.csv(peak.df, "~/ELNID001000000000535807/analysis/ATAC/hf_vs_lf_csaw_DA-windows_all.csv", quote=F, row.names=F)

sf <- counts$norm.factors
write.table(sf, "~/ELNID001000000000535807/analysis/ATAC/csaw_size_factors.txt",
            col.names=FALSE, row.names=FALSE)

###########################################

# Generate MA plot
library(ggplot2)

peak.df$sig <- "n.s."
peak.df$sig[peak.df$FDR < 0.05] <- "significant"

ggplot(data=peak.df,
       aes(x = avgLogCPM, y = rep.logFC, color = sig)) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = avgLogCPM, y = rep.logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)


keep <- results$table$logCPM > 0
final.peak.df <- peak.df[keep,]

#extract bed files for peaks (not differential)
# LF samples only (indices 1:2 based on your colnames: lf1, lf2, hf1, hf2)
lf_counts <- counts[, 1:2]

# Calculate average logCPM across LF replicates
lf_logcpm <- rowMeans(cpm(lf_counts, log=TRUE, prior.count=2))

# Define peaks as top 5-10% most accessible regions (adjust quantile as needed)
lf_peak_threshold <- quantile(lf_logcpm, 0.95)  # Top 5%
lf_peaks <- rowRanges(lf_counts)[lf_logcpm >= lf_peak_threshold]

# Convert to BED format and save
lf_bed <- data.frame(
  chrom = gsub("^chr", "", as.character(seqnames(lf_peaks))),
  start = as.integer(start(lf_peaks) - 1),  # BED is 0-based
  end = as.integer(end(lf_peaks))
)
write.table(lf_bed, "~/ELNID001000000000535807/analysis/ATAC/peaks_files/csaw/lf_peaks.bed",
            sep="\t", quote=F, col.names=F, row.names=F)

# HF samples only (indices 3:4)
hf_counts <- counts[, 3:4]
hf_logcpm <- rowMeans(cpm(hf_counts, log=TRUE, prior.count=2))
hf_peak_threshold <- quantile(hf_logcpm, 0.95)

hf_peaks <- rowRanges(hf_counts)[hf_logcpm >= hf_peak_threshold]
hf_bed <- data.frame(
  chrom = gsub("^chr", "", as.character(seqnames(hf_peaks))),
  start = as.integer(start(hf_peaks) - 1),
  end = as.integer(end(hf_peaks))
)
write.table(hf_bed, "~/ELNID001000000000535807/analysis/ATAC/peaks_files/csaw/hf_peaks.bed",
            sep="\t", quote=F, col.names=F, row.names=F)

