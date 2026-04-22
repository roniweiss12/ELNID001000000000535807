library(rtracklayer)
library(mclust)
library(GenomicRanges)

bw <- import("~/ELNID001000000000486355/analysis/dHICA/output_lf/predicted_bw/out-R-H3K4me3.bw")

# # use scores
# x <- bw$score
# 
# # fit 2-component Gaussian mixture
# fit <- Mclust(x, G = 2)
# 
# # get means
# means <- fit$parameters$mean
# 
# # threshold = midpoint between clusters
# thr <- mean(means)
# 
# thr
quantile(bw$score, c(0.95, 0.98, 0.99))

peaks <- bw[bw$score > thr]
peaks <- reduce(peaks)

export(peaks, "~/ELNID001000000000535807/analysis/dHICA/output_lf/predicted_bed/H3K27me3_peaks.bed")


#remove extra columns from bed file - keep only chr start end

library(data.table)

# directory with bed files
bed_dir <- "/net/beegfs/groups/tytgat/ELNID001000000000535807/analysis/dHICA/output_hf/predicted_bed"

# get all .bed files
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

# process each file
for (f in bed_files) {
  dt <- fread(f, header = FALSE)
  
  # keep only first 3 columns
  dt <- dt[, 1:3, with = FALSE]
  
  # overwrite (or change filename if you want)
  fwrite(dt, paste0(f,"3"), sep = "\t", col.names = FALSE)
}
