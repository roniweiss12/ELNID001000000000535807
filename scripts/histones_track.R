plot_gene_track_pg <- function(
    chr, region_start, region_end,
    atac_highFatigue_bw, atac_mediumFatigue_bw,
    bw_files,
    rna_hf_bw, rna_lf_bw,
    genome = "hg38",
    window = 50000,
    atac_y_top = 1,
    rna_y_top = 2.5
) {
  
  library(plotgardener)
  library(rtracklayer)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  
  # Expand region
  region_start <- max(region_start - window, 0)
  region_end <- region_end + window
  
  chrom <- paste0("chr", chr)
  
  ## --- PAGE SETUP
  pageCreate(width = 8, height = 20, default.units = "inches", xgrid = 0, ygrid = 0)
  
  y <- 0.5
  track_height <- 0.6
  gap <- 0.15
  
  ## --- Helper: BigWig plotting
  plot_bw_overlay <- function(bw_files, color, y, ymax) {
    for (bw in bw_files) {
      plotSignal(
        data = bw,
        chrom = chrom,
        chromstart = region_start,
        chromend = region_end,
        assembly = genome,
        x = 1, y = y,
        width = 6, height = track_height,
        linecolor = color,
        fill = color,
        alpha = 0.83,
        range = ymax
      )
    }
  }
  
  ## --- Genome axis
  plotGenomeLabel(
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = genome,
    x = 1, y = y,
    length = 6
  )
  y <- y + track_height 
  
  ## --- Gene track
  plotGenes(
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = genome,
    x = 1, y = y,
    width = 6,
    height = track_height
  )
  y <- y + track_height + gap
  
  ## --- ATAC HF
  plot_bw_overlay(atac_highFatigue_bw, "#ca0020", y, c(0, atac_y_top))
  plotText(label = "ATAC HF", fontcolor = "#ca0020", fontsize = 8,
           x = 1, y = y + gap, just = c("left"), default.units = "inches")
  y <- y + track_height
  
  ## --- ATAC LF
  plot_bw_overlay(atac_mediumFatigue_bw, "#f4a582", y, c(0, atac_y_top))
  plotText(label = "ATAC LF", fontcolor = "#f4a582", fontsize = 8,
           x = 1, y = y + gap, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  ## --- RNA HF
  plotSignal(
    data = rna_hf_bw,
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = genome,
    x = 1, y = y,
    width = 6, height = track_height,
    fill = "#ca0020", linecolor = "#ca0020",
    alpha = 0.83,
    range = c(0, rna_y_top)
  )
  plotText(label = "RNA HF", fontcolor = "#ca0020", fontsize = 8,
           x = 1, y = y + gap, just = c("left"), default.units = "inches")
  y <- y + track_height
  
  ## --- RNA LF
  plotSignal(
    data = rna_lf_bw,
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = genome,
    x = 1, y = y,
    width = 6, height = track_height,
    fill = "#f4a582", linecolor = "#f4a582",
    alpha = 0.83,
    range = c(0, rna_y_top)
  )
  plotText(label = "RNA LF", fontcolor = "#f4a582", fontsize = 8,
           x = 1, y = y + gap, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  ## --- Histone marks
  ##Active promoters
  plot_bw_overlay(bw_files[["H3K4me3"]][["high"]], "#006837", y, NULL)
  plotText(label = "H3K4me3 HF", fontcolor = "#006837", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K4me3"]][["medium"]], "#006837", y, NULL)
  plotText(label = "H3K4me3 LF", fontcolor = "#006837", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K9ac"]][["high"]], "#31a354", y, NULL)
  plotText(label = "H3K9ac HF", fontcolor = "#31a354", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K9ac"]][["medium"]], "#31a354", y, NULL)
  plotText(label = "H3K9ac LF", fontcolor = "#31a354", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  ##Active enhancers
  plot_bw_overlay(bw_files[["H3K27ac"]][["high"]], "#7fcdbb", y, NULL)
  plotText(label = "H3K27ac HF", fontcolor = "#7fcdbb", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K27ac"]][["medium"]], "#7fcdbb", y, NULL)
  plotText(label = "H3K27ac LF", fontcolor = "#7fcdbb", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K4me1"]][["high"]], "#41b6c4", y, NULL)
  plotText(label = "H3K4me1 HF", fontcolor = "#41b6c4", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K4me1"]][["medium"]], "#41b6c4", y, NULL)
  plotText(label = "H3K4me1 LF", fontcolor = "#41b6c4", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K4me2"]][["high"]], "#2c7fb8", y, NULL)
  plotText(label = "H3K4me2 HF", fontcolor = "#2c7fb8", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K4me2"]][["medium"]], "#2c7fb8", y, NULL)
  plotText(label = "H3K4me2 LF", fontcolor = "#2c7fb8", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K122ac"]][["high"]], "#253494", y, NULL)
  plotText(label = "H3K122ac HF", fontcolor = "#253494", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K122ac"]][["medium"]], "#253494", y, NULL)
  plotText(label = "H3K122ac LF", fontcolor = "#253494", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  ##Active gene bodies
  plot_bw_overlay(bw_files[["H3K36me3"]][["high"]], "#54278f", y, NULL)
  plotText(label = "H3K36me3 HF", fontcolor = "#54278f", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K36me3"]][["medium"]], "#54278f", y, NULL)
  plotText(label = "H3K36me3 LF", fontcolor = "#54278f", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  ##Repressed chromatin
  plot_bw_overlay(bw_files[["H3K27me3"]][["high"]], "#a50f15", y, NULL)
  plotText(label = "H3K27me3 HF", fontcolor = "#a50f15", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K27me3"]][["medium"]], "#a50f15", y, NULL)
  plotText(label = "H3K27me3 LF", fontcolor = "#a50f15", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K9me3"]][["high"]], "#de2d26", y, NULL)
  plotText(label = "H3K9me3 HF", fontcolor = "#de2d26", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H3K9me3"]][["medium"]], "#de2d26", y, NULL)
  plotText(label = "H3K9me3 LF", fontcolor = "#de2d26", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H4K20me1"]][["high"]], "#fb6a4a", y, NULL)
  plotText(label = "H4K20me1 HF", fontcolor = "#fb6a4a", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")
  y <- y + track_height + gap
  
  plot_bw_overlay(bw_files[["H4K20me1"]][["medium"]], "#fb6a4a", y, NULL)
  plotText(label = "H4K20me1 LF", fontcolor = "#fb6a4a", fontsize = 8,
           x = 1, y = y, just = c("left"), default.units = "inches")

  ## Remove guides
  pageGuideHide()
}
atac_highFatigue_bw <- c(atac_hf1_bw <- "~/ELNID001000000000535807/data/ATACseq/BigWig/03_0KC9_02SLTytgat_High-fatigue-1_ATAC_hs_i29-48_REP1.mLb.clN.bigWig",
                         atac_hf2_bw <- "~/ELNID001000000000535807/data/Q52867_ATACseq/BigWig/08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1.mLb.clN.bigWig",
                         atac_hf3_bw <- "~/ELNID001000000000535807/data/Q52867_ATACseq/BigWig/09_0KXR_02UCAUMC_High-fatigue-3_ATAC_i7-15_REP1.mLb.clN.bigWig")
atac_mediumFatigue_bw <- c(atac_lf1_bw <- "~/ELNID001000000000535807/data/ATACseq/BigWig/01_0KC8_02SLTytgat_Medium-fatigue-1_ATAC_hs_i28-47_REP1.mLb.clN.bigWig",
                        atac_lf2_bw <- "~/ELNID001000000000535807/data/Q52867_ATACseq/BigWig/05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1.mLb.clN.bigWig",
                        atac_lf3_bw <- "~/ELNID001000000000535807/data/Q52867_ATACseq/BigWig/06_0KXP_02UCAUMC_Medium-fatigue-3_ATAC_i5-521_REP1.mLb.clN.bigWig")
rna_hf_bw <- "~/ELNID001000000000535807/data/RNAseq/BIGWIG/chr_fixed/02_P-02SNTytgat_High-fatigue-1_RNA-seq_hs_iA-D8_filtered.bw"
rna_lf_bw <- "~/ELNID001000000000535807/data/RNAseq/BIGWIG/chr_fixed/01_P-02SNTytgat_Medium-fatigue-1_RNA-seq_hs_iA-C8_filtered.bw"

samples <- list(
  high = c(
    "03_0KC9_02SLTytgat_High-fatigue-1_ATAC_hs_i29-48_REP1",
    "08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1",
    "09_High-fatigue-3"
  ),
  medium = c(
    "01_0KC8_02SLTytgat_Medium-fatigue-1_ATAC_hs_i28-47_REP1",
    "05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1",
    "06_Medium-fatigue-3"
  )
)
marks <- c("H3K27ac","H3K4me1","H3K4me2","H3K4me3",
           "H3K122ac","H3K27me3","H3K36me3",
           "H3K9ac","H3K9me3","H4K20me1")
base_dir <- "~/ELNID001000000000535807/analysis/dHICA/output"

bw_files <- lapply(marks, function(mark) {
  list(
    high = file.path(base_dir, samples$high, "predicted_bw",
                     paste0("out-R-", mark, ".bw")),
    medium = file.path(base_dir, samples$medium, "predicted_bw",
                       paste0("out-R-", mark, ".bw"))
  )
})
names(bw_files) <- marks

png("~/ELNID001000000000535807/analysis/histone_overlap/geneTracks/ZSCAN18_HM.png", width = 8, height = 20, units = "in", res = 300)
plot_gene_track_pg(
  "19",58083838,  58118427,
  atac_highFatigue_bw, atac_mediumFatigue_bw,
  bw_files,
  rna_hf_bw, rna_lf_bw, 
  window = 10000, 
  atac_y_top = 3,
  rna_y_top = 500
)

dev.off()
