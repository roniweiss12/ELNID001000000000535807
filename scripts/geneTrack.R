plot_gene_track <- function(
    gene_name,
    atac_hf1_bw, atac_hf2_bw, atac_hf3_bw,
    atac_lf1_bw, atac_lf2_bw, atac_lf3_bw,
    atac_hf1_bed, atac_hf2_bed, atac_hf3_bed,
    atac_lf1_bed, atac_lf2_bed, atac_lf3_bed,
    ser_hf_bigwig_file, ser_mf_bigwig_file,
    ser_hf_bed, ser_mf_bed, 
    rna_hf_bw, 
    rna_lf_bw,
    genome = "hg38",
    edb = EnsDb.Hsapiens.v86,
    window = 50000, 
    ct_y_top = 0.05, atac_y_top = 1, rna_y_top = 1
) {
  library(Gviz)
  library(ensembldb)
  library(AnnotationFilter)
  library(EnsDb.Hsapiens.v86)
  library(rtracklayer)
  library(biomaRt)
  library(GenomicRanges)
  library(scales)  # for alpha transparency
  
  # --- 1. Get gene info
  gene_info <- genes(edb, filter = GeneNameFilter(gene_name))
  if (length(gene_info) == 0) return(NULL)
  if (length(gene_info) > 1) {
    warning("Multiple genes found. Using the first one.")
    gene_info <- gene_info[1]
  }
  
  chr <- as.character(seqnames(gene_info))
  region_start <- pmax((start(gene_info) - window), 0)
  region_end <- end(gene_info) + window
  
  # --- 2. Base tracks
  gtrack <- GenomeAxisTrack()
  tss_track <- AnnotationTrack(
    start = start(gene_info),
    end = end(gene_info),
    chromosome = chr,
    genome = genome,
    name = paste0("chr", chr),
    stacking = "squish",
    background.title = "#B7B7B7",
    strand = strand(gene_info)
  )
  
  # --- 3. Gene annotation track
  options(ucscChromosomeNames = TRUE)
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  grtrack <- BiomartGeneRegionTrack(
    genome = "GRCh38",
    chromosome = chr,
    start = region_start,
    end = region_end,
    transcriptAnnotation = "symbol",
    name = "Gene",
    biomart = mart,
    showId = TRUE,
    background.title = "#B7B7B7",
    collapseTranscripts = "meta",
    stacking = "squish"
  )
  
  # --- 4. Helper to import and build DataTrack
  import_bw <- function(file, chr, region_start, region_end, color, name, opacity = 0.3, ylim_top = 1) {
    if (!file.exists(file)) stop(paste("File not found:", file))
    gr <- rtracklayer::import(file, which = GRanges(paste0("chr", chr), IRanges(region_start, region_end)))
    strand(gr) <- "*"
    DataTrack(
      range = gr,
      genome = genome,
      chromosome = chr,
      type = "histogram",
      from = region_start,
      to = region_end,
      col.histogram = alpha(color, opacity),
      fill.histogram = alpha(color, opacity),
      window = 200,
      aggregation = "mean",
      ylim = c(0, ylim_top),
      background.title = "#B7B7B7",
      name = name
    )
  }
  
  # --- 5. BED tracks helper
  import_bed <- function(bed_file, chr, region_start, region_end, color = "blue", opacity = 0.5) {
    gr <- rtracklayer::import(bed_file)
    gr <- gr[seqnames(gr) == paste0("chr",chr)]  # filter by gene chromosome
    AnnotationTrack(
      start = region_start,
      end = region_end,
      range = gr,
      genome = genome,
      chromosome = chr,
      stacking = "squish",
      col = alpha(color, opacity),
      fill = alpha(color, opacity)
    )
  }
  
  # --- 6. Create overlay groups for ATAC-seq
  hf_atac_color <- "#ca0020"
  mf_atac_color <- "#f4a582"
  hf_ser_color <- "#0571b0"
  mf_ser_color <- "#92c5de"
  hf_rna_color <- "#1a9850"
  mf_rna_color <- "#91cf60"
  
  hf_tracks <- list(
    import_bw(atac_hf1_bw, chr, region_start, region_end, hf_atac_color, "ATAC - HF", ylim_top = atac_y_top),
    import_bw(atac_hf2_bw, chr, region_start, region_end, hf_atac_color, "ATAC - HF", ylim_top = atac_y_top),
    import_bw(atac_hf3_bw, chr, region_start, region_end, hf_atac_color, "ATAC - HF", ylim_top = atac_y_top)
  )
  mf_tracks <- list(
    import_bw(atac_lf1_bw, chr, region_start, region_end, mf_atac_color, "ATAC - LF", ylim_top = atac_y_top),
    import_bw(atac_lf2_bw, chr, region_start, region_end, mf_atac_color, "ATAC - LF", ylim_top = atac_y_top),
    import_bw(atac_lf3_bw, chr, region_start, region_end, mf_atac_color, "ATAC - LF", ylim_top = atac_y_top)
  )
  
  overlay_hf <- OverlayTrack(trackList = hf_tracks)
  overlay_mf <- OverlayTrack(trackList = mf_tracks)
  
  hf_atac_bed_tracks <- list(
    import_bed(atac_hf1_bed, chr, region_start, region_end, hf_atac_color),
    import_bed(atac_hf2_bed, chr, region_start, region_end, hf_atac_color),
    import_bed(atac_hf3_bed, chr, region_start, region_end, hf_atac_color)
  )
  mf_atac_bed_tracks <- list(
    import_bed(atac_lf1_bed, chr, region_start, region_end, mf_atac_color),
    import_bed(atac_lf2_bed, chr, region_start, region_end, mf_atac_color),
    import_bed(atac_lf3_bed, chr, region_start, region_end, mf_atac_color)
  )
  overlay_hf_atac_bed <- OverlayTrack(trackList = hf_atac_bed_tracks)
  overlay_mf_atac_bed <- OverlayTrack(trackList = mf_atac_bed_tracks)
  
  # --- 7. CUT&TAG BigWigs
  ser_hf <- import_bw(ser_hf_bigwig_file, chr, region_start, region_end, hf_ser_color, "H3Q5ser - HF", ylim_top = ct_y_top)
  ser_mf <- import_bw(ser_mf_bigwig_file, chr, region_start, region_end, mf_ser_color, "H3Q5ser - LF", ylim_top = ct_y_top)
  
  # --- 8. CUT&TAG BED tracks
  bed_ser_hf <- import_bed(ser_hf_bed, chr, region_start, region_end, color = hf_ser_color)
  bed_ser_mf <- import_bed(ser_mf_bed, chr, region_start, region_end, color = mf_ser_color)
  
  # --- 9. RNA-seq
  # hf_rna_tracks <- list(
  #   import_bw(rna_hf1_bw, chr, region_start, region_end, hf_color, "RNA-seq - HF", ylim_top = rna_y_top),
  #   import_bw(rna_hf2_bw, chr, region_start, region_end, hf_color, "RNA-seq - HF", ylim_top = rna_y_top),
  #   import_bw(rna_hf3_bw, chr, region_start, region_end, hf_color, "RNA-seq - HF", ylim_top = rna_y_top)
  # )
  # mf_rna_tracks <- list(
  #   import_bw(rna_lf1_bw, chr, region_start, region_end, mf_color, "RNA-seq - LF", ylim_top = rna_y_top),
  #   import_bw(rna_lf2_bw, chr, region_start, region_end, mf_color, "RNA-seq - LF", ylim_top = rna_y_top),
  #   import_bw(rna_lf3_bw, chr, region_start, region_end, mf_color, "RNA-seq - LF", ylim_top = rna_y_top)
  # )
  # 
  # overlay_hf_rna <- OverlayTrack(trackList = hf_rna_tracks)
  # overlay_mf_rna <- OverlayTrack(trackList = mf_rna_tracks)
  hf_rna <- import_bw(rna_hf_bw, chr, region_start, region_end, hf_rna_color, "RNA - HF", ylim_top = rna_y_top)
  mf_rna <- import_bw(rna_lf_bw, chr, region_start, region_end, mf_rna_color, "RNA - LF", ylim_top = rna_y_top)
  
  # --- 10. Plot all tracks
  plotTracks(
    list(
      gtrack, grtrack,
      ser_hf, bed_ser_hf, 
      ser_mf, bed_ser_mf,
      overlay_hf, overlay_hf_atac_bed, 
      overlay_mf,overlay_mf_atac_bed,
      hf_rna, mf_rna
    ),
    from = region_start,
    to = region_end,
    sizes = c(1, 1, 2, 0.5, 2, 0.5, 2, 0.5, 2, 0.5, 2, 2),
    main = bquote(plain(.(gene_name)))
  )
}
atac_hf1_bw <- "~/ELNID001000000000486355/ATACseq/BigWig/03_0KC9_02SLTytgat_High-fatigue-1_ATAC_hs_i29-48_REP1.mLb.clN.bigWig" 
atac_hf2_bw <- "~/ELNID001000000000535807/Q52867_ATACseq/BigWig/08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1.mLb.clN.bigWig" 
atac_hf3_bw <- "~/ELNID001000000000535807/Q52867_ATACseq/BigWig/09_0KXR_02UCAUMC_High-fatigue-3_ATAC_i7-15_REP1.mLb.clN.bigWig"
atac_lf1_bw <- "~/ELNID001000000000535807/ATACseq/BigWig/01_0KC8_02SLTytgat_Medium-fatigue-1_ATAC_hs_i28-47_REP1.mLb.clN.bigWig" 
atac_lf2_bw <- "~/ELNID001000000000535807/Q52867_ATACseq/BigWig/05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1.mLb.clN.bigWig" 
atac_lf3_bw <- "~/ELNID001000000000535807/Q52867_ATACseq/BigWig/06_0KXP_02UCAUMC_Medium-fatigue-3_ATAC_i5-521_REP1.mLb.clN.bigWig"
atac_hf1_bed <- "~/ELNID001000000000535807/ATACseq/PeakCalling/macs2/03_0KC9_02SLTytgat_High-fatigue-1_ATAC_hs_i29-48_REP1.mLb.clN_peaks.narrowPeak" 
atac_hf2_bed <- "~/ELNID001000000000535807/Q52867_ATACseq/PeakCalling/macs2/08_0KXQ_02UCAUMC_High-fatigue-2_ATAC_i6-522_REP1.mLb.clN_peaks.narrowPeak" 
atac_hf3_bed <- "~/ELNID001000000000535807/Q52867_ATACseq/PeakCalling/macs2/09_0KXR_02UCAUMC_High-fatigue-3_ATAC_i7-15_REP1.mLb.clN_peaks.narrowPeak"
atac_lf1_bed <- "~/ELNID001000000000535807/ATACseq/PeakCalling/macs2/01_0KC8_02SLTytgat_Medium-fatigue-1_ATAC_hs_i28-47_REP1.mLb.clN_peaks.narrowPeak" 
atac_lf2_bed <- "~/ELNID001000000000535807/Q52867_ATACseq/PeakCalling/macs2/05_0KXO_02UCAUMC_Medium-fatigue-2_ATAC_hs_i4-4_REP1.mLb.clN_peaks.narrowPeak" 
atac_lf3_bed <- "~/ELNID001000000000535807/Q52867_ATACseq/PeakCalling/macs2/06_0KXP_02UCAUMC_Medium-fatigue-3_ATAC_i5-521_REP1.mLb.clN_peaks.narrowPeak"
ser_hf_bigwig_file <- "~/ELNID001000000000535807/CUT_TAG/BigWig/03_0KCB_02SMTytgat_High-fatigue-1_H3Q5ser_hs_i702-508_R1.bigWig" 
ser_mf_bigwig_file <- "~/ELNID001000000000535807/CUT_TAG/BigWig/01_0KCA_02SMTytgat_Medium-fatigue-1_H3Q5ser_hs_i701-508_R1.bigWig"
ser_hf_bed_file <- "~/ELNID001000000000535807/CUT_TAG/PeakCalling/macs2/03_0KCB_02SMTytgat_High-fatigue-1_H3Q5ser_hs_i702-508_R1.macs2.peaks.cut.bed" 
ser_mf_bed_file <- "~/ELNID001000000000535807/CUT_TAG/PeakCalling/macs2/01_0KCA_02SMTytgat_Medium-fatigue-1_H3Q5ser_hs_i701-508_R1.macs2.peaks.cut.bed" 
rna_hf_bw <- "~/ELNID001000000000535807/RNAseq/BIGWIG/02_P-02SNTytgat_High-fatigue-1_RNA-seq_hs_iA-D8.bw" 
rna_lf_bw <- "~/ELNID001000000000535807/RNAseq/BIGWIG/01_P-02SNTytgat_Medium-fatigue-1_RNA-seq_hs_iA-C8.bw" 

plot_gene_track(
  "S100A13",
  atac_hf1_bw, atac_hf2_bw, atac_hf3_bw,
  atac_lf1_bw, atac_lf2_bw, atac_lf3_bw,
  atac_hf1_bed, atac_hf2_bed, atac_hf3_bed,
  atac_lf1_bed, atac_lf2_bed, atac_lf3_bed,
  ser_hf_bigwig_file, ser_mf_bigwig_file,
  ser_hf_bed_file, ser_mf_bed_file,
  rna_hf_bw, rna_lf_bw, 
  window = 1000, 
  ct_y_top = 0.03,
  atac_y_top = 0.4,
  rna_y_top = 150
)
