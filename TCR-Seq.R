
#scTCR-seq and scRNA-Seq data analysis and visualization

#Load necessary library packages
library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(Matrix)
library(cowplot)
theme_set(theme_cowplot())
library(scCustomize)
library(ggplot2)
library(viridis)
library(optparse)
library(MAST)
library(future)
library(qs)
library(dplyr)
library(patchwork)
library(ggrepel)
library(DoubletFinder)
library(limma)
library(clustree)
library(glue)
library(knitr)
library(BiocParallel)
library(writexl)
library(jsonlite)
library(vctrs)
library(conflicted)
library(qlcMatrix)
library(getopt)
library(dittoSeq)
library(parallel)
library(rhdf5)
library(harmony)
library(DropletUtils)
library(presto)
library(scRepertoire)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(SingleCellExperiment)
library(ggthemes)
library(dplyr)
library(RColorBrewer)



# Read the CSV files into data frames
# Assuming you have multiple annotation files
tcr_list <- list(
  TCR1_Control = read.csv("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/TCR_seq/TCR1_Control/filtered_contig_annotations.csv"),
  TCR2_Naive_T_cells = read.csv("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/TCR_seq/TCR2_Naive_T_cells/filtered_contig_annotations.csv"),
  TCR3_Adaptive_T_cells = read.csv("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/TCR_seq/TCR3_Adaptive_T_cells/filtered_contig_annotations.csv")
)


combined_tcr <- combineTCR(tcr_list, 
                           samples = c("TCR1_Control", "TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells"),
                           ID = c("TCR1", "TCR2", "TCR3"),  # Can use patient or sample IDs
                           filterMulti = TRUE)


# To reload:
#load("combined_TCR.RData")

#Loading Data into scRepertoire

TCR1 = read.csv("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/TCR_seq/TCR1_Control/filtered_contig_annotations.csv")
TCR2 = read.csv("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/TCR_seq/TCR2_Naive_T_cells/filtered_contig_annotations.csv")
TCR3 = read.csv("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/TCR_seq/TCR3_Adaptive_T_cells/filtered_contig_annotations.csv")

contig_list <- list(TCR1, TCR2, TCR3)

head(contig_list)
head(contig_list[[1]])

#Combining Contigs into Clones

combined_TCR <- combineTCR(contig_list, 
                           samples = c("TCR1_Control", "TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)


# Inspect first object to see structure
head(combined_TCR[[1]])
head(combined_TCR[[2]])
head(combined_TCR[[3]])

head(combined_TCR[[1]]$barcode)
head(combined_TCR[[2]]$barcode)
head(combined_TCR[[3]]$barcode)


# ============================================================
# scTCR-seq scRepertoire Visualization Pipeline
# - Uses your actual cloneCall columns: CTstrict / CTgene / CTaa / CTnt
# - Enforces ONE palette across ALL figures (strict/gene/aa/nt)
# ============================================================

# ----------------------------
# 0) Output directory
# ----------------------------
OUTDIR <- "figures_pub"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Publication theme (45° x-axis)
# ----------------------------
theme_pub_45 <- function(base_size = 12, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2, hjust = 0),
      axis.title    = element_text(face = "bold"),
      axis.text     = element_text(color = "black"),
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title  = element_text(face = "bold"),
      legend.text   = element_text(color = "black"),
      panel.grid    = element_blank(),
      plot.margin   = margin(8, 10, 8, 10)
    )
}

# ----------------------------
# 2) Palette (choose one)
# ----------------------------
PAL_NAME <- "okabe"  # "okabe" | "tableau10" | "dark2"

pal_okabe <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

pal_tableau10 <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"
)

expand_brewer <- function(n, brewer_name = "Dark2") {
  base <- brewer.pal(min(8, max(3, n)), brewer_name)
  colorRampPalette(base)(n)
}

make_pal <- function(n, pal_name = PAL_NAME) {
  if (pal_name == "okabe") {
    if (n <= length(pal_okabe)) return(pal_okabe[seq_len(n)])
    return(colorRampPalette(pal_okabe)(n))
  }
  if (pal_name == "tableau10") {
    if (n <= length(pal_tableau10)) return(pal_tableau10[seq_len(n)])
    return(colorRampPalette(pal_tableau10)(n))
  }
  expand_brewer(n, "Dark2")
}

# ----------------------------
# 3) SMART palette applier (prevents your warning + mismatched colors)
#    It reads plot levels and applies matching named colors,
#    overwriting any existing scales from scRepertoire.
# ----------------------------
apply_palette_smart <- function(p, prefer_levels = NULL, pal_name = PAL_NAME) {
  
  b <- ggplot_build(p)
  
  lev_fill  <- tryCatch(b$plot$scales$get_scales("fill")$get_limits(),  error = function(e) NULL)
  lev_color <- tryCatch(b$plot$scales$get_scales("colour")$get_limits(), error = function(e) NULL)
  
  lev <- lev_fill
  if (is.null(lev) || length(lev) == 0) lev <- lev_color
  
  # fallback order
  if ((is.null(lev) || length(lev) == 0) && !is.null(prefer_levels)) {
    lev <- prefer_levels
  }
  
  if (is.null(lev) || length(lev) == 0) return(p)
  
  cols <- make_pal(length(lev), pal_name)
  names(cols) <- lev
  
  # overwrite any existing scRepertoire scales
  p +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_color_manual(values = cols, drop = FALSE)
}

# ----------------------------
# 4) Finish plot = palette + theme (theme LAST)
# ----------------------------
finish_plot <- function(p, prefer_levels = NULL,
                        base_size = 12, base_family = "Helvetica") {
  
  p <- apply_palette_smart(p, prefer_levels = prefer_levels, pal_name = PAL_NAME)
  p <- p + theme_pub_45(base_size = base_size, base_family = base_family)
  p
}

# ----------------------------
# 5) Save helper (vector PDF)
# ----------------------------
save_pdf <- function(p, filename, width, height) {
  grDevices::cairo_pdf(filename, width = width, height = height)
  print(p)
  grDevices::dev.off()
}

# ============================================================
# 6) Define consistent levels for coloring (SAMPLES + TYPE)
#    These MUST match your combined_TCR$sample and combined_TCR$Type values.
# ============================================================
samples3 <- c("TCR1_Control", "TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells")
samples2_12 <- c("TCR1_Control", "TCR2_Naive_T_cells")
samples2_13 <- c("TCR1_Control", "TCR3_Adaptive_T_cells")
samples2_23 <- c("TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells")

# Your table shows Type = TCR1 etc. (edit if other values exist)
type_levels <- c("TCR1", "TCR2", "TCR3")


# ============================================================
# 7) BASIC CLONAL VISUALIZATIONS (Updated cloneCall columns)
# ============================================================

# ---- clonalQuant: STRICT (IMPORTANT FIX: color by sample) ----
p <- clonalQuant(combined_TCR,
                 cloneCall = "CTstrict",
                 chain = "both",
                 group.by = "sample",   # <-- critical to match gene plot colors
                 scale = TRUE)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalQuant_CTstrict_bySample.pdf"), 4.2, 3.5)

# ---- clonalQuant: GENE ----
p <- clonalQuant(combined_TCR,
                 cloneCall = "CTgene",
                 group.by = "sample",
                 scale = TRUE)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalQuant_CTgene_bySample.pdf"), 4.2, 3.5)


# ---- clonalAbundance ----
p <- clonalAbundance(combined_TCR, cloneCall = "CTgene", scale = FALSE)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalAbundance_CTgene_counts.pdf"), 5, 3)

p <- clonalAbundance(combined_TCR, cloneCall = "CTgene", scale = TRUE)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalAbundance_CTgene_scaled.pdf"), 5, 3)

p <- clonalAbundance(combined_TCR, cloneCall = "CTaa", scale = FALSE)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalAbundance_CTaa_counts.pdf"), 5, 3)

p <- clonalAbundance(combined_TCR, cloneCall = "CTaa", scale = TRUE)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalAbundance_CTaa_scaled1.pdf"), 5, 3)


# ---- clonalLength ----
p <- clonalLength(combined_TCR, cloneCall="CTaa", chain="both")
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalLength_CTaa_both.pdf"), 6, 4)

p <- clonalLength(combined_TCR, cloneCall="CTaa", chain="TRA", scale=TRUE)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalLength_CTaa_density_TRA.pdf"), 6, 4)

p <- clonalLength(combined_TCR, cloneCall="CTaa", chain="TRB", scale=TRUE)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalLength_CTaa_density_TRB.pdf"), 6, 4)


# ---- clonalCompare (alluvial) ----
p <- clonalCompare(combined_TCR, top.clones = 10, samples = samples3,
                   cloneCall="CTaa", graph="alluvial")
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalCompare_CTaa_alluvial_3samples.pdf"), 8, 6)

p <- clonalCompare(combined_TCR, top.clones = 10, samples = samples3,
                   cloneCall="CTgene", graph="alluvial")
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalCompare_CTgene_alluvial_3samples.pdf"), 8, 6)

# Pairwise
p <- clonalCompare(combined_TCR, top.clones = 10, samples = samples2_12,
                   cloneCall="CTaa", graph="alluvial")
p <- finish_plot(p, prefer_levels = samples2_12)
save_pdf(p, file.path(OUTDIR, "clonalCompare_CTaa_alluvial_TCR1_vs_TCR2.pdf"), 7, 7)

p <- clonalCompare(combined_TCR, top.clones = 10, samples = samples2_13,
                   cloneCall="CTaa", graph="alluvial")
p <- finish_plot(p, prefer_levels = samples2_13)
save_pdf(p, file.path(OUTDIR, "clonalCompare_CTaa_alluvial_TCR1_vs_TCR3.pdf"), 7, 7)

p <- clonalCompare(combined_TCR, top.clones = 10, samples = samples2_23,
                   cloneCall="CTaa", graph="alluvial")
p <- finish_plot(p, prefer_levels = samples2_23)
save_pdf(p, file.path(OUTDIR, "clonalCompare_CTaa_alluvial_TCR2_vs_TCR3.pdf"), 7, 7)


# ---- clonalScatter ----
p <- clonalScatter(combined_TCR, cloneCall="CTgene",
                   x.axis="TCR1_Control", y.axis="TCR2_Naive_T_cells",
                   dot.size="total", graph="proportion")
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalScatter_CTgene_TCR1_vs_TCR2.pdf"), 6, 4)

p <- clonalScatter(combined_TCR, cloneCall="CTaa",
                   x.axis="TCR1_Control", y.axis="TCR2_Naive_T_cells",
                   dot.size="total", graph="proportion")
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalScatter_CTaa_TCR1_vs_TCR2.pdf"), 6, 4)

p <- clonalScatter(combined_TCR, cloneCall="CTgene",
                   x.axis="TCR2_Naive_T_cells", y.axis="TCR3_Adaptive_T_cells",
                   dot.size="total", graph="proportion")
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalScatter_CTgene_TCR2_vs_TCR3.pdf"), 6, 4)

p <- clonalScatter(combined_TCR, cloneCall="CTgene",
                   x.axis="TCR1_Control", y.axis="TCR3_Adaptive_T_cells",
                   dot.size="total", graph="proportion")
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalScatter_CTgene_TCR1_vs_TCR3.pdf"), 6, 4)


# ============================================================
# 8) CLONAL SPACE / PROPORTIONS
# ============================================================

p <- clonalHomeostasis(combined_TCR, cloneCall="CTgene")
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalHomeostasis_CTgene.pdf"), 4, 3)

p <- clonalHomeostasis(combined_TCR, cloneCall="CTaa")
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalHomeostasis_CTaa.pdf"), 4, 3)

# Group by Type (uses your Type values)
p <- clonalHomeostasis(combined_TCR, group.by="Type", cloneCall="CTgene")
p <- finish_plot(p, prefer_levels = type_levels)
save_pdf(p, file.path(OUTDIR, "clonalHomeostasis_CTgene_byType.pdf"), 4, 3)

p <- clonalProportion(combined_TCR, cloneCall="CTgene")
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalProportion_CTgene.pdf"), 4, 4)

p <- clonalProportion(combined_TCR, cloneCall="CTaa")
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalProportion_CTaa.pdf"), 4, 4)

p <- clonalProportion(combined_TCR, cloneCall="CTnt",
                      clonalSplit = c(1,5,10,100,1000,10000))
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "clonalProportion_CTnt.pdf"), 4, 4)


# ============================================================
# 9) REPERTOIRE SUMMARIES
# ============================================================

p <- percentAA(combined_TCR, chain="TRB", aa.length=20)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "percentAA_TRB.pdf"), 6, 5.5)

p <- percentAA(combined_TCR, chain="TRA", aa.length=20)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "percentAA_TRA.pdf"), 6, 5.5)

p <- positionalEntropy(combined_TCR, chain="TRB", aa.length=20)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "positionalEntropy_TRB.pdf"), 6, 5)

p <- positionalEntropy(combined_TCR, chain="TRA", aa.length=20)
p <- finish_plot(p)
save_pdf(p, file.path(OUTDIR, "positionalEntropy_TRA.pdf"), 6, 5)


# ============================================================
# 10) POSITIONAL PROPERTY (Atchley factors)
#     (these are usually 2-sample comparisons)
# ============================================================

pair_cols <- setNames(make_pal(2, PAL_NAME), samples2_12)

p <- positionalProperty(combined_TCR[c("TCR1_Control","TCR2_Naive_T_cells")],
                        chain="TRB", aa.length=20, method="atchleyFactors") +
  scale_color_manual(values = pair_cols) +
  theme_pub_45()

save_pdf(p, file.path(OUTDIR, "positionalProperty_TRB_TCR1_vs_TCR2.pdf"), 8, 6)

p <- positionalProperty(combined_TCR[c("TCR1_Control","TCR2_Naive_T_cells")],
                        chain="TRA", aa.length=20, method="atchleyFactors") +
  scale_color_manual(values = pair_cols) +
  theme_pub_45()

save_pdf(p, file.path(OUTDIR, "positionalProperty_TRA_TCR1_vs_TCR2.pdf"), 8, 6)


# ----------------------------
# Additional pairwise comparisons (TRB + TRA)
# ----------------------------

pair_list <- list(
  TCR1_vs_TCR3 = c("TCR1_Control", "TCR3_Adaptive_T_cells"),
  TCR2_vs_TCR3 = c("TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells")
)

for (nm in names(pair_list)) {
  
  pair <- pair_list[[nm]]
  pair_cols <- setNames(make_pal(2, PAL_NAME), pair)
  
  # ---- TRB ----
  p_trb <- positionalProperty(combined_TCR[pair],
                              chain = "TRB",
                              aa.length = 20,
                              method = "atchleyFactors") +
    scale_color_manual(values = pair_cols) +
    theme_pub_45()
  
  save_pdf(p_trb, file.path(OUTDIR, paste0("positionalProperty_TRB_", nm, ".pdf")), 8, 6)
  
  # ---- TRA ----
  p_tra <- positionalProperty(combined_TCR[pair],
                              chain = "TRA",
                              aa.length = 20,
                              method = "atchleyFactors") +
    scale_color_manual(values = pair_cols) +
    theme_pub_45()
  
  save_pdf(p_tra, file.path(OUTDIR, paste0("positionalProperty_TRA_", nm, ".pdf")), 8, 6)
}



# ============================================================
# 11) GENE USAGE
# ============================================================

p <- vizGenes(combined_TCR, x.axis="TRBV", plot="barplot", summary.fun="proportion") +
  theme_pub_45()
save_pdf(p, file.path(OUTDIR, "vizGenes_TRBV.pdf"), 3, 6)

p <- vizGenes(combined_TCR, x.axis="TRBD", plot="barplot", summary.fun="proportion") +
  theme_pub_45()
save_pdf(p, file.path(OUTDIR, "vizGenes_TRBD.pdf"), 3, 5.5)

p <- vizGenes(combined_TCR, x.axis="TRAV", plot="barplot", summary.fun="proportion") +
  theme_pub_45()
save_pdf(p, file.path(OUTDIR, "vizGenes_TRAV.pdf"), 3, 6)

p <- vizGenes(combined_TCR[samples3],
              x.axis="TRBV", y.axis="TRBJ",
              plot="heatmap", summary.fun="percent") +
  theme_pub_45()
save_pdf(p, file.path(OUTDIR, "vizGenes_heatmap_TRBV_TRBJ.pdf"), 7, 3)

p <- vizGenes(combined_TCR[samples3],
              x.axis="TRAV", y.axis="TRAJ",
              plot="heatmap", summary.fun="percent") +
  theme_pub_45()
save_pdf(p, file.path(OUTDIR, "vizGenes_heatmap_TRAV_TRAJ.pdf"), 9, 3)


# ============================================================
# 12) DIVERSITY / RAREFACTION / SIZE / OVERLAP (FIXED COLORS)
# ============================================================

# ---- clonalDiversity (FORCE sample mapping) ----
p <- clonalDiversity(
  combined_TCR,
  cloneCall = "CTgene",
  group.by  = "sample"
)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalDiversity_CTgene_bySample.pdf"), 4, 4)


# ---- clonalRarefaction (sample-aware) ----
p <- clonalRarefaction(
  combined_TCR,
  plot.type    = 1,
  hill.numbers = 0,
  n.boots      = 2
)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalRarefaction_q0_type1.pdf"), 6, 4)


# ---- clonalSizeDistribution (clustered, still enforce palette) ----
p <- clonalSizeDistribution(
  combined_TCR,
  cloneCall = "CTaa",
  method    = "ward.D2"
)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalSizeDistribution_CTaa.pdf"), 5, 3)


# ---- clonalOverlap: Morisita (sample-aware) ----
p <- clonalOverlap(
  combined_TCR,
  cloneCall = "CTstrict",
  method    = "morisita"
)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalOverlap_CTstrict_morisita.pdf"), 4.8, 2.8)


# ---- clonalOverlap: Raw overlap (sample-aware) ----
p <- clonalOverlap(
  combined_TCR,
  cloneCall = "CTstrict",
  method    = "raw"
)
p <- finish_plot(p, prefer_levels = samples3)
save_pdf(p, file.path(OUTDIR, "clonalOverlap_CTstrict_raw.pdf"), 4.8, 2.8)


message("colors now consistent across all figures")






################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

# ============================================================
# 0) Output directory
# ============================================================
OUTDIR <- "figures_pub"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) Publication theme (45° x-axis)
# ============================================================
theme_pub_45 <- function(base_size = 12, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2, hjust = 0),
      axis.title    = element_text(face = "bold"),
      axis.text     = element_text(color = "black"),
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title  = element_text(face = "bold"),
      legend.text   = element_text(color = "black"),
      panel.grid    = element_blank(),
      plot.margin   = margin(8, 10, 8, 10)
    )
}

# ============================================================
# 2) Palette (choose one)
# ============================================================
PAL_NAME <- "okabe"  # "okabe" | "tableau10" | "dark2"

pal_okabe <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

pal_tableau10 <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"
)

expand_brewer <- function(n, brewer_name = "Dark2") {
  base <- brewer.pal(min(8, max(3, n)), brewer_name)
  colorRampPalette(base)(n)
}

make_pal <- function(n, pal_name = PAL_NAME) {
  if (pal_name == "okabe") {
    if (n <= length(pal_okabe)) return(pal_okabe[seq_len(n)])
    return(colorRampPalette(pal_okabe)(n))
  }
  if (pal_name == "tableau10") {
    if (n <= length(pal_tableau10)) return(pal_tableau10[seq_len(n)])
    return(colorRampPalette(pal_tableau10)(n))
  }
  expand_brewer(n, "Dark2")
}

# ============================================================
# 3) SAFER SMART palette applier
#    - Applies only to DISCRETE scales (prevents continuous-scale error)
# ============================================================
.is_discrete_limits <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  if (is.numeric(x)) return(FALSE)
  if (is.logical(x)) return(FALSE)
  TRUE
}

apply_palette_smart <- function(p, prefer_levels = NULL, pal_name = PAL_NAME) {
  
  b <- ggplot_build(p)
  
  lev_fill  <- tryCatch(b$plot$scales$get_scales("fill")$get_limits(),  error = function(e) NULL)
  lev_color <- tryCatch(b$plot$scales$get_scales("colour")$get_limits(), error = function(e) NULL)
  
  lev <- lev_fill
  if (!.is_discrete_limits(lev)) lev <- lev_color
  if (!.is_discrete_limits(lev) && !is.null(prefer_levels)) lev <- prefer_levels
  
  if (!.is_discrete_limits(lev)) return(p)
  
  cols <- make_pal(length(lev), pal_name)
  names(cols) <- lev
  
  # overwrite any existing scales
  p +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_color_manual(values = cols, drop = FALSE)
}

# ============================================================
# 4) Finish plot (palette + theme)
# ============================================================
finish_plot <- function(p, prefer_levels = NULL,
                        base_size = 12, base_family = "Helvetica") {
  p <- apply_palette_smart(p, prefer_levels = prefer_levels, pal_name = PAL_NAME)
  p <- p + theme_pub_45(base_size = base_size, base_family = base_family)
  p
}

# ============================================================
# 5) Save helper (vector PDF)
# ============================================================
save_pdf <- function(p, filename, width, height) {
  grDevices::cairo_pdf(filename, width = width, height = height)
  print(p)
  grDevices::dev.off()
}


# ============================================================
# 6) Consistent levels + LABELS for sample + Type
# ============================================================
samples3 <- c("TCR1_Control", "TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells")
samples2_12 <- c("TCR1_Control", "TCR2_Naive_T_cells")
samples2_13 <- c("TCR1_Control", "TCR3_Adaptive_T_cells")
samples2_23 <- c("TCR2_Naive_T_cells", "TCR3_Adaptive_T_cells")

type_levels <- c("TCR1", "TCR2", "TCR3")

# Pretty labels (EDIT if you want different wording)
sample_labels <- c(
  TCR1_Control          = "Control",
  TCR2_Naive_T_cells    = "Naive",
  TCR3_Adaptive_T_cells = "Adaptive"
)

type_labels <- c(
  TCR1 = "TCR1",
  TCR2 = "TCR2",
  TCR3 = "TCR3"
)

# Consistent SAMPLE colors (named!)
sample_cols <- make_pal(length(samples3), PAL_NAME)
names(sample_cols) <- samples3

# Helper: relabel facets for sample / CellType etc.
facet_labels <- function(var, labels_named_vec) {
  labeller::as_labeller(labels_named_vec, default = label_value)
}


# ============================================================
# 7) Consistent CellType colors (use YOUR Seurat palette)
#    - This assumes you already have `seurat_obj$CellType` and ct_cols defined
# ============================================================
seurat_obj <- readRDS("processed_seurat_object_final.rds")

celltypes <- levels(factor(seurat_obj$CellType))
ct_cols <- Seurat::DiscretePalette(length(celltypes), palette = "polychrome")
names(ct_cols) <- celltypes

# enforce order everywhere (Seurat + SCE metadata)
seurat_obj$CellType <- factor(seurat_obj$CellType, levels = celltypes)
Idents(seurat_obj) <- "CellType"

colData(sce)$CellType <- factor(colData(sce)$CellType, levels = celltypes)

# also enforce sample order + labels in SCE (VERY important)
colData(sce)$sample <- factor(colData(sce)$sample, levels = samples3)

# ============================================================
# 8) CloneSize palette (keep your 4-level order)
# ============================================================
clone_levels <- c(
  "Hyperexpanded (0.1 < X <= 1)",
  "Large (0.01 < X <= 0.1)",
  "Medium (0.001 < X <= 0.01)",
  "Small (1e-04 < X <= 0.001)"
)
sce$cloneSize <- factor(sce$cloneSize, levels = clone_levels)

clone_cols <- c(
  "Hyperexpanded (0.1 < X <= 1)" = "#D73027",
  "Large (0.01 < X <= 0.1)"      = "#FC8D59",
  "Medium (0.001 < X <= 0.01)"   = "#91BFDB",
  "Small (1e-04 < X <= 0.001)"   = "#4575B4"
)


# ============================================================
# 9) UPDATED FIGURES (matched sample labels + CellType colors)
# ============================================================
library(scater)
# ----------------------------
# A) clonalOverlay by CellType (publication)
# ----------------------------
p <- clonalOverlay(
  sce,
  reduction = "UMAP",
  freq.cutpoint = 0.01,
  bins = 10,
  label = "CellType"
)

p <- finish_plot(p, prefer_levels = celltypes) +
  scale_color_manual(values = ct_cols, drop = FALSE) +
  scale_fill_manual(values  = ct_cols, drop = FALSE) +
  labs(x = "UMAP1", y = "UMAP2")

save_pdf(p, file.path(OUTDIR, "clonalOverlay_byCellType.pdf"), 8, 6)

# ----------------------------
# B) UMAP colored by cloneSize (single panel)
# ----------------------------
p <- plotUMAP(sce, colour_by = "cloneSize") +
  scale_colour_manual(values = clone_cols, drop = TRUE, na.value = "grey80")

p <- finish_plot(p, prefer_levels = clone_levels) +
  labs(x = "UMAP1", y = "UMAP2")

save_pdf(p, file.path(OUTDIR, "UMAP_cloneSize_custom.pdf"), 8, 6)

# ----------------------------
# C) UMAP cloneSize faceted by SAMPLE (with pretty facet labels)
# ----------------------------
library(ggplot2)

umap_coords <- reducedDim(sce, "UMAP")

umap_df <- data.frame(
  UMAP_1   = umap_coords[, 1],
  UMAP_2   = umap_coords[, 2],
  sample   = factor(colData(sce)$sample, levels = samples3),
  cloneSize = factor(colData(sce)$cloneSize, levels = clone_levels)
)

umap_df <- umap_df[!is.na(umap_df$cloneSize), , drop = FALSE]

p_umap_facet <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = cloneSize)
) +
  geom_point(size = 0.5, alpha = 0.8) +
  facet_wrap(
    ~ sample,
    labeller = ggplot2::labeller(sample = sample_labels)
  ) +
  scale_color_manual(values = clone_cols, drop = TRUE) +
  labs(
    title = NULL,
    x = "UMAP1",
    y = "UMAP2"
  )

p_umap_facet <- finish_plot(p_umap_facet, prefer_levels = clone_levels)

save_pdf(
  p_umap_facet,
  file.path(OUTDIR, "UMAP_cloneSize_bySample.pdf"),
  width = 12,
  height = 5
)


# ----------------------------
# D) TSNE cloneSize faceted by SAMPLE (pretty labels)
# ----------------------------

tsne_coords <- reducedDim(sce, "TSNE")

tsne_df <- data.frame(
  TSNE_1   = tsne_coords[, 1],
  TSNE_2   = tsne_coords[, 2],
  sample   = factor(colData(sce)$sample, levels = samples3),
  cloneSize = factor(colData(sce)$cloneSize, levels = clone_levels)
)

tsne_df <- tsne_df[!is.na(tsne_df$cloneSize), , drop = FALSE]

p_tsne_facet <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, color = cloneSize)) +
  geom_point(size = 0.5, alpha = 0.8) +
  facet_wrap(
    ~ sample,
    labeller = ggplot2::labeller(sample = sample_labels)
  ) +
  scale_color_manual(values = clone_cols, drop = TRUE) +
  labs(title = NULL, x = "tSNE1", y = "tSNE2")

p_tsne_facet <- finish_plot(p_tsne_facet, prefer_levels = clone_levels)

save_pdf(
  p_tsne_facet,
  file.path(OUTDIR, "TSNE_cloneSize_bySample.pdf"),
  width = 12,
  height = 5
)


# ----------------------------
# E) clonalOverlay faceted by CellType / sample (scRepertoire)
#    NOTE: clonalOverlay returns ggplot; we can post-modify theme + facet labels
# ----------------------------
p <- clonalOverlay(
  sce,
  reduction = "UMAP",
  cutpoint = 1,
  bins = 10,
  facet.by = "CellType"
) + guides(color = "none")

p <- finish_plot(p, prefer_levels = celltypes)

save_pdf(p, file.path(OUTDIR, "clonalOverlay_facet_CellType_UMAP.pdf"), 12, 9)

p <- clonalOverlay(
  sce,
  reduction = "UMAP",
  cutpoint  = 1,
  bins      = 10,
  facet.by  = "sample"
) + guides(color = "none")

# Pretty facet strip labels (NO facet_labels helper)
p <- p + facet_wrap(
  ~ sample,
  labeller = ggplot2::labeller(sample = sample_labels)
)

p <- finish_plot(p, prefer_levels = samples3)

save_pdf(p, file.path(OUTDIR, "clonalOverlay_facet_sample_UMAP.pdf"), 12, 6)


# ----------------------------
# F) clonalOccupy (by CellType) with CellType x-axis 45° and your clone colors
# ----------------------------
p <- clonalOccupy(sce, x.axis = "CellType", label = FALSE) +
  scale_fill_manual(values = clone_cols) +
  labs(title = NULL, x = NULL, y = NULL)

p <- finish_plot(p, prefer_levels = clone_levels) +
  theme(legend.position = "right")

save_pdf(p, file.path(OUTDIR, "clonalOccup_byCellType.pdf"), 5, 3)

# ----------------------------
# G) clonalOccupy (by sample) with pretty sample labels
# ----------------------------
p <- clonalOccupy(sce, x.axis = "sample", label = FALSE) +
  scale_fill_manual(values = clone_cols) +
  scale_x_discrete(labels = sample_labels) +
  labs(title = NULL, x = NULL, y = NULL)

p <- finish_plot(p, prefer_levels = clone_levels) +
  theme(legend.position = "right")

save_pdf(p, file.path(OUTDIR, "clonalOccup_bySample.pdf"), 4, 3.5)

# ----------------------------
# H) alluvialClones (gene, colored by CellType) with SAME CellType colors
# ----------------------------
p <- alluvialClones(
  sce,
  cloneCall = "gene",
  y.axes = c("sample", "CellType"),
  color = "CellType"
) +
  scale_fill_manual(values = ct_cols, drop = FALSE) +
  scale_x_discrete(labels = sample_labels) +
  labs(title = NULL)

p <- finish_plot(p, prefer_levels = celltypes)

save_pdf(p, file.path(OUTDIR, "alluvialClones_bySample_gene_CellTypeColors.pdf"), 9, 6)





suppressPackageStartupMessages({
  library(scRepertoire)
  library(circlize)
  library(grDevices)
})

# Make sure these exist from your earlier code:
# - OUTDIR
# - samples3
# - sample_labels (named vector)
# - ct_cols (named vector of CellType colors)
# - colData(sce)$CellType is a factor with levels(ct_cols)

# Safety: enforce order + factor levels
colData(sce)$sample   <- factor(colData(sce)$sample, levels = samples3)
colData(sce)$CellType <- factor(colData(sce)$CellType, levels = names(ct_cols))

make_chord_by_sample <- function(sce, sample_id,
                                 group.by = "CellType",
                                 sector_cols = ct_cols,
                                 outdir = OUTDIR,
                                 width = 6, height = 6,
                                 title_cex = 1.0) {
  
  sub <- sce[, colData(sce)$sample == sample_id]
  
  # Circlize input (sectors will be levels of `group.by`)
  circles <- getCirclize(sub, group.by = group.by, proportion = TRUE)
  
  # Determine actual sector names present in the circles object
  # (row/col names are the sector labels)
  sectors <- union(rownames(circles), colnames(circles))
  sectors <- sectors[!is.na(sectors) & sectors != ""]
  
  # Build grid.col using your fixed palette; fallback to grey if missing
  grid.col <- sector_cols[sectors]
  names(grid.col) <- sectors
  grid.col[is.na(grid.col)] <- "grey80"
  
  # Pretty title
  nice_title <- if (!is.null(sample_labels[[sample_id]])) sample_labels[[sample_id]] else sample_id
  
  # Save
  outfile <- file.path(outdir, paste0("chord_diagram_", sample_id, "_", group.by, ".pdf"))
  grDevices::cairo_pdf(outfile, width = width, height = height)
  
  circos.clear()
  chordDiagram(
    x = circles,
    self.link = 1,
    grid.col = grid.col,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow"
  )
  
  # Add a clean title
  title(main = paste0("Chord diagram (", group.by, "): ", nice_title),
        cex.main = title_cex, font.main = 2)
  
  dev.off()
}

# ---- Run all 3 samples (sizes like you had) ----
make_chord_by_sample(sce, "TCR1_Control", width = 4,   height = 4,   title_cex = 1.1)
make_chord_by_sample(sce, "TCR2_Naive_T_cells", width = 4, height = 4, title_cex = 1.1)
make_chord_by_sample(sce, "TCR3_Adaptive_T_cells", width = 5, height = 5, title_cex = 1.1)

#conflicts_prefer(GenomicRanges::union)



















