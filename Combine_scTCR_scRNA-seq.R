#scRNA-seq data analysis

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
theme_set(theme_cowplot())
library(scCustomize)
library(ggplot2)
library(viridis)
library(patchwork)
library(DoubletFinder)
library(harmony)


#Combine the seurat_obj and combine_TCR
load("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/combined_TCR.RData")

head(combined_TCR[[1]])  # Look at first sample if it’s a list
head(combined_TCR[[2]])
head(combined_TCR[[3]])

seurat_obj <- readRDS("processed_seurat_object_final.rds")
head(seurat_obj)


#Combining Clones and Single-Cell Object
library(scater)


# Check the first 10 variable features before removal
VariableFeatures(seurat_obj)[1:10]
# Remove TCR VDJ genes
seurat_obj <- quietTCRgenes(seurat_obj)

# Check the first 10 variable features after removal
VariableFeatures(seurat_obj)[1:10]

head(seurat_obj)

#Making a Single-Cell Experiment object
#sce <- Seurat::as.SingleCellExperiment(seurat_obj)
sce <- as.SingleCellExperiment(seurat_obj)

head(sce)
head(colData(sce))


colData(sce)$ident <- NULL

# Check object class
class(sce)
class(combined_TCR)


all_tcr <- do.call(rbind, combined_TCR)
head(all_tcr)
class(all_tcr)

head(colnames(seurat_obj))
head(all_tcr$barcode)


sce <- combineExpression(
  all_tcr,
  sce,
  cloneCall = "gene",         # or "aa" or "nt" depending on analysis goals
  group.by = "sample",        # consistent with your metadata
  proportion = TRUE
)

head(sce)
head(colData(sce))
rownames(all_tcr)
colnames(sce)


all(rownames(all_tcr) %in% colnames(sce))
rownames(all_tcr) <- all_tcr$barcode
head(all_tcr)


table(colData(sce)$cloneSize, useNA = "ifany")


colnames(colData(sce))
table(colData(sce)$CellType)

reducedDimNames(sce)

library(SingleCellExperiment)

# Create a matrix of UMAP coords
umap_coords <- as.matrix(colData(sce)[, c("UMAP_1", "UMAP_2")])
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Store it properly as reducedDim
reducedDim(sce, "UMAP") <- umap_coords
reducedDimNames(sce)

colData(sce)$ident <- colData(sce)$CellType

colnames(colData(sce))
# should now include "ident"


#by cell type
pdf("figures/clonalOverlay_byCellType.pdf", width = 8, height = 5)
clonalOverlay(
  sce,
  reduction = "UMAP",
  freq.cutpoint = 0.01,
  bins = 10,
  label = "CellType"
) + 
  theme_classic()  # or try theme_void() for no background at all
dev.off()

table(colData(sce)$cloneSize, useNA = "always")

# How many distinct values?
table(colData(sce)$cloneSize)

# If you want to be explicit:
unique(colData(sce)$cloneSize)

# Make sure it's a factor (it already is, but safe)
sce$cloneSize <- factor(sce$cloneSize)

# Define the levels you actually see in the data
clone_levels <- c(
  "Hyperexpanded (0.1 < X <= 1)",
  "Large (0.01 < X <= 0.1)",
  "Medium (0.001 < X <= 0.01)",
  "Small (1e-04 < X <= 0.001)"
)

# Optional: enforce this level order
sce$cloneSize <- factor(sce$cloneSize, levels = clone_levels)

# 4 colors from inferno
color_vec <- hcl.colors(
  n = length(clone_levels),
  palette = "inferno",
  fixup   = TRUE
)

# Name the colors so ggplot matches by level name
names(color_vec) <- clone_levels
color_vec
p <- plotUMAP(
  sce,
  colour_by = "cloneSize"
)

p <- p +
  scale_colour_manual(
    values   = color_vec,
    drop     = TRUE,        # drop unused levels like Rare / None
    na.value = "grey80"     # color for NA cloneSize
  ) +
  theme_classic()

p

# Save to PDF
pdf("figures/UMAP_cloneSize_custom.pdf", width = 8, height = 6)
print(p)
dev.off()



# Define color vector (adjust based on number of cloneSize levels present)
# Extract UMAP coordinates
umap_coords <- reducedDim(sce, "UMAP")
umap_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  sample = colData(sce)$sample,
  cloneSize = colData(sce)$cloneSize
)

# Drop NAs if necessary
umap_df <- umap_df[!is.na(umap_df$cloneSize), ]

# Use only levels present
levels_now <- unique(umap_df$cloneSize)
#color_vec <- setNames(rev(hcl.colors(n = length(levels_now), palette = "inferno")), levels_now)

# Explicit color mapping
color_vec <- c(
  "Hyperexpanded (0.1 < X <= 1)"   = "#D73027",  # Deep Red
  "Large (0.01 < X <= 0.1)" = "#FC8D59",  # Orange
  "Medium (0.001 < X <= 0.01)"   = "#91BFDB",  # Teal
  "Small (1e-04 < X <= 0.001)" = "#4575B4"  # Navy Blue
)

p_umap_facet <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cloneSize)) +
  geom_point(size = 0.5, alpha = 0.8) +
  facet_wrap(~ sample) +
  scale_color_manual(values = color_vec) +
  theme_classic() +
  labs(title = "UMAP by Clone Size", x = "UMAP 1", y = "UMAP 2")

pdf("figures/UMAP_cloneSize_bySample.pdf", width = 12, height = 5)
print(p_umap_facet)
dev.off()

# TSNE faceted by sample
# Extract TSNE coordinates
tsne_coords <- reducedDim(sce, "TSNE")
tsne_df <- data.frame(
  TSNE_1 = tsne_coords[, 1],
  TSNE_2 = tsne_coords[, 2],
  sample = colData(sce)$sample,
  cloneSize = colData(sce)$cloneSize
)

tsne_df <- tsne_df[!is.na(tsne_df$cloneSize), ]

# Use same color vector
p_tsne_facet <- ggplot(tsne_df, aes(x = TSNE_1, y = TSNE_2, color = cloneSize)) +
  geom_point(size = 0.5, alpha = 0.8) +
  facet_wrap(~ sample) +
  scale_color_manual(values = color_vec) +
  theme_classic() +
  labs(title = "TSNE by Clone Size", x = "TSNE 1", y = "TSNE 2")

# Save to PDF
pdf("figures/TSNE_cloneSize_bySample.pdf", width = 12, height = 5)
print(p_tsne_facet)
dev.off()


levels(factor(colData(sce)$cloneSize))


#Visualizations for Single-Cell Objects
#clonalOverlay
colnames(seurat_obj@meta.data)

pdf("figures/clonalOverlay_byCellType1.pdf", width = 12, height = 9)
clonalOverlay(sce, 
              reduction = "UMAP", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "CellType") + 
  guides(color = "none")
dev.off()

pdf("figures/clonalOverlay_byCellType2.pdf", width = 12, height = 9)
clonalOverlay(sce, 
              reduction = "TSNE", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "CellType") + 
  guides(color = "none")
dev.off()


pdf("figures/clonalOverlay_bysample.pdf", width = 12, height = 6)
clonalOverlay(sce, 
              reduction = "UMAP", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "sample") + 
  guides(color = "none")
dev.off()

pdf("figures/clonalOverlay_bysample1.pdf", width = 12, height = 6)
clonalOverlay(sce, 
              reduction = "TSNE", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "sample") + 
  guides(color = "none")
dev.off()


#clonalNetwork
#Filtering Options for clonalNetwork()
#ggraph needs to be loaded due to issues with ggplot
library(ggraph)

#No Identity filter
pdf("figures/clonalNetwork_byCellType.pdf", width = 10, height = 8)
clonalNetwork(sce, 
              reduction = "UMAP", 
              group.by = "CellType",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")
dev.off()


#Examining Cluster 3 only
pdf("figures/clonalNetwork_byCellType1.pdf", width = 5, height = 4)
clonalNetwork(sce, 
              reduction = "UMAP", 
              group.by = "CellType",
              filter.identity = 3,
              cloneCall = "aa")
dev.off()


pdf("figures/clonalNetwork_byCellType_tsne.pdf", width = 6, height = 4)
clonalNetwork(sce, 
              reduction = "TSNE", 
              group.by = "CellType",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")
dev.off()

pdf("figures/clonalNetwork_byCellType_tsne1.pdf", width = 5, height = 4)
clonalNetwork(sce, 
              reduction = "TSNE", 
              group.by = "CellType",
              filter.identity = 3,
              cloneCall = "aa")
dev.off()



shared.clones <- clonalNetwork(sce, 
                               reduction = "UMAP", 
                               group.by = "CellType",
                               cloneCall = "aa", 
                               exportClones = TRUE)
head(shared.clones)
write.csv(shared.clones, file = "shared_clones.csv", row.names = FALSE)

saveRDS(sce, file = "sce_object.rds")

sce <- readRDS("sce_object.rds")
sce

#highlightClones
sce <- highlightClones(sce, cloneCall= "aa", sequence = c("NA_CASSLDRVEQYF", "CATAGSGGKLTL_NA", "NA_CAWSLVQSGEQYF", "CAVSKSTNTGKLTF_NA", "NA_CASSKEGGRDEQYF"))

library(scater)
pdf("figures/UMAP_highlight_highlightClones.pdf", width = 6, height = 6)
plotUMAP(sce, colour_by = "highlight") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_classic() +
  theme(plot.title = element_blank(), legend.position = "bottom")
dev.off()



# Extract UMAP coordinates and metadata
umap <- reducedDim(sce, "UMAP")
df <- data.frame(
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  sample = colData(sce)$sample,
  highlight = colData(sce)$highlight
)

pdf("figures/UMAP_highlight_bySample.pdf", width = 10, height = 5)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = highlight)) +
  geom_point(size = 0.6, alpha = 0.9) +
  facet_wrap(~ sample) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_classic() +
  theme(plot.title = element_blank(), legend.position = "bottom")
dev.off()



pdf("figures/TSNE_highlight_highlightClones.pdf", width = 5, height = 5)
plotTSNE(sce, colour_by = "highlight") +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_classic() +
  theme(plot.title = element_blank(), legend.position = "bottom")
dev.off()


# Extract tSNE coordinates and metadata
tsne <- reducedDim(sce, "TSNE")
df_tsne <- data.frame(
  TSNE_1 = tsne[, 1],
  TSNE_2 = tsne[, 2],
  sample = colData(sce)$sample,
  highlight = colData(sce)$highlight
)

# Save PDF: Faceted tSNE by sample
pdf("figures/TSNE_highlight_bySample.pdf", width = 10, height = 5)
ggplot(df_tsne, aes(x = TSNE_1, y = TSNE_2, color = highlight)) +
  geom_point(size = 0.6, alpha = 0.9) +
  facet_wrap(~ sample) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_classic() +
  theme(plot.title = element_blank(), legend.position = "bottom")
dev.off()




#clonalOccupy
unique(colData(sce)$cloneSize)

#clonalOccupy(sce, x.axis = "CellType_cluster")


# Run clonalOccupy() and add custom colors
pdf("figures/clonalOccup_byCellType.pdf", width = 8, height = 6)
clonalOccupy(sce, x.axis = "CellType") +
  scale_fill_manual(
    values = c(
      "Hyperexpanded (0.1 < X <= 1)"   = "#D73027",  # deep red
      "Large (0.01 < X <= 0.1)"        = "#FC8D59",  # orange
      "Medium (0.001 < X <= 0.01)"     = "#91BFDB",  # light teal
      "Small (1e-04 < X <= 0.001)"     = "#4575B4"   # navy blue
    )
  ) +
  theme_classic() +
  theme(legend.position = "right")
dev.off()

#above figure without numbers on the bars and x-axis 45 degree text
pdf("figures/clonalOccup_byCellType_new.pdf", width = 5, height = 3)
clonalOccupy(sce, x.axis = "CellType", label = FALSE) +
  scale_fill_manual(values = c(
    "Hyperexpanded (0.1 < X <= 1)"   = "#D73027",
    "Large (0.01 < X <= 0.1)"        = "#FC8D59",
    "Medium (0.001 < X <= 0.01)"     = "#91BFDB",
    "Small (1e-04 < X <= 0.001)"     = "#4575B4"
  )) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.x = element_blank()
  )
dev.off()


pdf("figures/clonalOccup_bySample.pdf", width = 4, height = 3.5)
clonalOccupy(sce, x.axis = "sample", label = FALSE) +
  scale_fill_manual(
    values = c(
      "Hyperexpanded (0.1 < X <= 1)"   = "#D73027",  # deep red
      "Large (0.01 < X <= 0.1)"        = "#FC8D59",  # orange
      "Medium (0.001 < X <= 0.01)"     = "#91BFDB",  # light teal
      "Small (1e-04 < X <= 0.001)"     = "#4575B4"   # navy blue
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )
dev.off()



pdf("figures/clonalOccup_bySample_percentage.pdf", width = 5.5, height = 3)
clonalOccupy(sce, 
             x.axis = "CellType", 
             proportion = TRUE, 
             label = FALSE)+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )
dev.off()


pdf("figures/clonalOccup_bySample_percentage1.pdf", width = 6, height = 3)
clonalOccupy(
  sce,
  x.axis = "CellType",
  proportion = TRUE,
  label = FALSE
) +
  scale_fill_manual(
    values = c(
      "Hyperexpanded (0.1 < X <= 1)"   = "#D73027",  # deep red
      "Large (0.01 < X <= 0.1)"        = "#FC8D59",  # orange
      "Medium (0.001 < X <= 0.01)"     = "#91BFDB",  # light teal
      "Small (1e-04 < X <= 0.001)"     = "#4575B4"   # navy blue
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )
dev.off()


#alluvialClones
# Simple colorblind-friendly palette with 7 colors
colorblind_vector <- hcl.colors(n = 7, palette = "inferno", fixup = TRUE)

pdf("figures/alluvialClones_bySample.pdf", width = 8, height = 6)
alluvialClones(sce, 
               cloneCall = "aa", 
               y.axes = c("sample", "CellType"), 
               color = c("NA_CASSLDRVEQYF", "CATAGSGGKLTL_NA", "NA_CAWSLVQSGEQYF", "CAVSKSTNTGKLTF_NA", "NA_CASSKEGGRDEQYF")) + 
  scale_fill_manual(values = c("grey", colorblind_vector[3]))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )
dev.off()


pdf("figures/alluvialClones_bySample_gene.pdf", width = 8, height = 6)
alluvialClones(sce, 
               cloneCall = "gene", 
               y.axes = c("sample", "CellType", "CellType_cluster"), 
               color = "CellType") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  ) 
dev.off()

#Updated above two figures for publication
# Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# ---- A) Dark palette for MANY CellTypes ----
# Start from Dark2 (very dark), then extend with darker HCL colors if needed
make_dark_pal <- function(n) {
  base <- brewer.pal(min(8, n), "Dark2")
  if (n <= length(base)) return(base[seq_len(n)])
  extra <- hcl.colors(n - length(base), palette = "Dark 3")  # darker than "Set3"
  c(base, extra)
}

# ---- B) Dark 2-level palette (Selected vs Other) ----
# Publication-safe light grey + dark red
dark_two <- c(
  "Other"    = "lightgrey",  # light grey
  "Selected" = "darkgreen"   # dark red (deep, non-pink)
)


#alluvialClones by sample (cloneCall = "aa") with darker “Selected” highlight
# Your selected clones (keep your exact vector)
selected_clones <- c(
  "NA_CASSLDRVEQYF",
  "CATAGSGGKLTL_NA",
  "NA_CAWSLVQSGEQYF",
  "CAVSKSTNTGKLTF_NA",
  "NA_CASSKEGGRDEQYF"
)

grDevices::cairo_pdf("figures/alluvialClones_bySample_dark.pdf", width = 9, height = 6)

p1 <- alluvialClones(
  sce,
  cloneCall = "aa",
  y.axes = c("sample", "CellType", "CellType_cluster"),
  color = selected_clones
) +
  # alluvialClones typically maps fill to "Selected/Other" (2 groups)
  scale_fill_manual(values = dark_two) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_blank(),
    panel.border = element_rect(color = "light grey", fill = NA, linewidth = 0.6)
  )

print(p1)
dev.off()

#alluvialClones by sample (cloneCall = "gene", color = "CellType") with dark CellType palette
# Determine how many CellTypes exist in your object (adjust accessor if needed)
n_ct <- length(unique(colData(sce)$CellType))  # or: sce$CellType if stored that way
celltype_cols <- make_dark_pal(n_ct)
names(celltype_cols) <- sort(unique(colData(sce)$CellType))

grDevices::cairo_pdf("figures/alluvialClones_bySample_gene_dark.pdf", width = 9, height = 6)

p2 <- alluvialClones(
  sce,
  cloneCall = "gene",
  y.axes = c("sample", "CellType", "CellType_cluster"),
  color = "CellType"
) +
  scale_fill_manual(values = celltype_cols) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p2)
dev.off()

#Build a stable CellType palette (same as your barplot)
library(ggplot2)
library(scales)

# Make sure CellType is a factor with your desired order (from your barplot figure)
celltype_levels <- c(
  "cnkCD8","cT","efCD4","efCD8","exCD8","iNKT","mCD8",
  "nCD4","nCD4/CD8","nCD8","nkCD8","pT","TregCD4"
)

# Pull CellTypes present in your object and keep them in that order
cts_present <- intersect(celltype_levels, sort(unique(colData(sce)$CellType)))
colData(sce)$CellType <- factor(colData(sce)$CellType, levels = celltype_levels)

# Define a consistent palette (edit these hex codes once to match your barplot exactly)
celltype_cols <- c(
  cnkCD8   = "#4D4D4D",
  cT       = "#1B9E77",
  efCD4    = "#7570B3",
  efCD8    = "#E7298A",
  exCD8    = "#D95F02",
  iNKT     = "#66A61E",
  mCD8     = "#A6761D",
  nCD4     = "#1F78B4",
  `nCD4/CD8` = "#6A3D9A",
  nCD8     = "#33A02C",
  nkCD8    = "#B15928",
  pT       = "#E6AB02",
  TregCD4  = "#FB9A99"
)

# Keep only those present
celltype_cols <- celltype_cols[names(celltype_cols) %in% cts_present]

#Apply the same colors to alluvialClones (cloneCall = "gene", color = "CellType")
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

grDevices::cairo_pdf("figures/alluvialClones_bySample_gene_CellTypeColors.pdf",
                     width = 9, height = 6)

p2 <- alluvialClones(
  sce,
  cloneCall = "gene",
  y.axes = c("sample", "CellType", "CellType_cluster"),
  color = "CellType"
) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

print(p2)
dev.off()


#graphics.off()

#getCirclize

library(circlize)
library(scales)

circles <- getCirclize(sce, 
                       group.by = "CellType")

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(sce$CellType)))
names(grid.cols) <- unique(sce$CellType)

#Graphing the chord diagram

pdf("figures/chord_diagram_bySample.pdf", width = 6, height = 6)
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()

head(sce$sample)

#TCR1_Control
# Step 1: Subset by sample
subset <- sce[, colData(sce)$sample == "TCR1_Control"]
head(subset)

# Step 2: Generate circlize input data
library(scRepertoire)
circles <- getCirclize(subset, group.by = "CellType", proportion = TRUE)

# Step 3: Define colors for each cell type
library(scales)
cell_types <- unique(colData(subset)$CellType_cluster)
grid.cols <- hue_pal()(length(cell_types))
names(grid.cols) <- cell_types

# Step 4: Draw chord diagram
library(circlize)
pdf("figures/chord_diagram_bySample_TCR1.pdf", width = 6, height = 6)
chordDiagram(
  circles,
  self.link = 1,
  grid.col = grid.cols,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow"
)
dev.off()

#conflicts_prefer(dplyr::filter)

#TCR2_Naive_T_cells
# Step 1: Subset by sample
subset <- sce[, colData(sce)$sample == "TCR2_Naive_T_cells"]
head(subset)

# Step 2: Generate circlize input data
library(scRepertoire)
circles <- getCirclize(subset, group.by = "CellType", proportion = TRUE)

# Step 3: Define colors for each cell type
library(scales)
cell_types <- unique(colData(subset)$CellType_cluster)
grid.cols <- hue_pal()(length(cell_types))
names(grid.cols) <- cell_types

# Step 4: Draw chord diagram
library(circlize)
pdf("figures/chord_diagram_bySample_TCR2.pdf", width = 6, height = 6)
chordDiagram(
  circles,
  self.link = 1,
  grid.col = grid.cols,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow"
)
dev.off()



#TCR3_Adaptive_T_cells
# Step 1: Subset by sample
subset <- sce[, colData(sce)$sample == "TCR3_Adaptive_T_cells"]
head(subset)

# Step 2: Generate circlize input data
library(scRepertoire)
circles <- getCirclize(subset, group.by = "CellType", proportion = TRUE)

# Step 3: Define colors for each cell type
library(scales)
cell_types <- unique(colData(subset)$CellType_cluster)
grid.cols <- hue_pal()(length(cell_types))
names(grid.cols) <- cell_types

# Step 4: Draw chord diagram
library(circlize)
pdf("figures/chord_diagram_bySample_TCR3.pdf", width = 7.5, height = 7.5)
chordDiagram(
  circles,
  self.link = 1,
  grid.col = grid.cols,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow"
)
dev.off()


#Quantifying Clonal Bias
#StartracDiversity


# STARTRAC diversity by sample, grouped by CellType
pdf("figures/StartracDiversity_bySample.pdf", width = 4, height = 5)
StartracDiversity(
  sce,
  type = "sample",           # comparison level: sample-wise diversity
  group.by = "CellType_cluster",     # grouping variable in metadata
  exportTable = FALSE        # set TRUE if you want the diversity table
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )

dev.off()


#Calculating a Single Index

# Calculate and plot only the clonal expansion index

pdf("figures/Single_Index_bySample.pdf", width = 4, height = 4.5)
StartracDiversity(sce, 
                  type = "sample", 
                  group.by = "CellType",
                  index = "expa")+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title = element_blank()
  )
dev.off()

#Pairwise Migration Analysis

# # Calculate pairwise migration between tissues

StartracDiversity(sce, 
                  type = "sample", 
                  group.by = "CellType",
                  index = "migr",
                  pairwise = "ident")


#clonalBias

pdf("figures/clonalBias.pdf", width = 7, height = 5)
clonalBias(sce, 
           cloneCall = "aa", 
           split.by = "sample", 
           group.by = "CellType",
           n.boots = 10, 
           min.expand =5)
dev.off()

pdf("figures/clonalBias1.pdf", width = 7, height = 5)
clonalBias(sce, 
           cloneCall = "aa", 
           split.by = "sample", 
           group.by = "CellType_cluster",
           n.boots = 10, 
           min.expand =5)
dev.off()


pdf("figures/clonalBias2.pdf", width = 7, height = 5)
clonalBias(sce, 
           cloneCall = "nt", 
           split.by = "sample", 
           group.by = "CellType",
           n.boots = 10, 
           min.expand =5)
dev.off()

pdf("figures/clonalBias3.pdf", width = 7, height = 5)
clonalBias(sce, 
           cloneCall = "nt", 
           split.by = "sample", 
           group.by = "CellType_cluster",
           n.boots = 10, 
           min.expand =5)
dev.off()


#Builds mouse C7 T-cell gene set
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scater)
  library(escape)
  library(msigdbr)
  library(dplyr)
  library(tibble)
  library(MatrixGenerics)
  library(matrixStats)
  library(BiocParallel)
  library(GSVA)
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(pheatmap)
})

set.seed(1)

# -----------------------------
# 0) Settings + helpers
# -----------------------------
outdir <- "escape_pub_figs"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

sample_col   <- "sample"
celltype_col <- "CellType"

pdf_hq <- function(file, width = 11, height = 8.5) {
  grDevices::cairo_pdf(filename = file, width = width, height = height, onefile = TRUE)
}
safe_dev_off <- function() while (!is.null(dev.list())) grDevices::dev.off()

pal_group <- function(n) grDevices::hcl.colors(n, palette = "Dark 3")
pal_cont  <- function(n = 256) grDevices::hcl.colors(n, palette = "inferno")

pick_expr_assay <- function(sce) {
  an <- assayNames(sce)
  if ("logcounts" %in% an) return("logcounts")
  if ("counts" %in% an) return("counts")
  stop("No 'logcounts' or 'counts' assay found. Provide normalized data.")
}

ensure_umap <- function(sce) {
  if ("UMAP" %in% reducedDimNames(sce)) return(sce)
  message("UMAP not found; computing PCA+UMAP.")
  sce <- scater::runPCA(sce, exprs_values = pick_expr_assay(sce), ncomponents = 30)
  sce <- scater::runUMAP(sce, dimred = "PCA", name = "UMAP")
  sce
}

stopifnot(sample_col %in% colnames(colData(sce)))
stopifnot(celltype_col %in% colnames(colData(sce)))

# -----------------------------
# 1) Load SCE
# -----------------------------
sce <- readRDS("sce_object.rds")
sce <- ensure_umap(sce)
expr_assay <- pick_expr_assay(sce)

# -----------------------------
# 2) Build T-cell gene sets (MSigDB C7, mouse) + overlap filter
# -----------------------------
mouse_c7 <- msigdbr(species = "Mus musculus", category = "C7")

tcell_sets <- mouse_c7 %>%
  dplyr::filter(grepl("T_CELL", gs_name, ignore.case = TRUE)) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::distinct()

tcell_gene_sets <- tcell_sets %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  tibble::deframe()

genes_in_data <- rownames(sce)
tcell_gene_sets2 <- lapply(tcell_gene_sets, function(g) intersect(g, genes_in_data))
tcell_gene_sets2 <- tcell_gene_sets2[lengths(tcell_gene_sets2) >= 5]
if (length(tcell_gene_sets2) < 2) stop("Too few gene sets after overlap filtering.")

# -----------------------------
# 3A) UCell scoring (stable) -> store in altExp
# -----------------------------
register(SerialParam())

ucell_scores <- escape::escape.matrix(
  sce,
  gene.sets = tcell_gene_sets2,
  method = "UCell",
  min.size = 5,
  BPPARAM = BiocParallel::SerialParam()
)

# escape returns CELLS x SETS in your case (59446 x 135)
ucell_cells_by_sets <- as.matrix(ucell_scores)
stopifnot(nrow(ucell_cells_by_sets) == ncol(sce))

ucell_sets_by_cells <- t(ucell_cells_by_sets)  # SETS x CELLS
colnames(ucell_sets_by_cells) <- colnames(sce)
rownames(ucell_sets_by_cells) <- colnames(ucell_cells_by_sets)

sce[["escape.UCell"]] <- SingleCellExperiment(
  assays = list(scores = ucell_sets_by_cells)
)

# -----------------------------
# 3B) Optional: try ssGSEA (new GSVA API) -> store in altExp if succeeds
# -----------------------------
try_ssgsea <- TRUE

if (try_ssgsea && requireNamespace("GSVA", quietly = TRUE)) {
  message("Attempting GSVA ssGSEA (new API). If it fails, we keep UCell only.")
  
  ssgsea_ok <- FALSE
  ssgsea_sets_by_cells <- NULL
  
  ssgsea_ok <- tryCatch({
    expr <- as.matrix(assay(sce, expr_assay))
    storage.mode(expr) <- "double"
    
    rn <- rownames(expr)
    keep <- !is.na(rn) & nzchar(rn) & !duplicated(rn)
    expr <- expr[keep, , drop = FALSE]
    
    gs <- lapply(tcell_gene_sets2, as.character)
    
    param <- GSVA::ssgseaParam(expr, gs, minSize = 5)
    ssgsea_mat <- GSVA::gsva(param)  # SETS x CELLS
    
    # align columns to sce
    ssgsea_sets_by_cells <- as.matrix(ssgsea_mat)
    ssgsea_sets_by_cells <- ssgsea_sets_by_cells[, colnames(sce), drop = FALSE]
    
    TRUE
  }, error = function(e) {
    message("GSVA ssGSEA failed: ", conditionMessage(e))
    FALSE
  })
  
  if (ssgsea_ok && !is.null(ssgsea_sets_by_cells)) {
    sce[["escape.ssGSEA"]] <- SingleCellExperiment(
      assays = list(scores = ssgsea_sets_by_cells)
    )
  }
}

# -----------------------------
# 4) Choose score source and top variable pathways
# -----------------------------
score_source <- if ("escape.ssGSEA" %in% altExpNames(sce)) "escape.ssGSEA" else "escape.UCell"
score_mat <- assay(altExp(sce, score_source), "scores")  # SETS x CELLS

path_var <- matrixStats::rowVars(score_mat)
top_paths <- names(sort(path_var, decreasing = TRUE))[1:min(12, length(path_var))]
message("Scoring source: ", score_source)
message("Top pathways: ", paste(head(top_paths, 6), collapse = ", "), " ...")

# -----------------------------
# 5) Publication-level figures -> one PDF
# -----------------------------
pdf_file <- file.path(outdir, paste0("ESCAPE_", score_source, "_publication_figures.pdf"))
pdf_hq(pdf_file, width = 11, height = 8.5)

# ---------- Prep UMAP df ----------
df_umap <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(df_umap) <- c("UMAP1", "UMAP2")
df_umap[[sample_col]]   <- as.factor(colData(sce)[[sample_col]])
df_umap[[celltype_col]] <- as.factor(colData(sce)[[celltype_col]])

# ---------- FIG 1: UMAP by sample and cell type ----------
p1a <- ggplot(df_umap, aes(UMAP1, UMAP2, color = .data[[sample_col]])) +
  geom_point(size = 0.15, alpha = 0.7) +
  theme_classic(base_size = 13) +
  labs(title = "UMAP by sample", color = sample_col) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_color_manual(values = pal_group(nlevels(df_umap[[sample_col]])))

p1b <- ggplot(df_umap, aes(UMAP1, UMAP2, color = .data[[celltype_col]])) +
  geom_point(size = 0.15, alpha = 0.7) +
  theme_classic(base_size = 13) +
  labs(title = "UMAP by cell type", color = celltype_col) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_color_manual(values = pal_group(nlevels(df_umap[[celltype_col]])))

print(p1a + p1b + patchwork::plot_layout(guides = "collect"))

# ---------- FIG 2: UMAP feature maps (top 6 pathways) ----------
top6 <- top_paths[1:min(6, length(top_paths))]
feature_plots <- lapply(top6, function(p) {
  df <- df_umap
  df$score <- as.numeric(score_mat[p, ])
  ggplot(df, aes(UMAP1, UMAP2, color = score)) +
    geom_point(size = 0.15, alpha = 0.85) +
    theme_classic(base_size = 12) +
    labs(title = p, color = score_source) +
    scale_color_gradientn(colours = pal_cont())
})
print(patchwork::wrap_plots(feature_plots, ncol = 3))

# ---------- FIG 3: Violin+box by sample (top 3 pathways) ----------
top3 <- top_paths[1:min(3, length(top_paths))]
for (p in top3) {
  df <- data.frame(
    score  = as.numeric(score_mat[p, ]),
    sample = as.factor(colData(sce)[[sample_col]])
  )
  
  p_violin <- ggplot(df, aes(sample, score, fill = sample)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.6, linewidth = 0.2) +
    theme_classic(base_size = 13) +
    labs(title = paste0(p, " (", score_source, ")"), x = sample_col, y = "Enrichment score") +
    scale_fill_manual(values = pal_group(nlevels(df$sample))) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  print(p_violin)
}

# ---------- FIG 4: Heatmap mean score per CellType (top 12) ----------
ct <- as.factor(colData(sce)[[celltype_col]])
mean_by_ct <- sapply(levels(ct), function(level) {
  cells <- which(ct == level)
  if (length(cells) == 0) return(rep(NA_real_, length(top_paths)))
  MatrixGenerics::rowMeans(score_mat[top_paths, cells, drop = FALSE])
})
rownames(mean_by_ct) <- top_paths

zmat <- t(scale(t(mean_by_ct)))
zmat[is.na(zmat)] <- 0

pheatmap::pheatmap(
  zmat,
  border_color = NA,
  clustering_method = "ward.D2",
  main = paste0("Pathway activity by cell type (", score_source, ", row z-score)"),
  fontsize_row = 8,
  fontsize_col = 9
)

# ---------- FIG 5: Correlation among pathways (top 12) ----------
cor_mat <- stats::cor(t(score_mat[top_paths, , drop = FALSE]), method = "spearman")
pheatmap::pheatmap(
  cor_mat,
  border_color = NA,
  clustering_method = "ward.D2",
  main = "Spearman correlation among top pathways",
  fontsize_row = 7,
  fontsize_col = 7
)

# ---------- FIG 6: Density distributions by sample (top 3) ----------
for (p in top3) {
  df <- data.frame(
    score  = as.numeric(score_mat[p, ]),
    sample = as.factor(colData(sce)[[sample_col]])
  )
  dens <- ggplot(df, aes(score, fill = sample)) +
    geom_density(alpha = 0.35, linewidth = 0.2) +
    theme_classic(base_size = 13) +
    labs(title = paste0("Score distribution: ", p), x = "Enrichment score", y = "Density") +
    scale_fill_manual(values = pal_group(nlevels(df$sample)))
  print(dens)
}

safe_dev_off()
message("Done. Wrote PDF: ", pdf_file)


################################################################################
#Clustering by Edit Distance
#clonalCluster: Cluster by Sequence Similarity

load("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/combined_TCR.RData")

head(combined_TCR[[1]])
head(combined_TCR[[2]])
head(combined_TCR[[3]])

head(combined_TCR[[1]]$barcode)
head(combined_TCR[[2]]$barcode)
head(combined_TCR[[3]]$barcode)

## ================================
## 0. Libraries
## ================================
library(tcrpheno)     # or scRepertoire, whichever has clonalCluster()
library(dplyr)
library(tibble)
library(ggplot2)
library(viridisLite)
library(ggrepel)
library(scRepertoire)
library(Seurat)
library(viridisLite)
library(ggrepel)
library(purrr)


## ================================
## 1. Load data
## ================================
# TCR metadata list
load("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/combined_TCR.RData")
# -> object: combined_TCR

# Seurat object
seurat_obj <- readRDS("processed_seurat_object_final.rds")
#saveRDS(seurat_obj, file = "TCR_combined_seurat_onj.rds")

colnames(seurat_obj@meta.data)
# should include: sample, CellType, barcode, UMAP_1, UMAP_2, etc.


## ================================
## 2. Run clonal clustering for TRA
##    (grouped by sample)
## ================================
scRep_example <- clonalCluster(
  combined_TCR,
  chain     = "TRA",
  sequence  = "aa",
  threshold = 0.85,
  group.by  = "sample"
)

# Quick inspect
head(scRep_example[[1]][, c("barcode", "TCR1", "TRA.Cluster")])
table(scRep_example[[1]]$TRA.Cluster)
table(scRep_example[[2]]$TRA.Cluster)
table(scRep_example[[3]]$TRA.Cluster)

#conflicts_prefer(dplyr::filter)

## ================================
## 3. Collapse list -> one TCR table
##    and keep barcode + TRA.Cluster
## ================================
tcr_meta <- do.call(rbind, scRep_example)

tcr_meta_small <- tcr_meta %>%
  select(barcode, TRA.Cluster) %>%
  filter(!is.na(TRA.Cluster))  # keep only labeled clones

head(tcr_meta_small)


## ================================
## 4. Join TRA.Cluster into Seurat
##    (remove old TRA.Cluster first)
## ================================
meta_df <- seurat_obj@meta.data %>%
  rownames_to_column("cell") %>%
  select(-any_of("TRA.Cluster"))  # drop old column if present

meta_joined <- meta_df %>%
  left_join(tcr_meta_small, by = "barcode") %>%
  column_to_rownames("cell")

head(meta_joined)

seurat_obj@meta.data <- meta_joined

# Check labeled vs unlabeled cells
table(seurat_obj$TRA.Cluster, useNA = "ifany")


## ================================
## 5. Prepare data for plotting
## ================================
meta <- seurat_obj@meta.data %>%
  mutate(cell = rownames(seurat_obj@meta.data))

# Background: all cells
bg_df <- meta
head(bg_df)

# Foreground: only cells with TRA.Cluster
fg_df <- meta %>%
  filter(!is.na(TRA.Cluster))
head(fg_df)

table(fg_df$TRA.Cluster)

# Colors for TRA clusters (deep viridis)
valid_clusters <- sort(unique(fg_df$TRA.Cluster))
num_clusters   <- length(valid_clusters)
head(num_clusters)

cluster_colors <- setNames(
  viridis(num_clusters, option = "C"),   # saturated palette
  valid_clusters
)

# Shapes for samples to show sample legend
valid_samples <- sort(unique(fg_df$sample))
shape_values  <- setNames(
  c(16, 17, 15)[seq_along(valid_samples)],  # one shape per sample
  valid_samples
)

# Cell-type label positions (centroid of each CellType on UMAP)
celltype_centroids <- bg_df %>%
  filter(!is.na(CellType)) %>%
  group_by(CellType) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )


## ================================
## 6. Final UMAP plot
## ================================
p_tra_umap <- ggplot() +
  # 6.1 Background: all cells in light grey
  geom_point(
    data  = bg_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey92",
    size  = 0.2,
    alpha = 0.3,
    stroke = 0
  ) +
  # 6.2 Foreground: TRA clusters, with shape = sample
  geom_point(
    data  = fg_df,
    aes(x = UMAP_1,
        y = UMAP_2,
        color = TRA.Cluster,
        shape = sample),
    size  = 1.8,
    alpha = 1,
    stroke = 0
  ) +
  # 6.3 Cell-type labels on UMAP
  geom_text_repel(
    data  = celltype_centroids,
    aes(x = UMAP_1, y = UMAP_2, label = CellType),
    size         = 3,
    color        = "black",
    max.overlaps = Inf,
    box.padding   = 0.3,
    point.padding = 0.2
  ) +
  # 6.4 Scales & legends
  scale_color_manual(values = cluster_colors) +
  scale_shape_manual(values = shape_values) +
  guides(
    color = guide_legend(order = 1, title = "TRA.Cluster"),
    shape = guide_legend(order = 2, title = "Sample")
  ) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.title      = element_blank(),
    legend.position = "right"
  )

pdf("figures/TCR_scSeu_obj_umap_combined_gene.pdf", width = 6, height = 4)
p_tra_umap
dev.off()


#Motif → cell-type tables

## ============================================================
## 0. Setup
## ============================================================
setwd("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/scRNA_seq/")

library(tcrpheno)     # or scRepertoire – clonalCluster() provider
library(scRepertoire)
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(viridisLite)
library(ggrepel)
library(tidyr)        # <-- NEW: for pivot_longer()

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

## ============================================================
## 1. Load data
## ============================================================
load("/Users/madhusaddala/Documents/MaSu/Ganji/TCR-seq/Data/combined_TCR.RData")
seurat_obj <- readRDS("processed_seurat_object_final.rds")

## ============================================================
## 2. TRA clonal clustering (by sample)
## ============================================================
scRep_example <- clonalCluster(
  combined_TCR,
  chain     = "TRA",
  sequence  = "aa",
  threshold = 0.85,
  group.by  = "sample"
)

## ============================================================
## 3. Collapse clonalCluster output & merge into Seurat
## ============================================================
tcr_meta <- do.call(rbind, scRep_example)

tcr_meta_small <- tcr_meta %>%
  select(barcode, TRA.Cluster) %>%
  filter(!is.na(TRA.Cluster)) %>%
  distinct()

meta_df <- seurat_obj@meta.data %>%
  rownames_to_column("cell") %>%
  select(-any_of("TRA.Cluster"))

meta_joined <- meta_df %>%
  left_join(tcr_meta_small, by = "barcode") %>%
  column_to_rownames("cell")

seurat_obj@meta.data <- meta_joined


## ============================================================
## 4. Prepare data for TRA UMAP
## ============================================================
meta <- seurat_obj@meta.data %>%
  mutate(cell = rownames(seurat_obj@meta.data))

bg_df <- meta
fg_df <- meta %>% filter(!is.na(TRA.Cluster))

celltype_centroids <- bg_df %>%
  filter(!is.na(CellType)) %>%
  group_by(CellType) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

valid_clusters <- sort(unique(fg_df$TRA.Cluster))
cluster_colors <- setNames(
  viridis(length(valid_clusters), option = "C"),
  valid_clusters
)

valid_samples <- sort(unique(fg_df$sample))
shape_values  <- setNames(
  c(16, 17, 15)[seq_along(valid_samples)],
  valid_samples
)

## ============================================================
## 5. TRA UMAP (same as your final figure)
## ============================================================
p_tra_umap <- ggplot() +
  geom_point(
    data  = bg_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey92",
    size  = 0.2,
    alpha = 0.3,
    stroke = 0
  ) +
  geom_point(
    data  = fg_df,
    aes(x = UMAP_1,
        y = UMAP_2,
        color = TRA.Cluster,
        shape = sample),
    size  = 1.8,
    alpha = 1,
    stroke = 0
  ) +
  geom_text_repel(
    data  = celltype_centroids,
    aes(x = UMAP_1, y = UMAP_2, label = CellType),
    size         = 3,
    color        = "black",
    max.overlaps = Inf,
    box.padding   = 0.3,
    point.padding = 0.2
  ) +
  scale_color_manual(values = cluster_colors) +
  scale_shape_manual(values = shape_values) +
  guides(
    color = guide_legend(order = 1, title = "TRA.Cluster"),
    shape = guide_legend(order = 2, title = "Sample")
  ) +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.title      = element_blank(),
    legend.position = "right"
  )

pdf("figures/TCR_TRACluster_UMAP.pdf", width = 6, height = 4)
print(p_tra_umap)
dev.off()

## ============================================================
## 6. Build per-cell TCR + CellType metadata (for motifs)
## ============================================================
# Collapse all VDJ into one table
vdj_all <- bind_rows(combined_TCR)

# Make a long table with one row per chain per cell:
# cdr3_aa1 (usually TRA) and cdr3_aa2 (usually TRB)
vdj_long <- vdj_all %>%
  select(barcode, sample, cdr3_aa1, cdr3_aa2) %>%
  pivot_longer(
    cols      = c(cdr3_aa1, cdr3_aa2),
    names_to  = "chain",
    values_to = "CDR3.aa"
  ) %>%
  filter(!is.na(CDR3.aa), CDR3.aa != "")

# Join with Seurat metadata
tcr_cell_meta <- seurat_obj@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, barcode, sample, CellType, CellType_cluster,
         TRA.Cluster, UMAP_1, UMAP_2) %>%
  left_join(vdj_long, by = c("barcode", "sample"))

## ============================================================
## 7. Motifs of interest
## ============================================================
motifs <- c(
  "CASSLDRVEQYF",
  "CATAGSGGKLTL",
  "CAWSLVQSGEQYF",
  "CAVSKSTNTGKLTF",
  "CASSKEGGRDEQYF"
)

## ============================================================
## 8. Motif summary table (motif × cell type × sample × clone)
## ============================================================
motif_table <- tcr_cell_meta %>%
  filter(CDR3.aa %in% motifs) %>%
  group_by(CDR3.aa, sample, CellType, CellType_cluster,
           TRA.Cluster, chain) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(CDR3.aa, desc(n_cells))

motif_table

write.csv(
  motif_table,
  file = "tables/TCR_motif_CellType_summary.csv",
  row.names = FALSE
)

## ============================================================
## 9. Detailed per-cell motif table
## ============================================================
motif_cells_long <- tcr_cell_meta %>%
  filter(CDR3.aa %in% motifs) %>%
  arrange(CDR3.aa, CellType, sample)

motif_cells_long

write.csv(
  motif_cells_long,
  file = "tables/TCR_motif_per_cell_details.csv",
  row.names = FALSE
)

## ============================================================
## 10. Per-motif UMAP figures WITH CELL-TYPE LABELS
## ============================================================
library(ggrepel)

for (motif in motifs) {
  
  # Background vs motif cells
  df <- tcr_cell_meta %>%
    mutate(flag = ifelse(CDR3.aa == motif, "Motif", "Other"))
  
  # Cell-type centroids *only for motif cells*
  motif_celltype_centroids <- tcr_cell_meta %>%
    filter(CDR3.aa == motif, !is.na(CellType)) %>%
    group_by(CellType) %>%
    summarise(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_motif <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    # grey background cells
    geom_point(
      data  = subset(df, flag == "Other"),
      color = "grey92",
      size  = 0.25,
      alpha = 0.3
    ) +
    # red motif cells
    geom_point(
      data  = subset(df, flag == "Motif"),
      color = "red",
      size  = 1.5,
      alpha = 0.9
    ) +
    # cell-type labels for motif cells
    geom_text_repel(
      data  = motif_celltype_centroids,
      aes(x = UMAP_1, y = UMAP_2, label = CellType),
      size         = 3,
      color        = "black",
      max.overlaps = Inf,
      box.padding   = 0.3,
      point.padding = 0.2
    ) +
    theme_bw() +
    ggtitle(paste("CDR3 motif:", motif)) +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(
    filename = paste0("figures/UMAP_motif_", motif, "_withCellTypes.pdf"),
    plot     = p_motif,
    width    = 4,
    height   = 4
  )
}


#############################################
#Single-cell Gene Set Enrichment Analysis

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(escape)
library(BiocParallel)
library(ggplot2)

## ================================
## 0. Libraries & input objects
## ================================
library(SingleCellExperiment)
library(Seurat)
library(escape)
library(BiocParallel)
library(dplyr)
library(ggplot2)

## If needed:
## sce <- readRDS("sce_object.rds")

## Create Seurat object from SCE (if not already done)
seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")

## Use CellType as active identity
Idents(seurat_obj) <- seurat_obj$CellType

## Quick check that sample / CellType exist
stopifnot(all(c("CellType", "sample") %in% colnames(seurat_obj[[]])))


## ================================
## 1. Define mouse T cell gene sets
## ================================
## Custom gene sets in mouse symbol format (e.g. Foxp3, Cd8a, Gzmb)
gene.sets.custom <- list(
  Cytolytic      = c("Gzmb", "Gzmk", "Prf1", "Ifng", "Nkg7"),
  Exhaustion     = c("Pdcd1", "Lag3", "Tigit", "Ctla4", "Havcr2"),
  Treg_signature = c("Foxp3", "Ctla4", "Ikzf2", "Il2ra"),
  CD8_signature  = c("Cd8a", "Cd8b1", "Gzmb", "Cxcr3"),
  CD4_signature  = c("Cd4", "Cxcr4", "Il7r", "Sell")
)

gene.sets <- gene.sets.custom   # keep it simple and species-correct

## Optional: sanity check overlap with your expression matrix
genes_in_data <- rownames(seurat_obj)
overlap_counts <- lapply(gene.sets, function(gs) length(intersect(gs, genes_in_data)))
print(overlap_counts)


## ================================
## 2. Run escape with GSVA (no UCell)
## ================================
## Use SerialParam so there is NO cluster / serialize error
bp_serial <- BiocParallel::SerialParam()

## Run GSVA enrichment per cell
seurat_obj <- runEscape(
  seurat_obj,
  method         = "GSVA",       # << ONLY GSVA, no UCell
  gene.sets      = gene.sets,
  groups         = 1000,
  min.size       = 3,
  new.assay.name = "escape.GSVA",
  BPPARAM        = bp_serial
)

## ================================
## 3. (Optional) Normalize enrichment
## ================================
## This rescales each pathway to make scores positive / comparable
seurat_obj <- performNormalization(
  input.data    = seurat_obj,
  assay         = "escape.GSVA",
  gene.sets     = gene.sets,
  make.positive = TRUE
)

## Normalized assay usually named "escape.GSVA_normalized"
DefaultAssay(seurat_obj) <- "escape.GSVA_normalized"


## ================================
## 4. Heatmap: pathway × CellType (and sample)
## ================================
## Simple heatmap by CellType
p_heat_celltype <- heatmapEnrichment(
  seurat_obj,
  group.by       = "CellType",
  gene.set.use   = "all",                     # all gene sets in gene.sets
  assay          = "escape.GSVA_normalized",
  scale          = TRUE,
  cluster.rows   = TRUE,
  cluster.columns= TRUE,
  palette        = "Spectral"
)

#dir.create("figures", showWarnings = FALSE)

# Modify heatmap plot before saving
p_heat_celltype2 <- p_heat_celltype +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title      = element_blank()
  )

# Save to PDF
pdf("figures/escape_GSVA_heatmap_CellType_mouse_Tcells.pdf",
    width = 6, height = 3)
print(p_heat_celltype2)
dev.off()


## Heatmap, stratified by sample (optional)
p_heat_celltype_sample <- heatmapEnrichment(
  seurat_obj,
  group.by       = "CellType",
  facet.by       = "sample",
  gene.set.use   = "all",
  assay          = "escape.GSVA_normalized",
  scale          = TRUE,
  palette        = "Spectral"
)

# Add the theme BEFORE saving
p_heat_celltype_sample2 <- p_heat_celltype_sample +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title      = element_blank()
  )

# Save the updated plot
pdf(
  "figures/escape_GSVA_heatmap_CellType_by_sample_mouse_Tcells.pdf",
  width = 10,
  height = 5
)
print(p_heat_celltype_sample2)
dev.off()


## ================================
## 5. Geyser / ridge / split plots
## ================================
## Choose representative gene sets
target_cyt <- "Cytolytic"
target_treg <- "Treg_signature"

## 5a. Geyser plot: Cytolytic score by CellType
gey_cyt <- geyserEnrichment(
  seurat_obj,
  assay    = "escape.GSVA_normalized",
  gene.set = target_cyt,
  group.by = "CellType",
  order.by = "mean"
)

# Modify geyser plot with theme improvements
gey_cyt2 <- gey_cyt +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title      = element_blank()
  )

# Save updated plot
pdf("figures/escape_geyser_Cytolytic_by_CellType.pdf",
    width = 6, height = 5)
print(gey_cyt2)
dev.off()


## 5b. Ridge plot: Treg_signature score by CellType
# Use the correct gene set name (with dash)
target_treg <- "Treg-signature"

# Make sure we're on the right assay
DefaultAssay(seurat_obj) <- "escape.GSVA_normalized"

gsva_mat <- GetAssayData(seurat_obj, assay = "escape.GSVA_normalized")
stopifnot(target_treg %in% rownames(gsva_mat))

# Ridge plot for Treg-signature by CellType
ridge_treg <- ridgeEnrichment(
  seurat_obj,
  assay    = "escape.GSVA_normalized",
  gene.set = target_treg,
  group.by = "CellType",
  add.rug  = TRUE,
  scale    = TRUE
)

# Optional styling
ridge_treg2 <- ridge_treg +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title      = element_blank()
  )

# Save to PDF
pdf("figures/escape_ridge_Treg_by_CellType.pdf",
    width = 7, height = 5)
print(ridge_treg2)
dev.off()

#If you later want ridge plots for other sets, just change:
#target_treg <- "Cytolytic"
# or "Exhaustion", "CD8-signature", "CD4-signature"



## 5c. Split plot: Cytolytic by CellType, split by sample
split_cyt <- splitEnrichment(
  seurat_obj,
  assay    = "escape.GSVA_normalized",
  gene.set = target_cyt,
  group.by = "CellType",
  split.by = "sample"
)

# Apply your preferred styling to the split Cytolytic plot
split_cyt2 <- split_cyt +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title      = element_blank()
  )

# Save to PDF
pdf("figures/escape_split_Cytolytic_CellType_sample.pdf",
    width = 7, height = 4)
print(split_cyt2)
dev.off()



## ================================
## 6. PCA on enrichment scores
## ================================
seurat_obj <- performPCA(
  seurat_obj,
  assay  = "escape.GSVA_normalized",
  n.dim  = 1:10,
  dimRed = "escape.PCA"
)

p_pca <- pcaEnrichment(
  seurat_obj,
  dimRed                    = "escape.PCA",
  x.axis                    = "PC1",
  y.axis                    = "PC2",
  group.by                  = "CellType",
  add.percent.contribution  = TRUE,
  display.factors           = TRUE,
  number.of.factors         = 10
)

pdf("figures/escape_PCA_CellType_mouse_Tcells.pdf",
    width = 7, height = 6)
print(p_pca)
dev.off()

#Motif → CellType summary table + plots
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

sce <- readRDS("sce_object.rds")

head(sce)
head(colData(sce))


if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
if (!dir.exists("tables"))  dir.create("tables",  recursive = TRUE)

motifs <- c(
  "CASSLDRVEQYF",
  "CATAGSGGKLTL",
  "CAWSLVQSGEQYF",
  "CAVSKSTNTGKLTF",
  "CASSKEGGRDEQYF"
)

# ---- set your CDR3 amino-acid column here ----
cdr3_col <- "CTaa"   # <- CHANGE if needed (e.g., "aa", "cdr3_aa", etc.)

stopifnot("CellType" %in% colnames(colData(sce)))
stopifnot(cdr3_col %in% colnames(colData(sce)))

md <- as.data.frame(colData(sce)) %>%
  mutate(
    CDR3_AA = as.character(.data[[cdr3_col]]),
    CellType = as.character(.data$CellType),
    sample  = if ("sample" %in% colnames(.)) as.character(.data$sample) else "sample1"
  ) %>%
  filter(!is.na(CDR3_AA), CDR3_AA != "", !is.na(CellType), CellType != "")

# Expand: one row per (cell × motif) for cells that contain that motif
df_hits <- tidyr::expand_grid(cell_id = seq_len(nrow(md)), motif = motifs) %>%
  mutate(
    CDR3_AA  = md$CDR3_AA[cell_id],
    CellType = md$CellType[cell_id],
    sample   = md$sample[cell_id],
    hit      = str_detect(CDR3_AA, fixed(motif))
  ) %>%
  filter(hit) %>%
  select(motif, sample, CellType)

#conflicts_prefer(dplyr::count)

# ---- Summary table: counts + within-motif proportions ----
tab_motif_celltype <- df_hits %>%
  count(motif, CellType, name = "n_cells") %>%
  group_by(motif) %>%
  mutate(
    motif_total = sum(n_cells),
    prop_within_motif = n_cells / motif_total
  ) %>%
  ungroup() %>%
  arrange(motif, desc(n_cells))

write.csv(tab_motif_celltype, "tables/motif_by_CellType_counts_new.csv", row.names = FALSE)

# ---- Summary table: counts + within-celltype proportions (optional) ----
tab_celltype_motif <- df_hits %>%
  count(CellType, motif, name = "n_cells") %>%
  group_by(CellType) %>%
  mutate(
    celltype_total = sum(n_cells),
    prop_within_celltype = n_cells / celltype_total
  ) %>%
  ungroup() %>%
  arrange(CellType, desc(n_cells))

write.csv(tab_celltype_motif, "tables/CellType_by_motif_counts_new.csv", row.names = FALSE)

#Plot 1 — Stacked barplot: for each motif, which CellTypes contribute?
# Optional: keep CellType order you use in barplots
celltype_levels <- c("cnkCD8","cT","efCD4","efCD8","exCD8","iNKT","mCD8",
                     "nCD4","nCD4/CD8","nCD8","nkCD8","pT","TregCD4")
tab_motif_celltype$CellType <- factor(tab_motif_celltype$CellType, levels = celltype_levels)

p_bar <- ggplot(tab_motif_celltype,
                aes(x = motif, y = prop_within_motif, fill = CellType)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(x = NULL, y = "CellType proportion within motif", fill = "CellType")

grDevices::cairo_pdf("figures/motif_CellType_stackedBar.pdf", width = 5, height = 4)
print(p_bar)
dev.off()


#Plot 2 — Heatmap-like tile: motif × CellType (proportions)
p_heat <- ggplot(tab_motif_celltype,
                 aes(x = CellType, y = motif, fill = prop_within_motif)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "grey90", high = "darkred",
                      labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(fill = "Proportion")

grDevices::cairo_pdf("figures/motif_by_CellType_heatmap.pdf", width = 6, height = 3)
print(p_heat)
dev.off()


#Plot 3 — Dot plot: size = #cells, color = proportion
p_dot <- ggplot(tab_motif_celltype,
                aes(x = CellType, y = motif)) +
  geom_point(aes(size = n_cells, color = prop_within_motif), alpha = 0.95) +
  scale_color_gradient(low = "grey70", high = "darkred",
                       labels = percent_format(accuracy = 1)) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(size = "# cells", color = "Proportion")

grDevices::cairo_pdf("figures/motif_by_CellType_dotplot.pdf", width = 4, height = 4)
print(p_dot)
dev.off()



