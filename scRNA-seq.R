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


#Load scRNA-seq Data and Create Seurat Object
# Define scRNA-seq file paths
sample_paths <- list(
  "TCR1-Control" = "/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/scRNA-seq/cellranget_out_transcriptome/TCR1-Control/outs/filtered_feature_bc_matrix",
  "TCR2-Naive_T_cells" = "/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/scRNA-seq/cellranget_out_transcriptome/TCR2-Naive_T_cells/outs/filtered_feature_bc_matrix",
  "TCR3-Adaptive_T_cells" = "/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/scRNA-seq/cellranget_out_transcriptome/TCR3-Adaptive_T_cells/outs/filtered_feature_bc_matrix"
)

# Load scRNA-seq data into Seurat objects
seurat_objects <- lapply(names(sample_paths), function(sample) {
  data <- Read10X(sample_paths[[sample]])
  seurat_obj <- CreateSeuratObject(counts = data, project = sample, min.cells = 3, min.features = 100)
  seurat_obj$sample <- sample
  return(seurat_obj)
})

# Merge Seurat objects into one dataset
seurat_merged <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = names(seurat_objects), project = "scRNA_TCR_Seq")


#Quality Control & Filtering
# Plot QC metrics
seurat_obj <- readRDS(file = "/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/scRNA-seq/seurat_merged_raw.rds")
head(seurat_obj)



# calculate the percentage of mitochondrial reads per cell
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-") #mt for mouse and rat and MT for humans only.


#Sample wise filtering 
#Step 1: Define Sample-Specific Thresholds
# Assuming sample metadata is stored in 'seurat_obj$sample'
# plot distributions of QC metrics, grouped by SampleID
pdf('figures/basic_qc.pdf', width=8, height=10)
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='sample', raster=FALSE,
  ncol = 2, pt.size=0)
dev.off()

table(seurat_obj@meta.data$sample)

#Step 2: Compute Sample-Specific Quantiles
#Calculate quantiles for each sample to define appropriate filtering thresholds dynamically.
# Function to compute filtering thresholds per sample
compute_thresholds <- function(seurat_obj, sample_col, feature) {
  return(seurat_obj@meta.data %>%
           group_by(!!sym(sample_col)) %>%
           summarize(
             min_val = quantile(!!sym(feature), 0.05),  # 5th percentile
             max_val = quantile(!!sym(feature), 0.95)   # 95th percentile
           ))
}


# Compute thresholds for each sample
nFeature_thresh <- compute_thresholds(seurat_obj, "sample", "nFeature_RNA")
nCount_thresh <- compute_thresholds(seurat_obj, "sample", "nCount_RNA")
percent_mt_thresh <- compute_thresholds(seurat_obj, "sample", "percent.mt")

print(nFeature_thresh)
print(nCount_thresh)
print(percent_mt_thresh)

#Apply Sample-Specific Filtering
# Create a filtered Seurat object
filtered_cells <- c()
#conflicts_prefer(dplyr::filter)
# Loop through each sample and apply sample-specific thresholds
for (sample_id in unique(seurat_obj$sample)) {
  
  # Get thresholds for this sample
  nFeature_min <- nFeature_thresh %>% filter(sample == sample_id) %>% pull(min_val)
  nFeature_max <- nFeature_thresh %>% filter(sample == sample_id) %>% pull(max_val)
  nCount_min <- nCount_thresh %>% filter(sample == sample_id) %>% pull(min_val)
  nCount_max <- nCount_thresh %>% filter(sample == sample_id) %>% pull(max_val)
  percent_mt_max <- percent_mt_thresh %>% filter(sample == sample_id) %>% pull(max_val)
  
  # Subset the cells for this sample
  sample_cells <- subset(seurat_obj, subset = sample == sample_id & 
                                        nFeature_RNA >= nFeature_min & 
                                        nFeature_RNA <= nFeature_max &
                                        nCount_RNA >= nCount_min & 
                                        nCount_RNA <= nCount_max &
                                        percent.mt <= percent_mt_max)
  
  # Collect cell names
  filtered_cells <- c(filtered_cells, colnames(sample_cells))
}

# Keep only the selected cells
seurat_obj_filtered <- subset(seurat_obj, cells = filtered_cells)

# Verify the filtered object
head(seurat_obj_filtered)

saveRDS(seurat_obj_filtered, file="/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/scRNA-seq/filter_seurat_obj.rds")

gc()
rm(list=ls())

seurat_obj <- readRDS("/dfs7/swaruplab/msaddala/test/ganji/TCR_seq/scRNA-seq/filter_seurat_obj.rds")
head(seurat_obj)

# apply filter
#seurat_obj <- subset(seurat_obj, nCount_RNA > 200 & nCount_RNA <= 10000 & 
#nFeature_RNA >= 10 & nFeature_RNA <= 6000 & percent.mt <= 10)

#seurat_obj <- subset(seurat_obj, sample != 'NA')#if any sample not good enough just excludes

# plot the number of cells in each sample post filtering
df <- as.data.frame(rev(table(seurat_obj$sample)))
colnames(df) <- c('sample', 'n_cells')
p <- ggplot(df, aes(y=n_cells, x=reorder(sample, -n_cells), fill=sample)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[cells])) + xlab('sample') +
  ggtitle(paste('Total cells post-filtering:', sum(df$n_cells))) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )
pdf('figures/basic_cells_per_sample_filtered_df.pdf', width=3, height=4)
print(p)
dev.off()

# plot distributions of QC metrics, grouped by SampleID
pdf('figures/basic_qc_filtered.pdf', width=6, height=7)
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='sample', raster=FALSE,
  ncol = 2, pt.size=0)
dev.off()

gc()
rm(list=ls())

# Check the assay structure
str(seurat_obj@assays$RNA)  # Adjust the assay name if it's different

#Join Data Layers
seurat_obj <- JoinLayers(seurat_obj)
# Set the future.globals.maxSize option
#options(future.globals.maxSize = 10 * 1024^3)

# log normalize data
# Normalize the data
# log normalize data

seurat_obj <- NormalizeData(seurat_obj)

#memory.limit(size = 160000)  # Set to available RAM in MB (e.g., 160GB)

# Use future for parallel processing
#plan("multicore", workers = 4)  # Adjust the number of workers accordingly

# Identify the 4000 most highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# scale data:
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

p <- LabelPoints(
  VariableFeaturePlot(seurat_obj),
  points = head(VariableFeatures(seurat_obj),10),
  repel = TRUE
) + theme(legend.position="bottom")

pdf('figures/basic_variable_features.pdf', width=6, height=6)
print(p)
dev.off()


#Join Data Layers
#seurat_obj <- JoinLayers(seurat_obj)

#Linear Dimensionality Reduction
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj),
  npcs=50
)

# plot the top genes contributing to the first 3 PCs
p <- VizDimLoadings(seurat_obj, dims = 1:3, reduction = "pca", ncol=3)

pdf('figures/basic_pca_loadings.pdf', width=10, height=5)
print(p)
dev.off()

# PCA scatter plot colored by SampleID
p <- DimPlot(seurat_obj, reduction = "pca", group.by='sample')

pdf('figures/basic_pca_scatter.pdf', width=7, height=6)
print(p)
dev.off()

# plot heatmaps for first 16 PCs
pdf('figures/basic_pca_heatmap.pdf', width=10, height=10)
DimHeatmap(seurat_obj, dims = 1:16, cells = 500, balanced = TRUE, ncol=4)
dev.off()

pdf('figures/basic_pca_elbow.pdf', width=6, height=3)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

#Clustering and non-linear dimensionality reduction

# KNN and clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2)

# non-linear reductions (UMAP & t-SNE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)


umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)
pdf('figures/basic_umap_clusters.pdf', width=5, height=5)
DimPlot(seurat_obj, reduction = "umap", group.by='seurat_clusters', label.size = 4, label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
dev.off()

#cluster wise
pdf('figures/basic_umap_clusters1.pdf', width=7, height=5)
DimPlot(seurat_obj, reduction = "umap", group.by='seurat_clusters', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

#sample wise
pdf('figures/basic_umap_sample.pdf', width=7, height=5)
DimPlot(seurat_obj, reduction = "umap", group.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

pdf('figures/basic_tsne_clusters.pdf', width=4, height=4)
DimPlot(seurat_obj, reduction = "tsne", group.by='seurat_clusters', label.size = 4, label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
dev.off()

#cluster wise
pdf('figures/basic_tsne_clusters1.pdf', width=5, height=5)
DimPlot(seurat_obj, reduction = "tsne", group.by='seurat_clusters', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

#sample wise
pdf('figures/basic_tsne_sample.pdf', width=6, height=5)
DimPlot(seurat_obj, reduction = "tsne", group.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()


#sample wise split
pdf('figures/basic_tsne_sample_split.pdf', width=12, height=5)
DimPlot(seurat_obj, reduction = "tsne", split.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()


#markers
Naïve_T_cells <- c('Ccr7', 'Sell', 'Il7r', 'Il2rg', 'Cd27')
CD8_cells <- c('Cd8a', 'Cd8b1', 'Cd3d', 'Cd3e', 'Cd3g')
CD4Tregs <- c('Ctla4', 'Il2ra', 'Foxp3')
CD4_cells <- c('Cd4', 'Lef1')
Effector_CD4_T_cells <- c('Il18r1', 'Tnfsf11')
Effector_CD8_T_cells <- c('Gzmb', 'Prf1', 'Klrd1', 'Tbx21', 'Cd69', 'Icos', 'Ccl3', 'Ccl4', 'Ccl5', 'Il2')
NK_cells <- c('Nkg7', 'Gzmm', 'Eomes', 'Gata3')
Exhausted_T_cells <- c("Pdcd1", "Lag3", "Havcr2")

features <- list("nTC" = Naïve_T_cells, "CD8" = CD8_cells, "CD4Tr" = CD4Tregs, "CD4" = CD4_cells, "eCD4" = Effector_CD4_T_cells, "eCD8" = Effector_CD8_T_cells, "NK" = NK_cells, "exTC" = Exhausted_T_cells)
head(features)


pdf('figures/new_markers_final.pdf', width=12, height=4.5)
DotPlot(object = seurat_obj, features=features, dot.scale=7, cols="RdBu", cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
dev.off()



# Harmony batch correction
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "sample", verbose = TRUE)

#Required: re-compute graph & clusters after Harmony
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Now run UMAP & t-SNE
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- RunTSNE(seurat_obj, reduction = "harmony", dims = 1:30)



pdf('figures/harmony_umap_clusters.pdf', width=5, height=5)
DimPlot(seurat_obj, reduction = "umap", group.by='seurat_clusters', label.size = 4, label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
dev.off()

#cluster wise
pdf('figures/harmony_umap_clusters1.pdf', width=5, height=5)
DimPlot(seurat_obj, reduction = "umap", group.by='seurat_clusters', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

#sample wise
pdf('figures/harmony_umap_sample.pdf', width=6, height=5)
DimPlot(seurat_obj, reduction = "umap", group.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

#sample wise split
pdf('figures/harmony_umap_sample_split.pdf', width=12, height=5)
DimPlot(seurat_obj, reduction = "umap", split.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()


#sample wise split
pdf('figures/harmony_umap_sample_split1.pdf', width=12, height=5)
DimPlot(seurat_obj, reduction = "umap", split.by='sample', label.size = 4, label=TRUE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()




pdf('figures/harmony_tsne_clusters.pdf', width=4, height=4)
DimPlot(seurat_obj, reduction = "tsne", group.by='seurat_clusters', label.size = 4, label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
dev.off()

#cluster wise
pdf('figures/harmony_tsne_clusters1.pdf', width=5, height=4)
DimPlot(seurat_obj, reduction = "tsne", group.by='seurat_clusters', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

#sample wise
pdf('figures/harmony_tsne_sample.pdf', width=6, height=4)
DimPlot(seurat_obj, reduction = "tsne", group.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()


#sample wise split
pdf('figures/harmony_tsne_sample_split.pdf', width=12, height=5)
DimPlot(seurat_obj, reduction = "tsne", split.by='sample', label.size = 4, label=FALSE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

#sample wise split
pdf('figures/harmony_tsne_sample_split1.pdf', width=12, height=5)
DimPlot(seurat_obj, reduction = "tsne", split.by='sample', label.size = 4, label=TRUE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()

Finally selected best markers
# Cluster wise
pdf('figures/Final_T-cell_dotPlot.pdf', width=17, height=7)
DotPlot(seurat_obj, features = c("Ccr7","Sell","Tcf7","Lef1","Cd3d","Cd3e","Cd3g",
                                 "Cd4","Cd8a","Cd8b1","Foxp3","Egr2","Nkg7","Gzmm","Gata3",
                                 "Il2ra","Il17rb","Il33","Tbx21","Prf1","Gzma","Gzmb",
                                 "Ifng","Il2","Klrg1","Il18r1","Tnfsf11","Havcr2",
                                 "Cd69","Ccl4","Ccl5","Ccr8","Pdcd1","Cxcr6","Ifngr1",
                                 "Cxcr3","Tnf","Batf","Areg","Ctsh",
                                 "Irf4","Il1r1","Klrb1","Tnfrsf4","Icos","Il21","Lag3",
                                 "Ctla4","Cd28","Ptprc"), dot.scale=9, cols="RdBu") + 
  RotatedAxis()
dev.off()


#Sample wise
pdf('figures/Final_T-cell_dotPlot_sample.pdf', width=17, height=6)
DotPlot(seurat_obj, features = c("Ccr7","Sell","Tcf7","Lef1","Cd3d","Cd3e","Cd3g",
                                 "Cd4","Cd8a","Cd8b1","Foxp3","Egr2","Nkg7","Gzmm","Gata3",
                                 "Il2ra","Il17rb","Il33","Tbx21","Prf1","Gzma","Gzmb",
                                 "Ifng","Il2","Klrg1","Il18r1","Tnfsf11","Havcr2",
                                 "Cd69","Ccl4","Ccl5","Ccr8","Pdcd1","Cxcr6","Ifngr1",
                                 "Cxcr3","Tnf","Batf","Areg","Ctsh",
                                 "Irf4","Il1r1","Klrb1","Tnfrsf4","Icos","Il21","Lag3",
                                 "Ctla4","Cd28","Ptprc"), group.by = "sample", dot.scale=9, cols="RdBu") + 
  RotatedAxis()
dev.off()



cluster_annotations <- list(
  '0' = 'nCD4',
  '1' = 'nCD4',
  '2' = 'nCD8',
  '3' = 'efCD4',
  '4' = 'TregCD4',
  '5' = 'nCD8',
  '6' = 'mCD8',
  '7' = 'TregCD4',
  '8' = 'nCD4',
  '9' = 'nCD4/CD8',
  '10' = 'cT',
  '11' = 'cnkCD8',
  '12' = 'nCD4',
  '13' = 'nCD4',
  '14' = 'exCD8',
  '15' = 'pT',
  '16' = 'cT',
  '17' = 'iNKT',
  '18' = 'nkCD8',
  '19' = 'efCD8'
)

seurat_obj@meta.data$CellType <- unlist(cluster_annotations[seurat_obj$seurat_clusters])
seurat_obj@meta.data$CellType_cluster <- paste0(seurat_obj$CellType, '-', seurat_obj$seurat_clusters)

pdf('figures/basic_umap_celltypes.pdf', width=11, height=8)
DimPlot(seurat_obj, reduction = "umap", group.by='CellType') +
  umap_theme + ggtitle('UMAP colored by cell type annotations')
dev.off()

pdf('figures/basic_umap_celltype_clusters.pdf', width=10, height=7)
DimPlot(seurat_obj, reduction = "umap", group.by='CellType_cluster', label=TRUE) +
  umap_theme + ggtitle('UMAP colored by cell type + cluster') + NoLegend()
dev.off()

pdf('figures/basic_umap_celltype_sample.pdf', width=10, height=8)
DimPlot(seurat_obj, reduction = "umap", group.by='sample') +
  umap_theme + ggtitle('UMAP colored by cell type + cluster')
dev.off()

#sample wise split
pdf('figures/basic_umap_celltype_sample_split.pdf', width=16, height=6)
DimPlot(seurat_obj, reduction = "umap", split.by='sample', label.size = 4, label=TRUE) +
  umap_theme + ggtitle('UMAP colored by seurat clusters')
dev.off()


pdf('figures/basic_tsne_celltypes.pdf', width=9, height=5)
DimPlot(seurat_obj, reduction = "tsne", group.by='CellType') +
  umap_theme + ggtitle('tSNE colored by cell type annotations')
dev.off()

pdf('figures/basic_tsne_celltype_clusters.pdf', width=8, height=6)
DimPlot(seurat_obj, reduction = "tsne", group.by='CellType_cluster', label=TRUE) +
  umap_theme + ggtitle('tSNE colored by cell type + cluster') + NoLegend()
dev.off()

pdf('figures/basic_tsne_celltype_sample.pdf', width=6, height=5)
DimPlot(seurat_obj, reduction = "tsne", group.by='sample') +
  umap_theme + ggtitle('tSNE colored by cell type + cluster')
dev.off()

#sample wise split
pdf('figures/basic_tsne_celltype_sample_split.pdf', width=16, height=6)
DimPlot(seurat_obj, reduction = "tsne", split.by='sample', label.size = 4, label=TRUE) +
  umap_theme + ggtitle('tSNE colored by seurat clusters')
dev.off()

#Marker genes for cell types
pdf('figures/Final_marker_T-cell_dotPlot.pdf', width=6, height=5)
DotPlot(seurat_obj, features=c("Cd3d", "Cd4", "Cd8a", "Ccr7", "Gzmb", "Nkg7"), dot.scale=9, group.by="CellType", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()
dev.off()


#Marker genes for cell types
pdf('figures/Final_marker_T-cell_dotPlot_new.pdf', width=6, height=5)
DotPlot(seurat_obj, features=c("Cd3d", "Cd4", "Cd8a", "Ccr7", "Gzmb", "Nkg7", "Ifng", "Prf1"), dot.scale=9, group.by="CellType", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()
dev.off()


# Set PDF output with higher resolution
pdf("figures/Final_marker_T-cell_dotPlot_new1.pdf", width = 6, height = 14, useDingb = FALSE)
# Generate dot plot
DotPlot(seurat_obj, 
        features = new_markers, 
        dot.scale = 8, 
        group.by = "CellType", 
        cols = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
  theme_minimal(base_size = 14) +  # Base font size for publication
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90")) +
  coord_flip() +
  RotatedAxis()
dev.off()


#Marker genes for cell types
pdf('figures/Final_marker_T-cell_dotPlot_new1a.pdf', width=6, height=5)
DotPlot(seurat_obj, features=c("Cd3d", "Cd4", "Cd8a", "Ccr7", "Gzmb", "Nkg7", "Ifng", "Prf1", "Klrg1"), dot.scale=9, group.by="CellType", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()
dev.off()



pdf('figures/new_markers_T-cell_dotPlot1b.pdf', width=19, height=7)
DotPlot(seurat_obj, features = new_markers, group.by = "CellType_cluster", dot.scale=9, cols="RdBu") + 
  RotatedAxis()
dev.off()

pdf('figures/new_markers_T-cell_dotPlot2.pdf', width=19, height=7)
DotPlot(seurat_obj, features = new_markers, group.by = "CellType", dot.scale=9, cols="RdBu") + 
  RotatedAxis()
dev.off()

saveRDS(seurat_obj, file='processed_seurat_object_final.rds')




