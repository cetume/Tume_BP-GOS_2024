#--------------------------------------------------------------------------------------
#
#    Load Herring snRNA-seq data and downsample
#
#--------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(Matrix)

##  Set Variables  --------------------------------------------------------------------
cat('\nLoad Herring snRNA-seq data ... \n')
counts <- toString(snakemake@input[['counts']])
cellMeta <- toString(snakemake@input[['cellMeta']])
geneMeta <- toString(snakemake@input[['geneMeta']])
#seurat_obj <- toString(snakemake@input[['seurat_obj']])
r_dir <- toString(snakemake@params[['r_dir']])

## Set log file  ----------------------------------------------------------------------
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = c("output", "message"))

## Load herring data  -----------------------------------------------------------------
counts <- readMM(counts)
cellMeta <- read.csv(cellMeta)
geneMeta <- read.csv(geneMeta)

# Set the rownames and colnames of matrix
rownames(counts)<-cellMeta$Barcode
colnames(counts)<-geneMeta$GeneName

seurat_herring <- CreateSeuratObject(counts = t(counts))

# Set the meta data
seurat_herring@meta.data<-cbind(cellMeta,seurat_herring@meta.data)
rownames(seurat_herring@meta.data)<-colnames(seurat_herring)

# Remove data frames no longer need
rm(counts, geneMeta, cellMeta)

## Filter out cell populations with <50 cells  ----------------------------------------

# Cell counts
cell_counts_lvl_2 <- as.data.frame(table(seurat_herring@meta.data$sub_clust)) 
colnames(cell_counts_lvl_2) <- c('sub_clust', 'cell_count_per_cluster')

# Add meta.data column - number of cells for each cell population (sub_clust)
sub_clusters <- as.data.frame(seurat_herring@meta.data$sub_clust)
colnames(sub_clusters) <- 'sub_clust'
sub_clusters_cell_count <- sub_clusters %>% inner_join(cell_counts_lvl_2, by = 'sub_clust')
sub_clusters_cell_count$cell_count_per_cluster <- as.numeric(as.character(sub_clusters_cell_count$cell_count_per_cluster))
as.data.frame(sub_clusters_cell_count)
rownames(sub_clusters_cell_count) = rownames(seurat_herring@meta.data)
seurat_herring <- AddMetaData(seurat_herring, sub_clusters_cell_count)

# Filter cell populations
seurat_herring <- subset(seurat_herring, subset = cell_count_per_cluster >= 50) #2 cell populations removed (Oligo-6 and Oligo-7)

## Downsample data - standardise cell count across cell populations -------------------

set.seed(123)

# Sub clusters - downsampled to 300 cells - 76/84 have >300 cells
seurat_herring_dwnSmpl_lvl2 <- seurat_herring
Idents(seurat_herring_dwnSmpl_lvl2) <- seurat_herring_dwnSmpl_lvl2@meta.data$sub_clust
seurat_herring_dwnSmpl_lvl2  <- subset(seurat_herring_dwnSmpl_lvl2, downsample = 300)
table(seurat_herring_dwnSmpl_lvl2@meta.data$sub_clust)

## Save RDS objects  ------------------------------------------------------------------
saveRDS(object = seurat_herring_dwnSmpl_lvl2, paste0(r_dir, 'seurat_herring_dwnSmpl_lvl2.rds'))
saveRDS(object = seurat_herring, paste0(r_dir, 'seurat_herring.rds'))

