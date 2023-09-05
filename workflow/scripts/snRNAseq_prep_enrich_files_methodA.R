#--------------------------------------------------------------------------------------
#
#    Prep enrichment test input files for MAGMA and SLDSR     
#
#--------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------
library(tdespec)
library(dplyr)
library(Seurat)
library(EWCE)

## Set variables  ---------------------------------------------------------------------
cat('\nPrepping enrichment test input files for MAGMA and SLDSR ... \n')
seurat_obj <- toString(snakemake@input[['seurat_obj']])
gene_coord <- toString(snakemake@input[['gene_coord']])
mhc_genes <- toString(snakemake@input[['mhc_genes']])
ctd_outdir <- toString(snakemake@params[['ctd_outdir']])
outdir <- toString(snakemake@params[['outdir']])
outfile <- toString(snakemake@output)
study_id <- toString(snakemake@params[['study_id']])
threads <- snakemake@threads

## Set log file  ----------------------------------------------------------------------
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = c("output", "message"))

## Report inputs  ---------------------------------------------------------------------
cat('\nVariables set to: \n\n')
tibble(Variable = c('seurat_obj', 'gene_coord', 'mhc_genes', 'ctd_outdir', 'outdir', 
                    'outfile', 'study_id', 'threads'),
       Value = c(seurat_obj, gene_coord, mhc_genes, ctd_outdir, outdir, outfile, study_id, 
                 threads))

## Load RDS objects  ------------------------------------------------------------------
cat('\nLoad RDS objects ... \n\n')
seurat_obj <- readRDS(seurat_obj)
mhc_genes_obj <- readRDS(mhc_genes) 
gene_coord_obj <- readRDS(gene_coord)

cat('\nSeurat obj loaded: ... \n\n')
seurat_obj
cat('\nSeurat metadata cols: ... \n\n')
glimpse(seurat_obj[[]])
cat('\nMHC obj loaded: ... \n\n')
head(mhc_genes_obj)
cat('\ngene_coord obj loaded: ... \n\n')
head(gene_coord_obj)

## Generate enrichment files for testing  ---------------------------------------------
raw_counts <- seurat_obj@assays$RNA@counts
raw_counts_no_mhc <- raw_counts[!(rownames(raw_counts) %in% mhc_genes_obj), ]
cat('\nMHC genes removed:', dim(raw_counts)[1] - dim(raw_counts_no_mhc)[1])

# Create annotations 
annotations <- as.data.frame(cbind(as.vector(rownames(seurat_obj@meta.data)),
                                   as.vector(seurat_obj$major_clust), 
                                   as.vector(seurat_obj$sub_clust)))
colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
rownames(annotations) <- NULL
annotLevels <- list(level1class = annotations$level1class, 
                    level2class = annotations$level2class)

# Normalize - this is optional, was not used in the original EWCE publication
cat('\nRunning SCT ... ', '\n\n')
options(future.globals.maxSize = 1000 * 1024^2) 
counts_sct <- EWCE::sct_normalize(raw_counts_no_mhc) 

cat('\ncounts_sct class: ', class(counts_sct), '\n\n')

#cat('\nRunning CPM ... ', '\n')
#counts_cpm <- edgeR::cpm(raw_counts_no_mhc)

#cat('\nCheck cpm counts ... ', '\n')
#print(counts_cpm[1:10,1:10])
#print(typeof(counts_cpm))
#counts_cpm <- as.matrix(counts_cpm)
#class(typeof(counts_cpm))

# cat('\nDropping uninformative genes raw ... ', '\n\n')
# DROP_GENES_RAW <- EWCE::drop_uninformative_genes(
#   exp = RAW_COUNTS_NO_MHC, 
#   input_species = "human",
#   output_species = "human",
#   level2annot = annotLevels$level2class) 

cat('\nDropping uninformative genes sct norm ... ', '\n\n')
drop_genes_sct <- EWCE::drop_uninformative_genes(
  exp = counts_sct, 
  input_species = "human",
  output_species = "human",
  level2annot = annotLevels$level2class) 

#cat('\nDropping uninformative genes cpm norm ... ', '\n\n')
#drop_genes_cpm <- EWCE::drop_uninformative_genes(
#  exp = counts_cpm, 
#  input_species = "human",
#  output_species = "human",
#  level2annot = annotLevels$level2class,
#  no_cores = threads) 

cat('\nGene counts:',
    '\n\nRAW:', dim(raw_counts)[1],
    '\nRAW_NO_MHC:', dim(raw_counts_no_mhc)[1],
    '\nCPM_DROP_GENES:', dim(drop_genes_sct)[1])

rm(raw_counts, raw_counts_no_mhc)

rm(counts_sct)

# cat('\nGene counts:',
#     '\n\nRAW:', dim(RAW_COUNTS)[1],
#     '\nRAW_NO_MHC:', dim(RAW_COUNTS_NO_MHC)[1],
#     '\nNO_NORM_DROP_GENES:', dim(DROP_GENES_RAW)[1],
#     '\nSCT_DROP_GENES:', dim( DROP_GENES_SCT)[1],
#     '\nCPM_DROP_GENES:', dim(DROP_GENES_CPM)[1])

# Create object - saves ctd obj to folder
cat('\nCreating CTD object ... \n\n')
ctd_path <- EWCE::generate_celltype_data(exp = drop_genes_sct, 
                                         annotLevels = annotLevels, 
                                         groupName = study_id,
                                         savePath = ctd_outdir,
                                         numberOfBins = 10)


create_enrich_test_files(ctd_path = ctd_path,
                         study_id = study_id,
                         gene_coordinates = gene_coord_obj,
                         levels = 2,
                         outdir = outdir,
                         sldsr_ext = c(100, 100),
                         magma_ext = c(35, 10))

file.create(outfile)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
