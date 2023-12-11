#--------------------------------------------------------------------------------------
#
#    Prep enrichment test input files for MAGMA and SLDSR
#
#--------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------

library(tidyr)
library(dplyr)
library(Seurat)
library(EWCE)
library(scTenifoldNet)
library(ggdendro)
library(readr)

## Set variables  ---------------------------------------------------------------------
cat('\nPrepping enrichment test input files for MAGMA and SLDSR ... \n')
seurat_obj <- toString(snakemake@input[['seurat_obj']])
gene_coord <- toString(snakemake@input[['gene_coord']])
protein_gene_coord <- toString(snakemake@input[['protein_gene_coord']])
mhc_genes <- toString(snakemake@input[['mhc_genes']])
ctd_outdir <- toString(snakemake@params[['ctd_outdir']])
outdir <- toString(snakemake@params[['outdir']])
outfile <- toString(snakemake@output)
study_id <- toString(snakemake@params[['study_id']])
level <- 2
threads <- snakemake@threads

## Set log file  ----------------------------------------------------------------------
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = c("output", "message"))

## Report inputs  ---------------------------------------------------------------------
cat('\nVariables set to: \n\n')
tibble(Variable = c('seurat_obj', 'gene_coord', 'protein_gene_coord', 'mhc_genes', 'ctd_outdir', 'outdir',
                    'outfile', 'study_id', 'threads'),
       Value = c(seurat_obj, gene_coord, protein_gene_coord, mhc_genes, ctd_outdir, outdir, outfile, study_id,
                 threads))

dir.create(ctd_outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(outdir,  recursive = TRUE, showWarnings = FALSE)

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

## Remove genes in MHC region and create annotation levels ----------------------------
raw_counts <- seurat_obj@assays$RNA@counts
raw_counts_no_mhc <- raw_counts[!(rownames(raw_counts) %in% mhc_genes_obj), ]
cat('\nMHC genes removed:', dim(raw_counts)[1] - dim(raw_counts_no_mhc)[1])

# Create annotations
annotations <- as.data.frame(cbind(as.vector(rownames(seurat_obj@meta.data)),
                                   as.vector(seurat_obj$sub_clust)))
colnames(annotations) <- c('cell_id', 'level2class')
rownames(annotations) <- NULL
annotLevels <- list(level1class = annotations$level2class,
                    level2class = annotations$level2class)

## Normalize - this is optional, was not used in the original EWCE publication --------
cat('\nRunning SCT ... ', '\n\n')
options(future.globals.maxSize = 1000 * 1024^2)
counts_sct <- EWCE::sct_normalize(raw_counts_no_mhc)

cat('\ncounts_sct class: ', class(counts_sct), '\n\n')

## Drop uninformative genes and create ctd object (saves to folder)--------------------
cat('\nDropping uninformative genes sct norm ... ', '\n\n')
drop_genes_sct <- EWCE::drop_uninformative_genes(
  exp = counts_sct,
  input_species = "human",
  output_species = "human",
  level2annot = annotLevels$level2class)

cat('\nGene counts:',
    '\n\nRAW:', dim(raw_counts)[1],
    '\nRAW_NO_MHC:', dim(raw_counts_no_mhc)[1],
    '\nSCT_DROP_GENES:', dim(drop_genes_sct)[1])

rm(raw_counts, raw_counts_no_mhc)

rm(counts_sct)

cat('\nCreating CTD object ... \n\n')
ctd_path <- EWCE::generate_celltype_data(exp = drop_genes_sct,
                                         annotLevels = annotLevels,
                                         groupName = study_id,
                                         savePath = ctd_outdir,
                                         numberOfBins = 10)

## Create enrichment files for MAGMA and LDSR -----------------------------------------

if (study_id == 'herring') {
   sub_dir <- 'herring/'
   magma_end <- paste0('_lvl', level)
} else {
   sub_dir <- 'herring_dwnSmpl/'
   magma_end <- ''
}

#MAGMA input files - 35UP_10DOWN
cat('\nCreating Enrichment files for', study_id, ' ... \n\n')

dir.create(paste0(outdir, sub_dir, 'MAGMA/'),  recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(outdir, sub_dir, 'LDSR/'),  recursive = TRUE, showWarnings = FALSE)

load(paste0(ctd_outdir, 'ctd_', study_id, '.rda'))
CELL_TYPES <- colnames(ctd[[level]]$specificity_quantiles)

MAGMA <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
      inner_join(gene_coord_obj) %>%
      pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
      filter(quantile == 10) %>%
      dplyr::select(cell_type, ensembl) %>%
      with(., split(ensembl, cell_type))

    for(i in names(MAGMA)) {

      cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
          file = paste0(outdir, sub_dir, '/MAGMA/', study_id, magma_end, '.txt'), sep = '', append = TRUE)

    }

#LDSR input files - 100UP_100DOWN
    LDSR <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
        inner_join(gene_coord_obj) %>%
        pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
        filter(quantile == 10) %>%
        mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
        dplyr::select(chr, start, end, ensembl, cell_type) %>%
        group_by(cell_type) %>%
        group_walk(~ write_tsv(.x[,1:4], paste0(outdir, sub_dir, '/LDSR/',
                                                .y$cell_type, '.lvl', level, '.100UP_100DOWN.bed'), col_names = FALSE))
                                                


##Additional analysis to run top 2000 genes rather than top 10% -----------------------

#if (study_id == 'herring') {
#   sub_dir <- 'herring/'
#   magma_end <-	paste0('_lvl', level)

#cat('\nCreating Enrichment files for top 2000 genes ... \n\n')
#load(paste0(ctd_outdir,	'ctd_',	study_id, '.rda'))
#CELL_TYPES <- colnames(ctd[[level]]$specificity_quantiles)

#MAGMA <- as_tibble(ctd[[level]]$specificity, rownames = 'hgnc') %>%
#    inner_join(gene_coord_obj) %>%
#    tidyr::pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'specificity') %>%
#    group_by(cell_type) %>%
#    top_n(n = 2000, wt = specificity) %>%
#    dplyr::select(cell_type, ensembl) %>%
#    with(., split(ensembl, cell_type))

#  for(i in names(MAGMA)) {

#    cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
#        file = paste0(outdir, sub_dir, 'MAGMA/', study_id, '_top2000', magma_end, '.txt')
#        , sep = '', append = TRUE)

#  }

#dir.create(paste0(outdir, sub_dir, 'LDSR_top2000/'),  recursive = TRUE, showWarnings = FALSE)

#  LDSR <- as_tibble(as.matrix(ctd[[level]]$specificity), rownames = 'hgnc') %>%
#    inner_join(gene_coord_obj) %>%
#    pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'specificity') %>%
#    mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
#    group_by(cell_type) %>%
#    top_n(n = 2000, wt = specificity) %>%
#    dplyr::select(chr, start, end, ensembl, cell_type) %>%
#    group_walk(~ write_tsv(.x[,1:4], paste0(outdir, sub_dir, 'LDSR_top2000/',
#                                            .y$cell_type, '.lvl', level,'.100UP_100DOWN.bed'), col_names = FALSE))
#
#  }

##Additional analysis to run protein-coding genes only (top 10%) ----------------------

#if (study_id == 'herring' ) {
#sub_dir <- 'herring/'
#magma_end <- paste0('_lvl', level)
#
#cat('\nCreating Enrichment files for protein-coding genes only ... \n\n')
#load(paste0(ctd_outdir, 'ctd_', study_id, '.rda'))
#CELL_TYPES <- colnames(ctd[[level]]$specificity_quantiles)
#
#MAGMA <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
#      inner_join(protein_gene_coord_obj) %>%
#      pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
#      filter(quantile == 10) %>%
#      dplyr::select(cell_type, ensembl) %>%
#      with(., split(ensembl, cell_type))

 #   for(i in names(MAGMA)) {

#      cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
#          file = paste0(outdir, sub_dir, '/MAGMA/', study_id, '_protein_coding', magma_end, '.txt'), sep = '', append = TRUE)

#    }

#dir.create(paste0(outdir, sub_dir, 'LDSR_protein_coding/'),  recursive = TRUE, showWarnings = FALSE)

#LDSR input files - 100UP_100DOWN
#    LDSR <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
#        inner_join(protein_gene_coord_obj) %>%
#        pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
#        filter(quantile == 10) %>%
#        mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
#        dplyr::select(chr, start, end, ensembl, cell_type) %>%
#        group_by(cell_type) %>%
#        group_walk(~ write_tsv(.x[,1:4], paste0(outdir, sub_dir, 'LDSR_protein_coding/',
#                                                .y$cell_type, '.lvl', level, '.100UP_100DOWN.bed'), col_names = FALSE))
#
#}

file.create(outfile)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

