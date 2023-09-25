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
gene_coord <- toString(snakemake@input[['gene_coord']])
protein_gene_coord <- toString(snakemake@input[['protein_gene_coord']])
ctd_object <- toString(snakemake@input[['ctd_object']])
outdir <- toString(snakemake@params[['outdir']])
outfile <- toString(snakemake@output)
study_id <- toString(snakemake@params[['study_id']])
#level <- toString(snakemake@params[['level']])
threads <- snakemake@threads

## Set log file  ----------------------------------------------------------------------
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = c("output", "message"))

## Report inputs  ---------------------------------------------------------------------
cat('\nVariables set to: \n\n')
tibble(Variable = c('gene_coord', 'ctd_object', 'outdir',
                    'outfile', 'study_id', 'threads'),
       Value = c(gene_coord, ctd_object, outdir, outfile, study_id,
                 threads))

dir.create(outdir,  recursive = TRUE, showWarnings = FALSE)

## Load RDS objects  ------------------------------------------------------------------
cat('\nLoad RDS objects ... \n\n')
gene_coord_obj <- readRDS(gene_coord)
protein_gene_coord_object <- readRDS(protein_gene_coord)

## Create enrichment files for MAGMA and LDSR -----------------------------------------

#Experimental function
#create_enrich_test_files(ctd_path = ctd_path,
#                        study_id = study_id,
#                         gene_coordinates = gene_coord_obj,
#                         levels = 2,
#                         outdir = outdir,
#                         sldsr_ext = c(100, 100),
#                         magma_ext = c(35, 10))

#if (study_id == 'herring_dwnSmpl_lvl2') {
#   levels <- 2
#} else if (study_id == 'herring_dwnSmpl_lvl1') {
#   levels <- 1
#} else {
#   levels <- c(1,2)
#}
#
#for(level in levels) {
#if (study_id == 'herring') {
#   sub_dir <- 'herring/'
#   magma_end <- paste0('_lvl', level)
#} else {
#   sub_dir <- 'herring_dwnSmpl/'
#   magma_end <- ''
#}
#
#MAGMA input files - 35UP_10DOWN
#cat('\nCreating Enrichment files for', study_id, ' ... \n\n')
#
#dir.create(paste0(outdir, sub_dir, 'MAGMA/'),  recursive = TRUE, showWarnings = FALSE)
#dir.create(paste0(outdir, sub_dir, 'LDSR/'),  recursive = TRUE, showWarnings = FALSE)
#
#load(ctd_object)
#CELL_TYPES <- colnames(ctd[[level]]$specificity_quantiles)
#
#MAGMA <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
#      inner_join(gene_coord_obj) %>%
#      pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
#      filter(quantile == 10) %>%
#      dplyr::select(cell_type, ensembl) %>%
#      with(., split(ensembl, cell_type))
#
#    for(i in names(MAGMA)) {
#
#      cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
#          file = paste0(outdir, sub_dir, '/MAGMA/', study_id, magma_end, '.txt'), sep = '', append = TRUE)
#
#    }
#
#LDSR input files - 100UP_100DOWN
#    LDSR <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
#        inner_join(gene_coord_obj) %>%
#        pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
#        filter(quantile == 10) %>%
#        mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
#        dplyr::select(chr, start, end, ensembl, cell_type) %>%
#        group_by(cell_type) %>%
#        group_walk(~ write_tsv(.x[,1:4], paste0(outdir, sub_dir, '/LDSR/',
#                                                .y$cell_type, '.lvl', level, '.100UP_100DOWN.bed'), col_names = FALSE))
#
# }

##Additional analysis to run top 2000 genes rather than top 10% -----------------------

#for(level in c(1,2)) {
#if (study_id == 'herring' ) {
#
#sub_dir <- 'herring/'
#magma_end <- paste0('_lvl', level)
#
#cat('\nCreating Enrichment files for top 2000 genes ... \n\n')
#load(ctd_object)
#CELL_TYPES <- colnames(ctd[[level]]$specificity_quantiles)
#
#MAGMA <- as_tibble(ctd[[level]]$specificity, rownames = 'hgnc') %>%
#    inner_join(gene_coord_obj) %>%
#    tidyr::pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'specificity') %>%
#    group_by(cell_type) %>%
#    top_n(n = 2000, wt = specificity) %>%
#    dplyr::select(cell_type, ensembl) %>%
#    with(., split(ensembl, cell_type))
#
#  for(i in names(MAGMA)) {
#
#    cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
#        file = paste0(outdir, sub_dir, 'MAGMA/', study_id, '_top2000', magma_end, '.txt')
#        , sep = '', append = TRUE)
#
#  }
#
#dir.create(paste0(outdir, sub_dir, 'LDSR_top2000_genes/'),  recursive = TRUE, showWarnings = FALSE)
#
#  LDSR <- as_tibble(as.matrix(ctd[[level]]$specificity), rownames = 'hgnc') %>%
#    inner_join(gene_coord_obj) %>%
#    pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'specificity') %>%
#    mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
#    group_by(cell_type) %>%
#    top_n(n = 2000, wt = specificity) %>%
#    dplyr::select(chr, start, end, ensembl, cell_type) %>%
#    group_walk(~ write_tsv(.x[,1:4], paste0(outdir, sub_dir, 'LDSR_top2000_genes/',
#                                            .y$cell_type, '.lvl', level,'.100UP_100DOWN.bed'), col_names = FALSE))
#
#  }
#}
#
##Additional analysis to run protein-coding genes only (top 10%) ----------------------

for(level in c(1,2)) {
if (study_id == 'herring' ) {

sub_dir <- 'herring/'
magma_end <- paste0('_lvl', level)

cat('\nCreating Enrichment files for protein-coding genes only ... \n\n')
load(ctd_object)
CELL_TYPES <- colnames(ctd[[level]]$specificity_quantiles)

MAGMA <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
      inner_join(protein_gene_coord_obj) %>%
      pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
      filter(quantile == 10) %>%
      dplyr::select(cell_type, ensembl) %>%
      with(., split(ensembl, cell_type))

    for(i in names(MAGMA)) {

      cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n",
          file = paste0(outdir, sub_dir, '/MAGMA/', study_id, '_protein_coding', magma_end, '.txt'), sep = '', append = TRUE)

    }

dir.create(paste0(outdir, sub_dir, 'LDSR_protein_coding/'),  recursive = TRUE, showWarnings = FALSE)

LDSR input files - 100UP_100DOWN
    LDSR <- as_tibble(as.matrix(ctd[[level]]$specificity_quantiles), rownames = 'hgnc') %>%
        inner_join(protein_gene_coord_obj) %>%
        pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
        filter(quantile == 10) %>%
        mutate(start = ifelse(start - 100000 < 0, 0, start - 100000), end = end + 100000) %>%
        dplyr::select(chr, start, end, ensembl, cell_type) %>%
        group_by(cell_type) %>%
        group_walk(~ write_tsv(.x[,1:4], paste0(outdir, sub_dir, 'LDSR_protein_coding/',
                                                .y$cell_type, '.lvl', level, '.100UP_100DOWN.bed'), col_names = FALSE))

 }








file.create(outfile)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
