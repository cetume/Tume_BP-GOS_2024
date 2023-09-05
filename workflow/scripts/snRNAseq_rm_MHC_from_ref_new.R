#--------------------------------------------------------------------------------------
#
#    Remove MHC genes from reference file
#
#--------------------------------------------------------------------------------------

##  Set Variables  --------------------------------------------------------------------
cat('\nRemove MHC genes from gene ref file ... \n')
ref_in <- toString(snakemake@input)
ref_noMHC <- toString(snakemake@params[['ref_noMHC']])
gene_coord_rds <- toString(snakemake@output[['gene_coord']])
mhc_genes_rds <- toString(snakemake@output[['mhc_genes']])

## Initialise R library  --------------------------------------------------------------
R.Version()
##  Load Packages  --------------------------------------------------------------------
library(dplyr)
library(readr)
library(biomaRt)

# Get MHC genes -----------------------------------------------------------------------
cat('Retrieving MHC genes ... \n')
mart_hg19 <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://grch37.ensembl.org')
mart_hg19 <- useDataset('hsapiens_gene_ensembl', mart_hg19)
mhc_genes <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
                   filters = c("chromosome_name","start","end"),
                   values = list(chromosome = "6", start = "28477797", end = "33448354"), #hg19 coords (https://doi.org/10.1038/s41431-019-0559-2)
                   mart = mart_hg19)
mhc_genes_uniq <- stringi::stri_remove_empty(unique(mhc_genes$external_gene_name), na_empty = FALSE)
cat('MHC genes:', length(mhc_genes_uniq), '\n')

## Generate gene coordinate reference file ---------------------------------------------
annot_lookup_hg19 <- as_tibble(
  getBM(
    mart = mart_hg19,
    attributes = c(
      'chromosome_name', 
      'start_position', 
      'end_position',       
      'gene_biotype', 
      'hgnc_symbol',
      'ensembl_gene_id',
      'strand',
      'external_gene_name'),
    uniqueRows = TRUE))
annot_lookup_hg19$strand <- ifelse(annot_lookup_hg19$strand == "-1", "-", "+")

# Identify and remove duplicates from reference file
annot_lookup_hg19 <- annot_lookup_hg19 %>% filter(!grepl("H", chromosome_name) & !grepl("G", chromosome_name)) #Remove patches
annot_lookup_hg19 <- annot_lookup_hg19[!duplicated(annot_lookup_hg19$ensembl_gene_id), ] #Remove duplicate ensembl_gene_ids - have identical gene co-ordinates

n_occur <- data.frame(table(annot_lookup_hg19$external_gene_name))
duplicates <- n_occur[n_occur$Freq > 1,]
colnames(duplicates) <- c("external_gene_name", "Freq")

annot_lookup_hg19_dup <- annot_lookup_hg19[annot_lookup_hg19$external_gene_name %in% duplicates$external_gene_name, ] #all duplicates
annot_lookup_hg19_dup_with_hgnc <- annot_lookup_hg19_dup %>% filter(hgnc_symbol != "") #duplicates with hgnc_symbols (mostly protein-coding)
duplicates_to_remove <- annot_lookup_hg19_dup[!annot_lookup_hg19_dup$ensembl_gene_id %in% annot_lookup_hg19_dup_with_hgnc$ensembl_gene_id, ] #remove duplicates with no hgnc_symbol
annot_lookup_hg19_filtered <- annot_lookup_hg19[!annot_lookup_hg19$ensembl_gene_id %in% duplicates_to_remove$ensembl_gene_id, ]
annot_lookup_hg19_filtered <- annot_lookup_hg19_filtered[!duplicated(annot_lookup_hg19_filtered$external_gene_name), ] #remove remaining duplicates (7 genes) 

# Remove MHC gene from gene coordinate reference file (and write to file)  ------------
cat('Removing MHC genes from gene coord ref (and writing) ...\n')
gene_coordinates <- annot_lookup_hg19_filtered %>% 
  dplyr::rename(chr = "chromosome_name", ensembl = "ensembl_gene_id", 
                start = 'start_position', end = 'end_position', hgnc = 'external_gene_name') %>% 
  dplyr::select(ensembl, chr, start, end, strand, hgnc) %>%
  filter(!hgnc %in% mhc_genes_uniq) %>% 
  write_tsv(ref_noMHC, col_names = FALSE) %>%
  mutate(chr = paste0("chr", chr)) %>%
  dplyr::select(chr, start, end, ensembl, hgnc)  

saveRDS(mhc_genes_uniq, mhc_genes_rds)
saveRDS(gene_coordinates, gene_coord_rds)

cat('Done.\n')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
