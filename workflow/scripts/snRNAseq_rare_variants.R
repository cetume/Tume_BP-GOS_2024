## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------
library(dplyr)
library(readr)
library(tibble)
library(EWCE)
library(ggplot2)

## Set variables  ---------------------------------------------------------------------
cat('\nRare variant analysis ... \n')

ctd_object <- toString(snakemake@input[['ctd_object']])
schema <- toString(snakemake@input[['schema']])
schema_genes <- toString(snakemake@input[['schema_genes']])
study_id <- toString(snakemake@params[['study_id']])
rare_dir <- toString(snakemake@params[['rare_dir']])
fig_dir <- toString(snakemake@params[['fig_dir']])
ctd_dir <- toString(snakemake@params[['ctd_dir']])
outfile <- toString(snakemake@output)

dir.create(rare_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

##Load data ---------------------------------------------------------------------------
cat('\nLoad schema and ctd data ... \n')

schema <- read_csv(schema) %>%
  head(32) %>%
  dplyr::select('Gene', 'Q meta') %>%
  rename(Q = 'Q meta', gene = 'Gene')

#schema genes - csv genes have ensemble encoding so need conversion file
schema_genes <- read_tsv(schema_genes) %>%
  left_join(schema) %>%
  arrange(Q) %>%
  rename(ensemble_gene = gene,
         gene = symbol)

#ctd specificity scores
load(ctd_object)

if(study_id == 'herring'){

levels <- c(1, 2)

} else if(study_id == 'herring_dwnSmpl_lvl2') {

levels <- 2

} else {

levels <- 1

}

for(level in levels){

ctd_specificity <- as.data.frame.matrix(ctd[[level]]$specificity)

##Run Wilcoxon Rank test --------------------------------------------------------------
cat('\nRunning Wilcoxon rank test ... \n')

STUDY_DF <- schema_genes
STUDY_GENES <- STUDY_DF$gene
wilcoxon_df <- data.frame()
i <- 1

#Load regional specificity scores

specificity_DF <- ctd_specificity %>%
      rownames_to_column(var = 'gene')

for (CELL_TYPE in colnames(specificity_DF)) {

 if (CELL_TYPE == 'gene') next

 cat('\n\nRunning wilcoxon test for:', CELL_TYPE)

specificity_cell <- data.frame(gene = specificity_DF$gene,
    cell_scores = specificity_DF[[CELL_TYPE]]) %>%
    mutate(study_status = ifelse(gene %in% STUDY_GENES, 'in_study', 'not_in_study'))

in_study_scores <- specificity_cell %>%
    filter(grepl('^in_study', study_status)) %>%
    pull(cell_scores)

not_in_study_scores <- specificity_cell %>%
    filter(grepl('^not_in_study', study_status)) %>%
    pull(cell_scores)

wilcox_result <- wilcox.test(in_study_scores, not_in_study_scores, alternative = 'greater',
                             paired = FALSE, data = specificity_cell)


wilcoxon_df <- rbind(wilcoxon_df, as.data.frame(t(c(CELL_TYPE, wilcox_result$p.value))))


wilcoxon_boxplot <- ggplot(specificity_cell, aes(x = study_status, y = cell_scores)) +
    geom_boxplot() + labs(title = CELL_TYPE) + theme_bw()
    jpeg(paste0(fig_dir, CELL_TYPE, '_lvl', level, '_boxplot.jpg'), width = 960, height = 960,
           units = "px", pointsize = 12)
    print(wilcoxon_boxplot)
    dev.off()
    cat(paste0('\nNumber of in sample:'), sum(specificity_cell$study_status == 'in_study'))

    }

colnames(wilcoxon_df) <- c('Category', 'P')
print(wilcoxon_df)

##Write Wilcoxon Rank test results ----------------------------------------------------
cat('\nWriting Wilcoxon rank test results ... \n')
if(study_id == 'herring') {

end <- paste0('_lvl', level)

} else {


end <- ''

}

wilcoxon_df <- wilcoxon_df %>%
    arrange(P) %>%
    mutate(BF = p.adjust(P, 'bonferroni', length(P))) %>%
    mutate(FDR = p.adjust(P, 'BH', length(P))) %>%
    write_tsv(paste0(rare_dir, 'wilcoxon_df_', study_id, end, '.txt'))

##Write specificity scores ------------------------------------------------------------
cat('\nWriting specificity scores ... \n')
write_tsv(specificity_DF, paste0(ctd_dir, 'specificity_df_', study_id, end, '.txt'))

}

file.create(outfile)
