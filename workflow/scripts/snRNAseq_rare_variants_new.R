## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------
library(dplyr)
library(readr)
library(tibble)
library(EWCE)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(reshape2)

## Set variables  ---------------------------------------------------------------------
cat('\nRare variant analysis ... \n')

ctd_object <- toString(snakemake@input[['ctd_object']])
schema <- toString(snakemake@input[['schema']])
schema_genes <- toString(snakemake@input[['schema_genes']])
study_id <- toString(snakemake@params[['study_id']])
rare_dir <- toString(snakemake@params[['rare_dir']])
boxplots_dir <- toString(snakemake@params[['boxplots_dir']])
fig_dir <- toString(snakemake@params[['fig_dir']])
ctd_dir <- toString(snakemake@params[['ctd_dir']])
outfile <- toString(snakemake@output)

dir.create(rare_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(boxplots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

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
    jpeg(paste0(boxplots_dir, CELL_TYPE, '_lvl', level, '_boxplot.jpg'), width = 960, height = 960,
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
    mutate(FDR = p.adjust(P, 'BH', length(P)))

write_tsv(wilcoxon_df, paste0(rare_dir, 'wilcoxon_df_', study_id, end, '.txt'))

##Write specificity scores ------------------------------------------------------------
cat('\nWriting specificity scores ... \n')
write_tsv(specificity_DF, paste0(ctd_dir, 'specificity_df_', study_id, end, '.txt'))

##Preparing rare variant data for plots -----------------------------------------------
cat('\nPreparing rare variant data for plots ... \n')

wilcoxon_df$P <- as.numeric(as.character(wilcoxon_df$P))

if (level == 1) {

BF_CORR <- 0.05/19
WIDTH <- 6

WILCOXON_PLOT <- ggplot(data = wilcoxon_df, aes(x = -log10(P), y = factor(Category, rev(levels(factor(Category)))))) +
      geom_bar(stat = "identity", color = 'black', fill = 'darkseagreen') +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      #ggtitle(TITLE) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 12),
            axis.title.y = element_text(colour = "#000000", size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 7),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 5)


} else {

BF_CORR <- 0.05/84
WIDTH <- 11

PLOT_DF <- wilcoxon_df %>%
           mutate(cell_type = case_when(Category %in% c("CCK_RELN", "CCK_SORCS1", "CCK_SYT6", "CGE_dev", "ID2_CSMD1", "ID2_dev", "LAMP5_CCK", "LAMP5_NDNF", "LAMP5_NOS1", "MGE_dev-1", "MGE_dev-2", "PV_dev", "PV_SCUBE3", "PV_SCUBE3_dev", "PV_SST", "PV_SULF1", "PV_SULF1_dev", "PV_WFDC2", "SST_ADGRG6", "SST_ADGRG6_dev", "SST_B3GAT2", "SST_BRINP3", "SST_CALB1", "SST_CALB1_dev", "SST_NPY", "SST_STK32A", "SST_TH", "VIP_ABI3BP", "VIP_ADAMTSL1", "VIP_CHRM2", "VIP_CRH", "VIP_dev", "VIP_DPP6", "VIP_HS3ST3A1", "VIP_KIRREL3", "VIP_PCDH20") ~ "IN",
                Category %in% c("L2_CUX2_LAMP5", "L3_CUX2_PRSS12", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME", "L5-6_THEMIS_CNR1", "L5-6_THEMIS_NTNG2", "L5-6_TLE4_HTR2C", "L5-6_TLE4_SCUBE1", "L5-6_TLE4_SORCS1") ~ "PN",
                Category %in% c("L2-3_CUX2_dev-1", "L2-3_CUX2_dev-2", "L2-3_CUX2_dev-3", "L2-3_CUX2_dev-4", "L2-3_CUX2_dev-5", "L2-3_CUX2_dev-6", "L2-3_CUX2_dev-fetal", "L2_CUX2_LAMP5_dev", "L4_RORB_dev-1", "L4_RORB_dev-2", "L4_RORB_dev-fetal", "L5-6_THEMIS_dev-1", "L5-6_THEMIS_dev-2", "L5-6_TLE4_dev", "PN_dev") ~ "PN_dev",
                Category %in% c("Astro_dev-1", "Astro_dev-2", "Astro_dev-3", "Astro_dev-4", "Astro_dev-5", "Astro_GFAP", "Astro_SLC1A2", "Astro_SLC1A2_dev", "Micro", "Micro_out", "Oligo-1", "Oligo-2", "Oligo-3", "Oligo-4", "Oligo-5", "Oligo_mat", "OPC", "OPC_dev", "OPC_MBP", "Vas_CLDN5", "Vas_PDGFRB", "Vas_TBX18", "BKGR_NRGN") ~ "Non-Neu/Poor-Quality")) %>%
           mutate(cell_type=factor(cell_type, levels = c("PN", "PN_dev", "IN", "Non-Neu/Poor-Quality")))

 WILCOXON_PLOT <- ggplot(data = PLOT_DF, aes(x = -log10(P), y = factor(Category, rev(levels(factor(Category)))))) +
      geom_bar(stat = "identity", color = 'black', fill = 'darkseagreen') +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      #ggtitle(TITLE) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 12),
            axis.title.y = element_text(colour = "#000000", size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 7),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 5) +
      facet_wrap(~cell_type, scales = "free")
      
}
      
##Produce final plots --------------------------------------------------------

cat('\nSave plots ... \n')

 jpeg(file = paste0(fig_dir, 'wilcoxon_', study_id, end, '_plot.jpeg'), units = "in", width = WIDTH, height = 11, res = 300)
 plot(WILCOXON_PLOT)
 dev.off()

}

file.create(outfile)
