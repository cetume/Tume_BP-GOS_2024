#--------------------------------------------------------------------------------------
#
#    Wilcoxon rank sum tests for fine-mapped prioritised genes 
#
#--------------------------------------------------------------------------------------
#Runs Wilcoxon rank sum tests for the specificity scores of the fine-mapped prioritised genes vs all other genes for each cell population
#Plots the -log10(P-values) for each cell population

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
library(ggforce)

## Set variables  ---------------------------------------------------------------------
cat('\nFinemapped variant analysis ... \n')

ctd_object <- "~/Desktop/Herring_snRNAseq_2023_pipeline/results/ctd_objects/ctd_herring.rda"
prioritised_genes <- "~/Desktop/Herring_snRNAseq_2023_pipeline/resources/prioritised_genes/GWAS_prioritised_genes_finemapped.csv"
study_id <- "herring"
prioritised_dir <- "~/Desktop/Herring_snRNAseq_2023_pipeline/resources/prioritised_genes/"
boxplots_dir <- "~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/boxplots/"
fig_dir <- "~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/"
ctd_dir <- "~/Desktop/Herring_snRNAseq_2023_pipeline/results/ctd_objects/"
output_dir <- "~/Desktop/Herring_snRNAseq_2023_pipeline/results/prioritised_genes/"
LEVEL <- 2

dir.create(prioritised_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(boxplots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

##Load data ---------------------------------------------------------------------------
cat('\nLoad prioritised gene and ctd data ... \n')

prioritised_genes <- read_csv(prioritised_genes)

#ctd specificity scores
load(ctd_object)

ctd_specificity <- as.data.frame.matrix(ctd[[LEVEL]]$specificity)

##Run Wilcoxon Rank test --------------------------------------------------------------
cat('\nRunning Wilcoxon rank test ... \n')

STUDY_DF <- prioritised_genes
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
    filter(cell_scores > 0) %>%
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

write_tsv(wilcoxon_df, paste0(output_dir, 'wilcoxon_df_', study_id, end, '.txt'))

##Preparing rare variant data for plots -----------------------------------------------
cat('\nPreparing prioritised genes data for plots ... \n')

wilcoxon_df$P <- as.numeric(as.character(wilcoxon_df$P))

BF_CORR <- 0.05/84

PLOT_DF <- wilcoxon_df %>%
  mutate(cell_type = case_when(Category %in% c("CCK_RELN", "CCK_SORCS1", "CCK_SYT6", "ID2_CSMD1", "LAMP5_CCK", "LAMP5_NDNF", "LAMP5_NOS1", "PV_SCUBE3", "PV_SST", "PV_SULF1", "PV_WFDC2", "SST_ADGRG6", "SST_B3GAT2", "SST_BRINP3", "SST_CALB1", "SST_NPY", "SST_STK32A", "SST_TH", "VIP_ABI3BP", "VIP_ADAMTSL1", "VIP_CHRM2", "VIP_CRH", "VIP_DPP6", "VIP_HS3ST3A1", "VIP_KIRREL3", "VIP_PCDH20") ~ "Mature Inhibitory Neurons",
                               Category %in% c("CGE_dev", "ID2_dev", "MGE_dev-1", "MGE_dev-2", "PV_dev", "PV_SCUBE3_dev", "PV_SULF1_dev", "SST_ADGRG6_dev", "SST_CALB1_dev", "VIP_dev") ~ "Developing Inhibitory Neurons",
                               Category %in% c("BKGR_NRGN", "L2_CUX2_LAMP5", "L3_CUX2_PRSS12", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME", "L5-6_THEMIS_CNR1", "L5-6_THEMIS_NTNG2", "L5-6_TLE4_HTR2C", "L5-6_TLE4_SCUBE1", "L5-6_TLE4_SORCS1") ~ "Mature Excitatory Neurons",
                               Category %in% c("L2-3_CUX2_dev-1", "L2-3_CUX2_dev-2", "L2-3_CUX2_dev-3", "L2-3_CUX2_dev-4", "L2-3_CUX2_dev-5", "L2-3_CUX2_dev-6", "L2-3_CUX2_dev-fetal", "L2_CUX2_LAMP5_dev", "L4_RORB_dev-1", "L4_RORB_dev-2", "L4_RORB_dev-fetal", "L5-6_THEMIS_dev-1", "L5-6_THEMIS_dev-2", "L5-6_TLE4_dev", "PN_dev") ~ "Developing Excitatory Neurons",
                               Category %in% c("Astro_dev-1", "Astro_dev-2", "Astro_dev-3", "Astro_dev-4", "Astro_dev-5", "Astro_GFAP", "Astro_SLC1A2", "Astro_SLC1A2_dev", "Micro", "Micro_out", "Oligo-1", "Oligo-2", "Oligo-3", "Oligo-4", "Oligo-5", "Oligo_mat", "OPC", "OPC_dev", "OPC_MBP", "Vas_CLDN5", "Vas_PDGFRB", "Vas_TBX18") ~ "Non-Neuronal")) %>%
  mutate(cell_type=factor(cell_type, levels = c("Mature Excitatory Neurons", "Developing Excitatory Neurons", "Mature Inhibitory Neurons", "Developing Inhibitory Neurons", "Non-Neuronal")))

PLOT_DF$Category <- gsub("_", "-", PLOT_DF$Category)

WILCOXON_PLOT <- ggplot(data = PLOT_DF, aes(x = -log10(P), y = factor(Category, rev(levels(factor(Category)))))) +
      geom_bar(stat = "identity", color = 'black', fill = 'tan2') +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      #ggtitle(TITLE) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 16),
            axis.title.y = element_text(colour = "#000000", size = 16),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            strip.text = element_text(size=14, face = 'bold'),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 4) +
  facet_col(facets = vars(cell_type), scales = "free", space = "free", strip.position = "top")
      
##Produce final plots --------------------------------------------------------

cat('\nSave plots ... \n')

 jpeg(file = paste0(fig_dir, 'wilcoxon_prioritised_genes', study_id, end, '_plot.jpeg'), units = "in", width = 7.2, height = 22, res = 300)
 plot(WILCOXON_PLOT)
 dev.off()
