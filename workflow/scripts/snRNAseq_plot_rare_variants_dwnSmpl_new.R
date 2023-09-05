## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(reshape2)

## Set variables  ---------------------------------------------------------------------
cat('\nPlotting rare variant analysis ... \n')

rare_variant_data <- toString(snakemake@input)
fig_dir <- toString(snakemake@params[['fig_dir']])
level <- toString(snakemake@params[['level']])

##Preparing rare variant data ---------------------------------------------------------
cat('\nPreparing rare variant data ... \n')

wilcoxon_df <- read_tsv(rare_variant_data)

if (level == 1) {

BF_CORR <- 0.05/19

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
      
##Produce final plots for lvl1 --------------------------------------------------------

 jpeg(file = paste0(fig_dir, 'wilcoxon_herring_dwnSmpl_lvl', level, '_plot.jpeg'), units = "in", width = 6, height = 11, res = 300)
 plot(WILCOXON_PLOT)
 dev.off()

} else {

BF_CORR <- 0.05/84

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
      

##Produce final plots for lvl2 --------------------------------------------------------

 jpeg(file = paste0(fig_dir, 'wilcoxon_herring_dwnSmpl_lvl', level, '_plot.jpeg'), units = "in", width = 11, height = 11, res = 300)
 plot(WILCOXON_PLOT)
 dev.off()
 
 }
