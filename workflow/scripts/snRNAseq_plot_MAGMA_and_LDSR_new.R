#--------------------------------------------------------------------------------------
#
#    Plot MAGMA and SLDSR results
#
#--------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
print(R.version)
.libPaths()

##  Load Packages  --------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(janitor)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(reshape2)

## Set variables  ---------------------------------------------------------------------
cat('\nPlotting MAGMA and SLDSR results ... \n')

magma <- toString(snakemake@input[['magma']])
ldsr <- toString(snakemake@input[['ldsr']])
GWAS <- toString(snakemake@params[['GWAS']])
study_id <- toString(snakemake@params[['study_id']])
level <- toString(snakemake@params[['level']])
fig_dir <- toString(snakemake@params[['fig_dir']])

## Report inputs  ---------------------------------------------------------------------
#cat('\nVariables set to: \n\n')
#tibble(Variable = c('magma', 'ldsr', 'GWAS', 'level', 'fig_dir', 'study_id'),
#       Value = c(magma, ldsr, GWAS, level, fig_dir, study_id))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

##Preparing MAGMA data ----------------------------------------------------------------
cat('\nPreparing MAGMA data ... \n')

MAGMA_DF <- read.table(magma, header = FALSE) %>%
        janitor::row_to_names(row_number = 1) %>%
        mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
        mutate(MAGMA = -log10(as.numeric(P))) %>%
        dplyr::select(VARIABLE, MAGMA) %>%
        dplyr::rename(Category = VARIABLE)

##Preparing LDSR data -----------------------------------------------------------------
cat('\nPreparing LDSR data ... \n')

LDSR_FULL_DF <- read_tsv(ldsr) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          filter(str_detect(Category, paste0('lvl', level, '.100UP_100DOWN'))) %>%
          separate(Category, into=c('Category', 'Window'), sep = '\\.', extra = "merge")

LDSR_DF <- LDSR_FULL_DF %>%
          dplyr::select(Category, LDSR)

##Plot MAGMA and LDSR barplot for lvl1 ------------------------------------------------
cat('\nCreate MAGMA and LDSR joint plots ... \n')

if (level == 1) {

BF_CORR <- 0.05/19
WIDTH <- 6

PLOT_DF <- left_join(MAGMA_DF, LDSR_DF,
                         by = 'Category') %>% reshape2::melt()


MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))),
                                                  fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(GWAS) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 12),
            axis.title.y = element_text(colour = "#000000", size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 7),
            legend.title = element_blank(),
            legend.position = "top") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8)

 #assign(paste0(GWAS, '_magma_ldsr_', study_id, '_lvl', level, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv)

##Plot average of MAGMA and LDSR for lvl1 ---------------------------------------------
cat('\nCreate MAGMA and LDSR average plots ... \n')

PLOT_mean <- PLOT_DF %>% pivot_wider(names_from = variable, values_from = value)
PLOT_mean$mean <- rowMeans(PLOT_mean[,c('MAGMA', 'LDSR')])

PLOT_mean <- PLOT_mean %>% mutate(COLOUR = ifelse(MAGMA > -log10(BF_CORR) & LDSR > -log10(BF_CORR), "Both",
                                                    ifelse(MAGMA > -log10(BF_CORR) & LDSR > -log10(0.05), "MAGMA  (LDSR P < 0.05)",
                                                           ifelse(LDSR > -log10(BF_CORR) & MAGMA > -log10(0.05), "LDSR  (MAGMA P < 0.05)", "None"))))

colour_table <- tibble(
  COLOUR = c("Both", "MAGMA  (LDSR P < 0.05)", "None", "LDSR  (MAGMA P < 0.05)"),
  Code = c("#00BA38", "yellow", "lightgrey", "#00B0F6")
  )

PLOT_mean$COLOUR <- factor(PLOT_mean$COLOUR, levels = colour_table$COLOUR)

    MAGMA_LDSR_MEAN_PLOT <- ggplot(data = PLOT_mean, aes(x = mean, y = factor(Category, rev(levels(factor(Category)))), fill = COLOUR)) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
      scale_fill_manual(values = colour_table$Code, drop = FALSE, name = "Significant") +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(GWAS) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 12),
            axis.title.y = element_text(colour = "#000000", size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 7),
            legend.position = "top") +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8)

#assign(paste0(GWAS, '_magma_ldsr_mean_', study_id, 'lvl', level, '_plot'), MAGMA_LDSR_MEAN_PLOT, envir = .GlobalEnv)

} else {

BF_CORR <- 0.05/84
WIDTH <- 11

PLOT_DF <- left_join(MAGMA_DF, LDSR_DF,
                         by = 'Category') %>% reshape2::melt() %>%
                         mutate(cell_type = case_when(Category %in% c("CCK_RELN", "CCK_SORCS1", "CCK_SYT6", "CGE_dev", "ID2_CSMD1", "ID2_dev", "LAMP5_CCK", "LAMP5_NDNF", "LAMP5_NOS1", "MGE_dev-1", "MGE_dev-2", "PV_dev", "PV_SCUBE3", "PV_SCUBE3_dev", "PV_SST", "PV_SULF1", "PV_SULF1_dev", "PV_WFDC2", "SST_ADGRG6", "SST_ADGRG6_dev", "SST_B3GAT2", "SST_BRINP3", "SST_CALB1", "SST_CALB1_dev", "SST_NPY", "SST_STK32A", "SST_TH", "VIP_ABI3BP", "VIP_ADAMTSL1", "VIP_CHRM2", "VIP_CRH", "VIP_dev", "VIP_DPP6", "VIP_HS3ST3A1", "VIP_KIRREL3", "VIP_PCDH20") ~ "IN",
                               Category %in% c("L2_CUX2_LAMP5", "L3_CUX2_PRSS12", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME", "L5-6_THEMIS_CNR1", "L5-6_THEMIS_NTNG2", "L5-6_TLE4_HTR2C", "L5-6_TLE4_SCUBE1", "L5-6_TLE4_SORCS1") ~ "PN",
                               Category %in% c("L2-3_CUX2_dev-1", "L2-3_CUX2_dev-2", "L2-3_CUX2_dev-3", "L2-3_CUX2_dev-4", "L2-3_CUX2_dev-5", "L2-3_CUX2_dev-6", "L2-3_CUX2_dev-fetal", "L2_CUX2_LAMP5_dev", "L4_RORB_dev-1", "L4_RORB_dev-2", "L4_RORB_dev-fetal", "L5-6_THEMIS_dev-1", "L5-6_THEMIS_dev-2", "L5-6_TLE4_dev", "PN_dev") ~ "PN_dev",
                               Category %in% c("Astro_dev-1", "Astro_dev-2", "Astro_dev-3", "Astro_dev-4", "Astro_dev-5", "Astro_GFAP", "Astro_SLC1A2", "Astro_SLC1A2_dev", "Micro", "Micro_out", "Oligo-1", "Oligo-2", "Oligo-3", "Oligo-4", "Oligo-5", "Oligo_mat", "OPC", "OPC_dev", "OPC_MBP", "Vas_CLDN5", "Vas_PDGFRB", "Vas_TBX18", "BKGR_NRGN") ~ "Non-Neu/Poor-Quality")) %>%
                         mutate(cell_type=factor(cell_type, levels = c("PN", "PN_dev", "IN", "Non-Neu/Poor-Quality")))


    MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))),
                                                  fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(GWAS) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 12),
            axis.title.y = element_text(colour = "#000000", size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 7),
            legend.title = element_blank(),
            legend.position = "top") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8) +
      facet_wrap(~cell_type, scales = "free")

#assign(paste0(GWAS, '_magma_ldsr_', study_id, '_lvl', level, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv)

##Plot average of MAGMA and LDSR ------------------------------------------------------
cat('\nCreate MAGMA and LDSR average plots ... \n')

PLOT_mean <- PLOT_DF %>% pivot_wider(names_from = variable, values_from = value)
PLOT_mean$mean <- rowMeans(PLOT_mean[,c('MAGMA', 'LDSR')])

PLOT_mean <- PLOT_mean %>% mutate(COLOUR = ifelse(MAGMA > -log10(BF_CORR) & LDSR > -log10(BF_CORR), "Both",
                                                    ifelse(MAGMA > -log10(BF_CORR) & LDSR > -log10(0.05), "MAGMA  (LDSR P < 0.05)",
                                                           ifelse(LDSR > -log10(BF_CORR) & MAGMA > -log10(0.05), "LDSR  (MAGMA P < 0.05)", "None"))))

colour_table <- tibble(
  COLOUR = c("Both", "MAGMA  (LDSR P < 0.05)", "None", "LDSR  (MAGMA P < 0.05)"),
  Code = c("#00BA38", "yellow", "lightgrey", "#00B0F6")
  )

PLOT_mean$COLOUR <- factor(PLOT_mean$COLOUR, levels = colour_table$COLOUR)

    MAGMA_LDSR_MEAN_PLOT <- ggplot(data = PLOT_mean, aes(x = mean, y = factor(Category, rev(levels(factor(Category)))), fill = COLOUR)) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
      scale_fill_manual(values = colour_table$Code, drop = FALSE, name = "Significant") +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(GWAS) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 12),
            axis.title.y = element_text(colour = "#000000", size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 7),
            legend.position = "top") +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8) +
      facet_wrap(~cell_type, scales = "free")
      
}

#assign(paste0(GWAS, '_magma_ldsr_mean_', study_id, 'lvl', level, '_plot'), MAGMA_LDSR_MEAN_PLOT, envir = .GlobalEnv)

##Produce final plots -----------------------------------------------------------------

end <- paste0('_lvl', level)

jpeg(file = paste0(fig_dir, GWAS, '_magma_ldsr_mean_', study_id, end, '_plot.jpeg'), units = "in", width = WIDTH, height = 11, res = 300)
plot(MAGMA_LDSR_MEAN_PLOT)
dev.off()

jpeg(file = paste0(fig_dir, GWAS, '_magma_ldsr_', study_id, end, '_plot.jpeg'), units = "in", width = WIDTH, height = 11, res = 300)
plot(MAGMA_LDSR_PLOT)
dev.off()
