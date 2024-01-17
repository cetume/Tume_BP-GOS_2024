#--------------------------------------------------------------------------------------
#
#    Plot MAGMA and SLDSR conditional results
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(ggsignif)
library(cowplot)
library(reshape2)
library(ggplot2)
library(dplyr)

## Set variables  ---------------------------------------------------------------------
cat('\nPlotting MAGMA and SLDSR conditional results ... \n')

FIG_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/'
dir.create(paste0(FIG_DIR),  recursive = TRUE, showWarnings = FALSE)
BF_CORR <- 0.05/84
GWAS <- 'SCZ'

## Load and prepare MAGMA data --------------------------------------------------------
cat('\nPreparing MAGMA data ... \n')

magma_dev <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma_conditional/snRNAseq.herring_all_sig_condition_L4_RORB_dev-2.lvl2.magma.35UP_10DOWN.gsa.out'

MAGMA_DF_dev <- read.table(magma_dev, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  dplyr::select(VARIABLE, MAGMA) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category == 'L4_RORB_LRRK1' | Category == 'L4_RORB_MET') %>%
  mutate(Category = ifelse(Category == 'L4_RORB_LRRK1', 'L4_RORB_LRRK1_vs_L4_RORB_dev-2', 'L4_RORB_MET_vs_L4_RORB_dev-2'))

magma_LRRK1 <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma_conditional/snRNAseq.herring_all_sig_condition_L4_RORB_LRRK1.lvl2.magma.35UP_10DOWN.gsa.out'

MAGMA_DF_LRRK1 <- read.table(magma_LRRK1, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  dplyr::select(VARIABLE, MAGMA) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category == 'L4_RORB_dev-2' | Category == 'L4_RORB_MET') %>%
  mutate(Category = ifelse(Category == 'L4_RORB_dev-2', 'L4_RORB_dev-2_vs_L4_RORB_LRRK1', 'L4_RORB_MET_vs_L4_RORB_LRRK1'))

magma_MET <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma_conditional/snRNAseq.herring_all_sig_condition_L4_RORB_MET.lvl2.magma.35UP_10DOWN.gsa.out'

MAGMA_DF_MET <- read.table(magma_MET, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  dplyr::select(VARIABLE, MAGMA) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category == 'L4_RORB_dev-2' | Category == 'L4_RORB_LRRK1') %>%
  mutate(Category = ifelse(Category == 'L4_RORB_dev-2', 'L4_RORB_dev-2_vs_L4_RORB_MET', 'L4_RORB_LRRK1_vs_L4_RORB_MET'))

MAGMA_DF <- merge(MAGMA_DF_dev, MAGMA_DF_LRRK1, all = TRUE) %>%
  merge(MAGMA_DF_MET, all = TRUE) 

## Load and prepare LDSR data ---------------------------------------------------------
cat('\nPreparing LDSR data ... \n')

ldsr <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/LDSR_part_herit/baseline_v1.2/herring_conditional/snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv'

LDSR_FULL_DF <- read_tsv(ldsr) %>%
  mutate(SLDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  mutate(Category = ifelse(Category == 'L4_RORB_dev-2.lvl2.100UP_100DOWN_vs_L4_RORB_LRRK1.lvl2.100UP_100DOWN', 'L4_RORB_dev-2_vs_L4_RORB_LRRK1',
                    ifelse(Category == 'L4_RORB_dev-2.lvl2.100UP_100DOWN_vs_L4_RORB_MET.lvl2.100UP_100DOWN', 'L4_RORB_dev-2_vs_L4_RORB_MET',
                    ifelse(Category == 'L4_RORB_LRRK1.lvl2.100UP_100DOWN_vs_L4_RORB_dev-2.lvl2.100UP_100DOWN', 'L4_RORB_LRRK1_vs_L4_RORB_dev-2',
                    ifelse(Category == 'L4_RORB_LRRK1.lvl2.100UP_100DOWN_vs_L4_RORB_MET.lvl2.100UP_100DOWN', 'L4_RORB_LRRK1_vs_L4_RORB_MET',
                    ifelse(Category == 'L4_RORB_MET.lvl2.100UP_100DOWN_vs_L4_RORB_dev-2.lvl2.100UP_100DOWN', 'L4_RORB_MET_vs_L4_RORB_dev-2', 'L4_RORB_MET_vs_L4_RORB_LRRK1'))))))
  
LDSR_DF <- LDSR_FULL_DF %>%
  dplyr::select(Category, SLDSR)

## Plot data -------------------------------------------------------------------------

PLOT_DF <- left_join(MAGMA_DF, LDSR_DF,
                     by = 'Category') %>% reshape2::melt() 

PLOT_DF$Category <- gsub("_", "-", PLOT_DF$Category)

MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))),
                                              fill = variable, group = rev(variable))) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 12),
        axis.title.y = element_text(colour = "#000000", size = 12),
        axis.text.x = element_text(colour = "#000000", size = 10, vjust = 0.5),
        axis.text.y = element_text(colour = "#000000", size = 7),
        legend.title = element_blank(),
        legend.position = "right") +
  xlab(expression(-log[10](P))) +
  ylab('Conditional analysis') +
  xlim(0, 6) 

PLOT_mean <- PLOT_DF %>% pivot_wider(names_from = variable, values_from = value) 
PLOT_mean$mean <- rowMeans(PLOT_mean[,c('MAGMA', 'SLDSR')]) 

PLOT_mean <- PLOT_mean %>% mutate(COLOUR = ifelse(MAGMA > -log10(BF_CORR) & SLDSR > -log10(BF_CORR), "Both",
                                                  ifelse(MAGMA > -log10(BF_CORR), "MAGMA",
                                                  ifelse(SLDSR > -log10(BF_CORR), "SLDSR", "None")))) 

colour_table <- tibble(
COLOUR = c("Both", "MAGMA", "SLDSR", "None"),
Code = c("#00BA38", "yellow", "#00B0F6", "lightgrey") 
)

PLOT_mean$COLOUR <- factor(PLOT_mean$COLOUR, levels = colour_table$COLOUR) 

MAGMA_LDSR_MEAN_PLOT <- ggplot(data = PLOT_mean, aes(x = mean, y = factor(Category, rev(levels(factor(Category)))), fill = COLOUR)) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  scale_fill_manual(values = colour_table$Code, drop = FALSE, name = "Significant at 
Bonferroni threshold") +
  geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        axis.text.x = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y = element_text(colour = "#000000", size = 12.5),
        title = element_text(colour = "#000000", size = 16),
        legend.text = element_text(colour = "#000000", size = 14),
        strip.text = element_text(size=14, face = 'bold')) +
  xlab(expression(Mean -log[10](P))) +
  ylab('Conditional analysis') +
  xlim(0, 4.5) 

## Save plots -------------------------------------------------------------------------

#jpeg(file = paste0(FIG_DIR, 'SCZ_magma_ldsr_herring_conditional_plot.jpeg'), units = "in", width = 11, height = 7, res = 300) 
#plot(MAGMA_LDSR_PLOT)
#dev.off() 

jpeg(file = paste0(FIG_DIR, 'SCZ_magma_ldsr_herring_mean_conditional_plot.jpeg'), units = "in", width = 11, height = 7, res = 300) 
plot(MAGMA_LDSR_MEAN_PLOT)
dev.off()

