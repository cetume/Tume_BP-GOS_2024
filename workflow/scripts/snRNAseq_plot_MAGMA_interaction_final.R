#--------------------------------------------------------------------------------------
#
#    Plot MAGMA interaction results (betas)
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggsignif)
library(cowplot)
library(reshape2)
library(ggplot2)
library(ggforce)
library(R.utils)

## Set variables  ---------------------------------------------------------------------

FIG_DIR <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/figures/'
dir.create(paste0(FIG_DIR),  recursive = TRUE, showWarnings = FALSE)

## Load data  -------------------------------------------------------------------------

interaction <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma/snRNAseq_SCZ.GO_term_genes_inter.magma.35UP_10DOWN.gsa.out'
GO_LRRK1 <- '~/Desktop/Herring_snRNAseq_2023_pipeline/results/magma/snRNAseq_SCZ.GO_term_genes_LRRK1.magma.35UP_10DOWN.gsa.out'

## Preparing data  --------------------------------------------------------------------

MAGMA_DF_GO <- read.table(GO_LRRK1, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  dplyr::select(VARIABLE, BETA, BETA_STD) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category != "L4_RORB_LRRK1") %>%
  mutate(GO_Term = ifelse(Category == 'GO:0007399', 'Nervous system development', 
                        ifelse(Category == 'GO:0042391', 'Regulation of membrane potential',
                               ifelse(Category == 'GO:0043269', 'Regulation of ion transport', 'Trans-synaptic signaling')))) 

MAGMA_DF_GO <- MAGMA_DF_GO %>%
  mutate(Category = case_when(Category %in% c("GO:0007399", "GO:0042391", "GO:0043269", "GO:0099537") ~ "GO Term"))

MAGMA_DF_LRRK1 <- read.table(GO_LRRK1, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  dplyr::select(VARIABLE, BETA, BETA_STD) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category == "L4_RORB_LRRK1") %>%
  crossing(data.frame(GO_Term = c('Nervous system development', 'Regulation of membrane potential', 'Regulation of ion transport', 'Trans-synaptic signaling')))

MAGMA_DF_inter <- read.table(interaction, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  dplyr::select(FULL_NAME, MODEL, BETA, BETA_STD) %>%
  dplyr::rename(Category = FULL_NAME) %>%
  mutate(MODEL = ifelse(MODEL == '1', 'Nervous system development', 
                    ifelse(MODEL == '2', 'Regulation of membrane potential',
                    ifelse(MODEL == '3', 'Regulation of ion transport', 'Trans-synaptic signaling'))))

colnames(MAGMA_DF_inter) <- c('Category', 'GO_Term', 'BETA', 'BETA_STD')

MAGMA_DF_inter <- MAGMA_DF_inter %>%
            mutate(Category = case_when(Category %in% c("GO:0007399", "GO:0042391", "GO:0043269", "GO:0099537") ~ "GO Term_vs_L4_RORB_LRRK1",
                                        Category %in% "L4_RORB_LRRK1" ~ "L4_RORB_LRRK1_vs_GO Term",
                                        Category %in% c("L4_RORB_LRRK1-GO:0007399", "L4_RORB_LRRK1-GO:0042391", "L4_RORB_LRRK1-GO:0043269", "L4_RORB_LRRK1-GO:0099537") ~ "L4_RORB_LRRK1",
                                        Category %in% c("INTERACT::GO:0007399::L4_RORB_LRRK1", "INTERACT::GO:0042391::L4_RORB_LRRK1", "INTERACT::GO:0043269::L4_RORB_LRRK1", "INTERACT::GO:0099537::L4_RORB_LRRK1") ~ "Interaction")) 

MAGMA_DF <- merge(MAGMA_DF_LRRK1, MAGMA_DF_inter, all = TRUE) %>% merge(MAGMA_DF_GO, all = TRUE)

MAGMA_DF$BETA <- as.numeric(as.character(MAGMA_DF$BETA))
MAGMA_DF$BETA_STD <- as.numeric(as.character(MAGMA_DF$BETA_STD))

MAGMA_DF$Category <- gsub("_", "-", MAGMA_DF$Category)

levels <- c("L4-RORB-LRRK1", "L4-RORB-LRRK1-vs-GO Term", "GO Term", "GO Term-vs-L4-RORB-LRRK1", "Interaction")

## Creating plots  --------------------------------------------------------------------

MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(y = BETA, x = factor(Category, levels))) +
  geom_bar(stat = "identity", color = 'black', fill = "yellow", position = "dodge") +
  geom_errorbar(aes(ymin=BETA-BETA_STD, ymax=BETA+BETA_STD), width = 0.2) +
  theme_bw() +
  theme(plot.margin = unit(c(1, 2.5, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 15),
        axis.title.y = element_text(colour = "#000000", size = 15),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5, hjust = 0, angle = -45),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        strip.text = element_text(size=14, face = 'bold')) +
  ylab('Beta') +
  xlab('Test') +
  ylim(0, 0.35) +
  facet_row(facets = vars(GO_Term), scales = "free", space = "free", strip.position = "top", labeller = label_wrap_gen())

## Save joint plots  ------------------------------------------------------------------

jpeg(file = paste0(FIG_DIR,'SCZ_magma_GO_term_LRRK1_interaction_herring_lvl2_plot.jpeg'), units = "in", width = 15, height = 9, res = 300)
plot(MAGMA_PLOT)
dev.off()
