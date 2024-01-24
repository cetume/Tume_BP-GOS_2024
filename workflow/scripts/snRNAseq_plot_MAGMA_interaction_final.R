#--------------------------------------------------------------------------------------
#
#    Plot MAGMA interaction results (betas)
#
#--------------------------------------------------------------------------------------

#GO terms includes are: GO:0003008~system process, GO:0050877~neurological system process, GO:0007267~cell-cell signaling, GO:0099537~trans-synaptic signaling, GO:0030001~metal ion transport, GO:0007399~nervous system development, GO:0048666~neuron development, GO:0034765~regulation of ion transmembrane transport, GO:0043269~regulation of ion transport, GO:0042391~regulation of membrane potential

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

MAGMA_DF_LRRK1 <- read.table(GO_LRRK1, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  dplyr::select(VARIABLE, BETA, SE) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category == "L4_RORB_LRRK1") %>%
  crossing(data.frame(GO_Term = c("Nervous system development", "Regulation of membrane potential", "Regulation of ion transmembrane transport", 
                                  "System process", "Cell-cell signaling", "Metal ion transport", "Regulation of ion transport", "Neuron development", 
                                  "Neurological system process", "Trans-synaptic signaling")))

MAGMA_DF_GO <- read.table(GO_LRRK1, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
  dplyr::select(VARIABLE, BETA, SE) %>%
  dplyr::rename(Category = VARIABLE) %>%
  filter(Category != "L4_RORB_LRRK1") %>%
  mutate(GO_Term = ifelse(Category == 'GO:0007399', 'Nervous system development', 
                   ifelse(Category == 'GO:0042391', 'Regulation of membrane potential',
                   ifelse(Category == 'GO:0034765', 'Regulation of ion transmembrane transport', 
                   ifelse(Category == 'GO:0003008', 'System process', 
                   ifelse(Category == 'GO:0007267', 'Cell-cell signaling',
                   ifelse(Category == 'GO:0030001', 'Metal ion transport',        
                   ifelse(Category == 'GO:0043269', 'Regulation of ion transport',         
                   ifelse(Category == 'GO:0048666', 'Neuron development',         
                   ifelse(Category == 'GO:0050877', 'Neurological system process', 'Trans-synaptic signaling'))))))))))

MAGMA_DF_GO <- MAGMA_DF_GO %>%
  mutate(Category = case_when(Category %in% c("GO:0003008", "GO:0050877", "GO:0007267", "GO:0099537", "GO:0030001", "GO:0007399", "GO:0048666", "GO:0034765", "GO:0043269", "GO:0042391") ~ "GO Term"))

MAGMA_DF_inter <- read.table(interaction, header = FALSE) %>%
  janitor::row_to_names(row_number = 1) %>%
  dplyr::rename(Category = FULL_NAME) %>%
  mutate(GO_Term = ifelse(MODEL == '1', 'Nervous system development', 
                   ifelse(MODEL == '2', 'Regulation of membrane potential',
                   ifelse(MODEL == '3', 'Regulation of ion transmembrane transport', 
                   ifelse(MODEL == '4', 'Trans-synaptic signaling',
                   ifelse(MODEL == '5', 'System process',       
                   ifelse(MODEL == '6', 'Cell-cell signaling',      
                   ifelse(MODEL == '7', 'Metal ion transport',     
                   ifelse(MODEL == '8', 'Regulation of ion transport', 
                   ifelse(MODEL == '9', 'Neuron development', 'Neurological system process')))))))))) %>% 
  dplyr::select(Category, GO_Term, BETA, SE)

MAGMA_DF_inter <- MAGMA_DF_inter %>%
  mutate(Category = case_when(grepl("INTERACT::", Category) ~ "Interaction",
                              grepl("GO:", Category) ~ "GO Term-cond-L4-RORB-LRRK1",
                              grepl("L4_RORB_LRRK1", Category) ~ "L4-RORB-LRRK1-cond-GO Term"))

df_list <- list(MAGMA_DF_LRRK1, MAGMA_DF_GO, MAGMA_DF_inter)
MAGMA_DF <- df_list %>% reduce(full_join)

MAGMA_DF$BETA <- as.numeric(as.character(MAGMA_DF$BETA))
MAGMA_DF$SE <- as.numeric(as.character(MAGMA_DF$SE))

MAGMA_DF$Category <- gsub("_", "-", MAGMA_DF$Category)

levels <- c("L4-RORB-LRRK1", "L4-RORB-LRRK1-cond-GO Term", "GO Term", "GO Term-cond-L4-RORB-LRRK1", "Interaction")

## Creating plots  --------------------------------------------------------------------

MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(y = BETA, x = factor(Category, levels))) +
  geom_bar(stat = "identity", color = 'black', fill = "red", position = "dodge") +
  geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0.2) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 3, 0.5, 0.5), "cm"),
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
  #ylim(0, 0.325) +
  facet_wrap(~factor(GO_Term, c('Nervous system development', 'Neuron development', 'System process', 'Neurological system process', 'Metal ion transport', 
                                'Regulation of ion transport', 'Regulation of ion transmembrane transport', 'Trans-synaptic signaling', 'Cell-cell signaling',
                                'Regulation of membrane potential')), nrow = 2, scales = "fixed", strip.position = "top", labeller = label_wrap_gen())

## Save joint plots  ------------------------------------------------------------------

jpeg(file = paste0(FIG_DIR,'SCZ_magma_GO_term_LRRK1_interaction_herring_lvl2_plot_all.jpeg'), units = "in", width = 15, height = 11, res = 300)
plot(MAGMA_PLOT)
dev.off()
