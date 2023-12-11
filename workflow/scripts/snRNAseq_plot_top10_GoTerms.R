#--------------------------------------------------------------------------------------
#
#    Plot GO Term results
#
#--------------------------------------------------------------------------------------

##Load packages  ----------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(janitor)
library(readxl)

##  Set variables  --------------------------------------------------------------------
DATA_DIR <- 'Herring_snRNAseq_2023_pipeline/resources/GO_Term_data/'
FIG_DIR <- 'Herring_snRNAseq_2023_pipeline/results/figures/'
CELL_TYPES <- c("L4_RORB_LRRK1", "L4_RORB_dev_2")


##  Load data  ------------------------------------------------------------------------
GO_DATA <- read_excel(paste0(DATA_DIR, "GOTerm_top10_data.xlsx"), skip = 2, col_names = F) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) # Get rid of NAs

colnames(GO_DATA) <- c("Term", "L4_RORB_LRRK1-FE", "L4_RORB_LRRK1-FDR", "L4_RORB_dev_2-FE", "L4_RORB_dev_2-FDR")

GO_DATA$Term <- c("GO:0007399~nervous system development",
                  "GO:0048666~neuron development",
                  "GO:0022008~neurogenesis",
                  "GO:0048699~generation of neurons",
                  "GO:0048468~cell development",
                  "GO:0048588~developmental cell growth",
                  "GO:0030182~neuron differentiation",
                  "GO:0048667~cell morphogenesis involved in neuron differentiation",
                  "GO:0000902~cell morphogenesis",
                  "GO:0031175~neuron projection development",
                  "GO:0048812~neuron projection morphogenesis",
                  "GO:0050804~modulation of synaptic transmission",
                  "GO:0099537~trans-synaptic signaling",
                  "GO:0099536~synaptic signaling",
                  "GO:0098916~anterograde trans-synaptic signaling",
                  "GO:0007267~cell-cell signaling",
                  "GO:0048167~regulation of synaptic plasticity")

# Create factor for y-axis order
ORDERED_LIST <- GO_DATA %>% dplyr::select(Term) 
  pull() %>%
  as_factor()


# Prepare data
GO_DATA <- GO_DATA %>%
  pivot_longer(-Term) %>%
  separate(name, into = c("cell_type", "Score"), '-') %>%
  pivot_wider(names_from = Score, values_from = value) %>%
  mutate(cell_type = ifelse(cell_type == "L4_RORB_dev_2", "L4_RORB_dev-2", cell_type))

colnames(GO_DATA) <- c("Term", "cell_type", "Fold Enrichment", "FDR")

GO_DATA$FDR <- as.numeric(as.character(GO_DATA$FDR))

# Plot data
GO_plot <- ggplot(data = GO_DATA, aes(y = factor(Term, level = rev(ORDERED_LIST)), x = cell_type,
                                      color = -log10(FDR), size = `Fold Enrichment`)) +
  geom_point() +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        #      panel.grid.major = element_blank(),
        #    panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 12, vjust = -0.5),
        axis.title.y = element_text(colour = "#000000", size = 12),
        axis.text.x  = element_text(colour = "#000000", size = 10, vjust = 0.8, hjust = 0, angle = -45),
        axis.text.y  = element_text(colour = "#000000", size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(colour = "#000000")) +
  ylab("Term") +
  xlab("Cell type") +
  scale_radius(limits = c(1, 9), range = c(1,10))

## Save plot --------------------------------------------------------------------------
tiff(paste0(FIG_DIR, "GoTerms_Top10.tiff"), height = 18, width = 19, units='cm',
     compression = "lzw", res = 300)
print(GO_plot)
dev.off()
