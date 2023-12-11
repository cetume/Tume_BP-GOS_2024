#--------------------------------------------------------------------------------------
#
#    Plot expression of SCHEMA genes
#
#--------------------------------------------------------------------------------------

##Load packages  ----------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(Seurat)
library(ggplot2)
library(tidyverse)

##  Set variables  --------------------------------------------------------------------
R_DIR <- 'Herring_snRNAseq_2023_pipeline/results/R_objects/'
FIG_DIR <- 'Herring_snRNAseq_2023_pipeline/results/figures/'

##  Load and normalise data  ----------------------------------------------------------
seurat.herring <- readRDS(paste0(R_DIR, 'seurat_herring.rds'))
gc()
seurat.herring <- NormalizeData(seurat.herring)

# Prepare data - creating and adding new meta data column (cell_type) to seurat object

CellsMeta = seurat.herring@meta.data
head(CellsMeta)

CellsMeta <- CellsMeta %>%
  mutate(cell_type = case_when(sub_clust %in% c("CCK_RELN", "CCK_SORCS1", "CCK_SYT6", "ID2_CSMD1", "LAMP5_CCK", "LAMP5_NDNF", "LAMP5_NOS1", "PV_SCUBE3", "PV_SST", "PV_SULF1", "PV_WFDC2", "SST_ADGRG6", "SST_B3GAT2", "SST_BRINP3", "SST_CALB1", "SST_NPY", "SST_STK32A", "SST_TH", "VIP_ABI3BP", "VIP_ADAMTSL1", "VIP_CHRM2", "VIP_CRH", "VIP_DPP6", "VIP_HS3ST3A1", "VIP_KIRREL3", "VIP_PCDH20") ~ "Mature Inhibitory Neurons",
                               sub_clust %in% c("CGE_dev", "ID2_dev", "MGE_dev-1", "MGE_dev-2", "PV_dev", "PV_SCUBE3_dev", "PV_SULF1_dev", "SST_ADGRG6_dev", "SST_CALB1_dev", "VIP_dev") ~ "Developing Inhibitory Neurons",
                               sub_clust %in% c("BKGR_NRGN", "L2_CUX2_LAMP5", "L3_CUX2_PRSS12", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME", "L5-6_THEMIS_CNR1", "L5-6_THEMIS_NTNG2", "L5-6_TLE4_HTR2C", "L5-6_TLE4_SCUBE1", "L5-6_TLE4_SORCS1") ~ "Mature Excitatory Neurons",
                               sub_clust %in% c("L2-3_CUX2_dev-1", "L2-3_CUX2_dev-2", "L2-3_CUX2_dev-3", "L2-3_CUX2_dev-4", "L2-3_CUX2_dev-5", "L2-3_CUX2_dev-6", "L2-3_CUX2_dev-fetal", "L2_CUX2_LAMP5_dev", "L4_RORB_dev-1", "L4_RORB_dev-2", "L4_RORB_dev-fetal", "L5-6_THEMIS_dev-1", "L5-6_THEMIS_dev-2", "L5-6_TLE4_dev", "PN_dev") ~ "Developing Excitatory Neurons",
                               sub_clust %in% c("Astro_dev-1", "Astro_dev-2", "Astro_dev-3", "Astro_dev-4", "Astro_dev-5", "Astro_GFAP", "Astro_SLC1A2", "Astro_SLC1A2_dev", "Micro", "Micro_out", "Oligo-1", "Oligo-2", "Oligo-3", "Oligo-4", "Oligo-5", "Oligo_mat", "OPC", "OPC_dev", "OPC_MBP", "Vas_CLDN5", "Vas_PDGFRB", "Vas_TBX18") ~ "Non-Neuronal"))

seurat.herring <- AddMetaData(seurat.herring, CellsMeta)

seurat.herring <- SetIdent(seurat.herring, value = seurat.herring@meta.data$sub_clust)

my_levels <- c("BKGR_NRGN", "L2_CUX2_LAMP5", "L3_CUX2_PRSS12", "L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME", "L5-6_THEMIS_CNR1", "L5-6_THEMIS_NTNG2", "L5-6_TLE4_HTR2C", "L5-6_TLE4_SCUBE1", "L5-6_TLE4_SORCS1", "L2-3_CUX2_dev-1", "L2-3_CUX2_dev-2", "L2-3_CUX2_dev-3", "L2-3_CUX2_dev-4", "L2-3_CUX2_dev-5", "L2-3_CUX2_dev-6", "L2-3_CUX2_dev-fetal", "L2_CUX2_LAMP5_dev", "L4_RORB_dev-1", "L4_RORB_dev-2", "L4_RORB_dev-fetal", "L5-6_THEMIS_dev-1", "L5-6_THEMIS_dev-2", "L5-6_TLE4_dev", "PN_dev", "CCK_RELN", "CCK_SORCS1", "CCK_SYT6", "ID2_CSMD1", "LAMP5_CCK", "LAMP5_NDNF", "LAMP5_NOS1", "PV_SCUBE3", "PV_SST", "PV_SULF1", "PV_WFDC2", "SST_ADGRG6", "SST_B3GAT2", "SST_BRINP3", "SST_CALB1", "SST_NPY", "SST_STK32A", "SST_TH", "VIP_ABI3BP", "VIP_ADAMTSL1", "VIP_CHRM2", "VIP_CRH", "VIP_DPP6", "VIP_HS3ST3A1", "VIP_KIRREL3", "VIP_PCDH20", "CGE_dev", "ID2_dev", "MGE_dev-1", "MGE_dev-2", "PV_dev", "PV_SCUBE3_dev", "PV_SULF1_dev", "SST_ADGRG6_dev", "SST_CALB1_dev", "VIP_dev", "Astro_dev-1", "Astro_dev-2", "Astro_dev-3", "Astro_dev-4", "Astro_dev-5", "Astro_GFAP", "Astro_SLC1A2", "Astro_SLC1A2_dev", "Micro", "Micro_out", "Oligo-1", "Oligo-2", "Oligo-3", "Oligo-4", "Oligo-5", "Oligo_mat", "OPC", "OPC_dev", "OPC_MBP", "Vas_CLDN5", "Vas_PDGFRB", "Vas_TBX18")

legend <- c("Mature Excitatory Neurons", "Developing Excitatory Neurons", "Mature Inhibitory Neurons", "Developing Inhibitory Neurons", "Non-Neuronal")

## Plot data --------------------------------------------------------------------------
seurat.herring@active.ident <- factor(x = seurat.herring@active.ident, levels = my_levels)
seurat.herring@meta.data$cell_type <- factor(x = seurat.herring@meta.data$cell_type, levels = legend)
PLOT <- VlnPlot(seurat.herring, features = c("CACNA1G", "CUL1", "GRIN2A", "GRIA3", "HERC1", "RB1CC1", "SETD1A", "SP4", "TRIO", "XPO7"), flip = TRUE, pt.size = FALSE, stack = TRUE, raster = FALSE, split.by = "cell_type") + theme(legend.text=element_text(size=15),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y  = element_text(size = 10),
        axis.text.x  = element_text(size = 10),
        legend.position = "top") + xlab("Cell type") + scale_fill_manual(values = c('#F8766D', '#7CAE00', '#00BFC4', '#C77CFF', 'yellow'))


## Save plot --------------------------------------------------------------------------
jpeg(file = paste0(FIG_DIR, 'Expression_SCHEMA_genes_10_plot.jpeg'), units = "in", width = 19, height = 10, res = 300)
  plot(PLOT)
  dev.off()
