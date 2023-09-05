# set a default repository for packages                         
local({r <- getOption('repos')                                                                                                                                
r['CRAN'] <- 'http://www.stats.bris.ac.uk/R/'                                                                                                          
options(repos=r)                    
}) 

list.of.cran.packages <- c("tidyverse", "Matrix", "readxl", "Seurat", "cowplot", "scCustomize", "rmarkdown", "data.table", "BiocManager", "ggdendro", "ggplot2")

# Install packages from list.of.cran.packages
#for (pkg in list.of.cran.packages){
#  if(!eval(bquote(require(.(pkg))))){
#    eval(bquote(install.packages(.(pkg), dependencies=T)))
#  }
#}

# Load packages from list.of.cran.packages
for (pkg in list.of.cran.packages){
  eval(bquote(library(.(pkg))))
}

list.of.bioc.packages <- c("EWCE", "AnnotationDbi", "org.Hs.eg.db", "scuttle", "biomaRt"))


#for (pkg in list.of.bioc.packages){
#  if(!eval(bquote(require(.(pkg))))){
#    eval(bquote(BiocManager::install(.(pkg), dependencies=T)))
#  }
#}

# Load packages from list.of.bioc.packages
for (pkg in list.of.bioc.packages){
  eval(bquote(library(.(pkg))))
}

#Set variables
RES_DIR <- '/scratch/c.c1837163/Herring_snRNAseq_2023/results/'
R_DIR <- paste0(RES_DIR, 'R_objects/')
CTD_DIR <- paste0(RES_DIR, 'ctd_objects/')


#Get MHC genes
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
mhc_genes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                   filters = c("chromosome_name","start","end"), 
                   values = list(chromosome = "6", start = "28510120", end = "33480577"),
                   mart = mart)
mhc_genes_uniq <- stringi::stri_remove_empty(unique(mhc_genes$hgnc_symbol), na_empty = FALSE)
cat('\n\nMHC genes:', length(mhc_genes_uniq), '\n')


#Create ctd object

for(ROBJ in "") { #c("", "_dwnSmpl_lvl1", "_dwnSmpl_lvl2")) {
  
  cat('\n\nCreating CTD object for:', paste0('seurat_herring', ROBJ))
  
  SEURAT_OBJ <- readRDS(paste0(R_DIR, 'seurat_herring', ROBJ, '.rds'))
  
  #Remove MHC genes
  RAW_COUNTS <- SEURAT_OBJ@assays$RNA@counts
  RAW_COUNTS_NO_MHC <- RAW_COUNTS[!(rownames(RAW_COUNTS) %in% mhc_genes_uniq), ]
  cat('\nMHC genes removed:', dim(RAW_COUNTS)[1] - dim(RAW_COUNTS_NO_MHC)[1])
  
  #Create annotations
  annotations <- as.data.frame(cbind(as.vector(rownames(SEURAT_OBJ@meta.data)),
                                     as.vector(SEURAT_OBJ@meta.data$major_clust), 
                                     as.vector(SEURAT_OBJ@meta.data$sub_clust)))
  colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
  rownames(annotations) <- NULL
  annotLevels <- list(level1class = annotations$level1class,
                      level2class = annotations$level2class)

  #Normalise - optional 
  #cat('\nRunning SCT ... ', '\n\n')
  #COUNTS_SCT <- EWCE::sct_normalize(RAW_COUNTS_NO_MHC) #Takes a very long time
  
  cat('\nRunning CPM ... ', '\n')
  COUNTs_CPM <- edgeR::cpm(RAW_COUNTS_NO_MHC)
  
  #Drop uninformative genes - excludes genes with very low or sporadic gene expression levels
  #cat('\nDropping uninformative genes raw ... ', '\n\n')
  #  DROP_GENES_RAW <- EWCE::drop_uninformative_genes(
  #  exp = RAW_COUNTS_NO_MHC, 
  #  input_species = "human",
  #  output_species = "human",
  #  level2annot = annotLevels$level2class) 
  
  #cat('\nDropping uninformative genes sct norm ... ', '\n\n')
  #  DROP_GENES_SCT <- EWCE::drop_uninformative_genes(
  #  exp = COUNTS_SCT, 
  #  input_species = "human",
  #  output_species = "human",
  #  level2annot = annotLevels$level2class) 
  
  cat('\nDropping uninformative genes cpm norm ... ', '\n\n')
  DROP_GENES_CPM <- EWCE::drop_uninformative_genes(
    exp = COUNTs_CPM, 
    input_species = "human",
    output_species = "human",
    level2annot = annotLevels$level2class,
    no_cores = 4) 
  
  cat('\nGene counts:',
      '\n\nRAW:', dim(RAW_COUNTS)[1],
      '\nRAW_NO_MHC:', dim(RAW_COUNTS_NO_MHC)[1],
      #      '\nNO_NORM_DROP_GENES:', dim(DROP_GENES_RAW)[1],
      #     '\nSCT_DROP_GENES:', dim( DROP_GENES_SCT)[1],
      '\nCPM_DROP_GENES:', dim(DROP_GENES_CPM)[1])
  
  #Create object - saves ctd obj to folder
  cat('\nCreating CTD object ... \n\n')
  ctd <- EWCE::generate_celltype_data(exp = DROP_GENES_CPM, 
                                      annotLevels = annotLevels, 
                                      groupName = paste0('herring', ROBJ),
                                      savePath = paste0(CTD_DIR),
                                      numberOfBins = 10)
}