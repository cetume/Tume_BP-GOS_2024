save_dir='../resources/raw_data/herring_2022/'

from pathlib import Path
Path(save_dir+"/data_for_R").mkdir(parents=True, exist_ok=True)
print(save_dir+"/data_for_R")

import scanpy as sc
from scipy import io

adata = sc.read(save_dir+'/RNA-all_full-counts-and-downsampled-CPM.h5ad')

#Save count matrix
io.mmwrite(save_dir+'/data_for_R/herring_counts.mtx',adata.X)

#Save metadata files
cell_meta=adata.obs.copy()
cell_meta['Barcode']=cell_meta.index
gene_meta=adata.var.copy()
gene_meta['GeneName']=gene_meta.index

cell_meta.to_csv(save_dir+'/data_for_R/herring_counts_cellMeta.csv',index=None)
gene_meta.to_csv(save_dir+'/data_for_R/herring_counts_geneMeta.csv',index=None)

