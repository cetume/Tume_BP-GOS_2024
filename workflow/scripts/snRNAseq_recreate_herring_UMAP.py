import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import math
import seaborn as sns
import os
sc.settings.verbosity = 0
sc.logging.print_header()
sns.set_context("paper")

import warnings
warnings.simplefilter( action="ignore", category=FutureWarning)

# Set paths to directories for data and figures
fig_path = "../results/figures/"

# Read in downsampled count matrices post whole-tissue clustering
adata = sc.read( "../resources/raw_data/herring_2022/RNA-all_full-counts-and-downsampled-CPM.h5ad")
adata

# Create and save UMAP
ax = sc.pl.umap(adata, color='sub_clust', size=2.5, sort_order=False, show=False)
plt.gcf().set_size_inches(5, 5)
plt.savefig( f"{fig_path}herring_UMAP.png", format='png', dpi=300, bbox_inches='tight')
