import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import pyscenic
import sys, getopt


#set workin dir
wdir='' 
os.chdir( wdir )

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 20

# read unfiltered data from a loom file
adata = sc.read_loom('data/Pancan_obj.loom')

# compute the number of genes per cell (computes ‘n_genes' column)
sc.pp.filter_cells(adata, min_genes=0)

#
cutoff=20
sc.pp.filter_genes(adata, min_cells=cutoff )

nCountsPerGene = np.sum(adata.X, axis=0)
nCellsPerGene = np.sum(adata.X>0, axis=0)
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

# create basic row and column attributes for the loom file:
row_attrs = {
   "Gene": np.array(adata.var_names) ,
}
col_attrs = {
   "CellID": np.array(adata.obs_names) ,
   "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
   "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( "data/Pancan_obj.v2.loom", adata.X.transpose(), row_attrs, col_attrs)

