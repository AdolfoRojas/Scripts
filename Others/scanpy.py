#!/usr/bin/env python3
import scanpy as sc
import numpy as np
import pandas as pd

sc.settings.set_figure_params(dpi=200, facecolor='white')
results_file = 'write/101010.h5ad'  # the file that will store the analysis results

#data_MICE_scRNA = sc.read_10x_h5('tutorial/combinar/run_count_101010/outs/raw_feature_bc_matrix.h5')
data_MICE_scRNA = sc.read_10x_mtx('tutorial/combinar/run_count_101010/outs/filtered_feature_bc_matrix', var_names='gene_symbols', cache=True)

data_MICE_scRNA.var_names_make_unique()
sc.pl.highest_expr_genes(data_MICE_scRNA, n_top=20,save='_myfigure.png')
data_MICE_scRNA.var['mt'] = data_MICE_scRNA.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(data_MICE_scRNA, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(data_MICE_scRNA, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save='_myfigure.png')


sc.pp.filter_cells(data_MICE_scRNA, min_genes=200)
sc.pp.filter_genes(data_MICE_scRNA, min_cells=3)

import os

os.system("scp (archivo) leticia@:(ruta)/(archivo")