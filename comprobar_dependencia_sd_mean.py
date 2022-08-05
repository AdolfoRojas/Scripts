#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
datos = pd.read_csv("/media/run-projects/Adolfo/Datos_tesis/2do_Objetivo/analisis_muestras_y_DE/Gene-level_mRNAs_DESeq_norm_matrix_WT.csv", sep = ",")
#datos = pd.read_csv("/media/run-projects/Adolfo/Datos_tesis/2do_Objetivo/analisis_muestras_y_DE/2_mRNAs_gene-level_TMM_TCGA_matrix_with_gene_id.tab", sep = "\t")
muestras = pd.read_csv("/media/run-projects/Adolfo/Datos_tesis/2do_Objetivo/analisis_muestras_y_DE/2_sample_annot_corregido.tab", sep ="\t")
datos = datos.loc[:,datos.columns.isin(muestras["SampleName"])]

#datos = np.log10(datos)

x = datos.mean(axis=1).to_numpy()
y = datos.std(axis=1).to_numpy()

plt.scatter(x, y, alpha=0.5)
plt.show()