#!/usr/bin/env python3
import os
import sys
import pandas as pd
os.chdir('/media/run-projects/Adolfo/Datos_tesis/2do_Objetivo')

gencode = pd.read_csv("gencode_genes.v38.annotation.tab",sep="\t")
gencode = gencode[["gene_id","gene_type"]]
interacciones = pd.read_csv("final_interaction.tab",sep="\t")
element1 = interacciones[["element1","int_type"]].copy()
element2 = interacciones[["element2","int_type"]].copy()
elementos_unicos = element1.append(element2.rename(columns={"element2":"element1"})).drop_duplicates()

elementos_unicos_gencode = elementos_unicos.loc[elementos_unicos["element1"].isin(gencode["gene_id"])].drop_duplicates()
elementos_unicos_gencode = elementos_unicos_gencode.merge(gencode, left_on="element1",right_on="gene_id").drop_duplicates()
elementos_unicos_gencode.loc[elementos_unicos_gencode["gene_type"] =="lncRNA"].groupby("int_type").count()
elementos_unicos_gencode.loc[elementos_unicos_gencode["gene_type"] =="protein_coding"].groupby("int_type").count()
elementos_unicos_gencode.loc[elementos_unicos_gencode["gene_type"] =="miRNA"].groupby("int_type").count()




elementos_unicos_no_gencode = elementos_unicos.loc[~elementos_unicos["element1"].isin(gencode["gene_id"])]
elementos_unicos_no_gencode.loc[elementos_unicos_no_gencode["element1"].str.contains("hsa-")].groupby("int_type").count()

elementos_unicos_no_gencode = elementos_unicos_no_gencode.loc[~elementos_unicos_no_gencode["element1"].str.contains("hsa-")]
