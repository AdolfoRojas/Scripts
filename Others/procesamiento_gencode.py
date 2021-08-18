#!/usr/bin/env python3
from gtfparse import read_gtf
df = read_gtf("gencode.v38.annotation.gtf")
#df = df[df["source"] != "HAVANA"]
len(df["gene_id"].unique())
len(df["gene_id"].str.split('.', n=-1, expand = True)[0].unique())
df["gene_id"] = df["gene_id"].str.split('.', n=-1, expand = True)[0]
conversion = df[["gene_name", "gene_id"]]
conversion = df[["gene_name", "gene_id"]].drop_duplicates()
df_genes = df[df["feature"] == "gene"]




df_genes[df_genes["gene_name"]=="ZNF197"]
df_genes.columns
df_genes.to_csv("gencode_genes.v38.annotation.tab", sep="\t",index=False)
