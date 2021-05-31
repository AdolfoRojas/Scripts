#!/usr/bin/env python3
import pandas as pd
from pandasql import *
print("Iniciando")
freq = pd.read_csv("Scripts/SABE1171.hg38.VQSR-AS.biallelic.recodecleanFBFD.Anno3.afreq", sep="\t",header=0)
print("Archivo de frecuencias cargado")
len(freq["ALT_FREQS"])
freq = freq.dropna(subset=['ALT_FREQS'])
freq = freq[freq["ALT_FREQS"] >= 0.01]
print("Filtro por MAF terminado")
len(freq["ALT_FREQS"])
Map = pd.read_csv("Scripts/SABE1171.hg38.VQSR-AS.biallelic.recodecleanFBFD.Anno3.bim", sep="\t",header=None)
print("Archivo de coordenadas cargado")
Map = Map[[1,3]]
variant_info = freq.merge(Map, left_on="ID", right_on= 1).drop(columns= [1]).rename(columns = {3: "Position", "#CHROM":"CHROM"})
variant_info["CHROM"] = 'chr' + variant_info["CHROM"]
variant_info = variant_info[["CHROM","ID","ALT_FREQS","Position"]]
del(Map)
del(freq)
Ref = pd.read_csv("Tesis/Human_lncRNAs.bed12", sep="\t",header=None)
Ref = Ref[[0,1,2,3]]
Ref = Ref.rename(columns = {0:"CHROM", 1:"Start", 2:"End", 3:"lncRNA"})
print("Referencia de lncRNAs cargada")

pysqldf = lambda q: sqldf(q, globals())
q = """
 select
    a.lncRNA
    , b.ID
    , b.ALT_FREQS 
 from
    Ref a
 inner join
    variant_info b
        on a.CHROM = b.CHROM
  where
    a.Start <= b.Position
    and a.End >= b.Position
 """
Full_list = pysqldf(q)
Summary_list = pd.DataFrame({"IDS":Full_list["lncRNA"].value_counts().keys().tolist(), "elements":Full_list["lncRNA"].value_counts().tolist()})
Full_list.to_csv("Trabajillos/lncRNAs/lncRNA_var_SABE1171_full_list.tab", sep="\t",header=True, index=False)
Summary_list.to_csv("Trabajillos/lncRNAs/lncRNA_var_SABE1171_summary_list.tab", sep="\t",header=True, index=False)