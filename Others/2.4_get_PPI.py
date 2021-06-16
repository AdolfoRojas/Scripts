#!/usr/bin/env python3
##################################################################
import os
from Bio import Entrez
import re 
import requests ## python -m pip install requests
import pandas as pd
##################################################################
info  = pd.read_csv("9606.protein.info.v11.0.txt", sep="\t", header=0)
info["protein_external_id"] = info["protein_external_id"].str.split('.', n=-1, expand = True)[1]
interactions_STRING =pd.read_csv("interactions_above_400_combined_score.txt", sep = " ", header = 0)
interactions_STRING = interactions_STRING[interactions_STRING["combined_score"] > 400] #### seleccionar las interacciones de alta confiabilidad
interactions_STRING["protein1"] = interactions_STRING["protein1"].str.split('.', n=-1, expand = True)[1]
interactions_STRING["protein2"] = interactions_STRING["protein2"].str.split('.', n=-1, expand = True)[1]
interactions_number = pd.DataFrame(interactions_STRING["protein1"].value_counts()).rename(columns={"protein1": "interactions_number_over_400"})
info = info.merge(interactions_number, right_on = interactions_number.index, left_on= "protein_external_id") 
protein_IDs = pd.read_csv("protein_IDs_gencode_v35.txt", sep=";", header=None).rename(columns={0: "gene_id", 1: "transcript_id", 2: "protein_external_id"})
protein_IDs["gene_id"] = protein_IDs["gene_id"].str.split('\t', n=-1, expand = True)[8].str.split('"', n=-1, expand = True)[1].str.split('.', n=-1, expand = True)[0]
protein_IDs["transcript_id"] = protein_IDs["transcript_id"].str.split('"', n=-1, expand = True)[1].str.split('.', n=-1, expand = True)[0]
protein_IDs["protein_external_id"] = protein_IDs["protein_external_id"].str.split('"', n=-1, expand = True)[1].str.split('.', n=-1, expand = True)[0]
expressed = pd.read_csv("breast_expressed_protein_coding_transcripts.txt", sep="\t", header=0).drop(columns = {"gene_type","transcript_type"})  ### 53840 rows x 4 columns
expressed_protein_IDs = expressed.merge(protein_IDs, on = ["transcript_id","gene_id"])   #### 53840 rows x 5 columns
string = expressed_protein_IDs.merge(info, on = "protein_external_id").drop(columns = {"preferred_name"})  #####  14835 rows x 7 columns
not_protein_id_merge = pd.DataFrame(list(set(expressed_protein_IDs["transcript_id"]) - set(string["transcript_id"]))).rename(columns={0: "transcript_id"})
not_protein_id_merge = not_protein_id_merge.merge(expressed_protein_IDs, on = 'transcript_id').drop(columns = {"protein_external_id"})
info = info.sort_values(by='interactions_number_over_400', ascending=False).drop_duplicates(subset=['preferred_name'])
info = info.merge(protein_IDs, on = "protein_external_id").drop(columns = {"transcript_id"}).sort_values(by='interactions_number_over_400', ascending=False).drop_duplicates(subset=['gene_id'])
isoforms = not_protein_id_merge.merge(info, on = "gene_id").drop(columns = {"preferred_name"}).drop_duplicates()
all = string.append(isoforms)
not_protein_id_merge = pd.DataFrame(list(set(expressed_protein_IDs["transcript_id"]) - set(all["transcript_id"]))).rename(columns={0: "transcript_id"})
not_protein_id_merge = not_protein_id_merge.merge(expressed_protein_IDs, on = 'transcript_id').drop(columns = {"protein_external_id"})
all.duplicated(subset=['transcript_name']).value_counts()
gene2ensembl  = pd.read_csv("human.name_2_string.tsv", sep="\t", header=None, skiprows = [0]).drop(columns = {0}).rename(columns={1: "gene_name", 2: "protein_external_id"})
gene2ensembl["protein_external_id"] = gene2ensembl["protein_external_id"].str.split('.', n=-1, expand = True)[1]
not_protein_id_merge = not_protein_id_merge.merge(gene2ensembl, on = "gene_name")
info  = pd.read_csv("9606.protein.info.v11.0.txt", sep="\t", header=0)
info["protein_external_id"] = info["protein_external_id"].str.split('.', n=-1, expand = True)[1]
info = info.merge(interactions_number, right_on = interactions_number.index, left_on= "protein_external_id") 
isoforms = not_protein_id_merge.merge(info, on = "protein_external_id").drop(columns = {"preferred_name"}).drop_duplicates()
all = all.append(isoforms)

not_protein_id_merge = pd.DataFrame(list(set(expressed_protein_IDs["transcript_id"]) - set(all["transcript_id"]))).rename(columns={0: "transcript_id"})
not_protein_id_merge = not_protein_id_merge.merge(expressed_protein_IDs, on = 'transcript_id').drop(columns = {"protein_external_id", "gene_name"})
STRING_query_gene_ids = not_protein_id_merge["gene_id"].drop_duplicates()
print(len(STRING_query_gene_ids))
count = 0
final_df = pd.DataFrame()
Interactors= {}
for Gene_id in STRING_query_gene_ids:
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "interaction_partners" 
    request_url = "/".join([string_api_url, output_format, method]) ## Construct URL
    my_genes = Gene_id
    interaction = []
    params = {}        
    params = {"identifiers" : my_genes,"species" : 9606,"required_score" : 10}              
    response = requests.post(request_url, data=params) ## Call STRING
    count+=1
    print(count)
    for line in response.text.strip().split("\n")[1:2]:
        try:
            l = line.strip().split("\t")                   
            query_name = l[2]                
            print(query_name)
            Interactors[Gene_id]= query_name
        except:
            pass
            print("Gene ID no encontrados en STRING")

df = pd.DataFrame.from_dict(Interactors, orient='index').rename(columns={0: "gene_name"})
df.index.name = 'gene_id'
df.reset_index(inplace=True)
df = df.merge(gene2ensembl, on = "gene_name")
not_protein_id_merge = not_protein_id_merge.merge(df, on = "gene_id")
isoforms = not_protein_id_merge.merge(info, on = "protein_external_id").drop(columns = {"preferred_name"})
all = all.append(isoforms)### 53502
print(all.duplicated(subset=['transcript_id']).value_counts())
print(len(expressed["transcript_id"]))
all.to_csv("4_PPI_info.tab", sep="\t",header=True, index=True)
os.system("scp 4_PPI_info.tab adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_2/4_PPI_info.tab")