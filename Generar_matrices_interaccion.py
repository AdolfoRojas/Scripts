#!/usr/bin/env python3
#from abc import ABCMeta
import os
import pandas as pd
import sys
from pandas.core.reshape.merge import merge
from sklearn import preprocessing
from rpy2.robjects import r
r('library(bigsnpr)')
r('load("2_entorno_validacion_cruzada.RData")')
r('best_auc_pred <- pred_sp')
r('write.table(info_snp,file = "info_snp_tab", quote = F, col.names=T, row.names= F, sep = "\t")')
print(sys.argv[1])
A = "element1"
B = "element2"
C = "int_type"

#info = pd.read_csv(sys.argv[1], sep="\t", header=0)
info = pd.read_csv("final_interaction.tab", sep="\t", header=0)
#out_file = "matrix_" + sys.argv[1]
info = info[[A, B, C]]

###########################
#     Valores             #
info[C] = info[C].replace(['ceRNA'],1)
info[C] = info[C].replace(['PPI'],1)
info[C] = info[C].replace(['TF-Target'],1)
info[C] = info[C].replace(['miRNA-mRNA'],1)
info[C] = info[C].replace(to_replace='\w+',value = 1, regex=True) # Redes co-expresion
###########################
info3 = pd.crosstab(info[A],info[B],values= info[C], aggfunc='sum').fillna(0) #,values= info["int_type"]).fillna(0)
print(info3)
info2 = info3.sum(axis=1)
VEP_file = pd.read_csv("VEP_GRCh38.tsv", sep="\t", header= 0)
VEP_file = VEP_file[["#Uploaded_variation","Gene","IMPACT","BIOTYPE","DISTANCE","CADD_PHRED"]]
gencode = pd.read_csv("gencode_genes.v38.annotation.tab", sep="\t", header=0).drop(['level','hgnc_id','ccdsid', 'exon_id', 'transcript_name', 'transcript_support_level', 'havana_transcript', 'transcript_id', 'transcript_type', 'frame', 'score', 'havana_gene', 'ont', 'exon_number', 'protein_id'], axis=1)
GTEx_expressed = pd.read_csv("GTEx_TFs_PPI/2_mRNAs_gene-level_TMM_GTEx_matrix_with_geneID.tab", sep="\t", header= 0,index_col = 0)
GTEx_expressed_genes = list(GTEx_expressed.index.values)

afected_genes_5_kb = list(VEP_file["Gene"].unique())  #### genes afectados segun 
afected_genes_5_kb_GTEx_not_expressed = list(set(afected_genes_5_kb) - set(GTEx_expressed_genes)) #### genes afectados segun VEP con muy baja expresion en tejido mamario segun GTEx o miRNA
afected_genes_5_kb_GTEx_not_expressed_df = pd.DataFrame(afected_genes_5_kb_GTEx_not_expressed)
afected_genes_5_kb_GTEx_not_expressed_df = afected_genes_5_kb_GTEx_not_expressed_df.merge(gencode, left_on=0, right_on="gene_id")
afected_genes_5_kb_GTEx_not_expressed_df["gene_type"].value_counts()
afected_genes_5_kb_GTEx_expressed = list(set(afected_genes_5_kb) - set(afected_genes_5_kb_GTEx_not_expressed)) #### genes afectados segun VEP que se expresan en tejido mamario segun GTEx

id_not_found = list() ## lista de genes a los cuales no se les encontro interaccion 
first_interacion_lvl_score = pd.DataFrame(columns = ["Gene","Score"])       
second_interacion_lvl_score = pd.DataFrame(columns = ["Gene","Score2"])  

for afected_gene in afected_genes_5_kb_GTEx_expressed:    
    try:
        efecto_ponderado = info2[afected_gene]        
        interacion_lvl_score_loop = pd.DataFrame([[afected_gene, efecto_ponderado]],columns = ["Gene","Score"]) 
        first_interacion_lvl_score = first_interacion_lvl_score.append(interacion_lvl_score_loop, ignore_index=True)
        asd = info3.loc[[afected_gene]][info3.loc[[afected_gene]].gt(0)].dropna(axis = 1).columns
        for interactor_1st_grade in asd:
            try:
                efecto_ponderado += info2[interactor_1st_grade] * info3.loc[afected_gene,interactor_1st_grade]
            except KeyError:                
                pass        
        interacion_lvl_score_loop2 = pd.DataFrame([[afected_gene, efecto_ponderado]],columns = ["Gene","Score2"]) 
        second_interacion_lvl_score = second_interacion_lvl_score.append(interacion_lvl_score_loop2, ignore_index=True)      
    except KeyError:        
        id_not_found.append(afected_gene)
        efecto_ponderado = 0
        interacion_lvl_score_loop = pd.DataFrame([[afected_gene, efecto_ponderado]],columns = ["Gene","Score"]) 
        first_interacion_lvl_score = first_interacion_lvl_score.append(interacion_lvl_score_loop, ignore_index=True)
        interacion_lvl_score_loop2 = pd.DataFrame([[afected_gene, efecto_ponderado]],columns = ["Gene","Score2"]) 
        second_interacion_lvl_score = second_interacion_lvl_score.append(interacion_lvl_score_loop2, ignore_index=True) 
        pass    
    print(str(round(len(second_interacion_lvl_score)/len(afected_genes_5_kb_GTEx_expressed)*100, 2))+"% Completado")

len(afected_genes_5_kb)
len(id_not_found) 
len(first_interacion_lvl_score)

result = pd.merge(first_interacion_lvl_score, second_interacion_lvl_score, on=["Gene"])
result2 = pd.merge(VEP_file, result, on=["Gene"], how= "left")
result2 = result2.fillna(value={"Score":0,"Score2":0})
result2 = result2.drop_duplicates()
result2[["DISTANCE","CADD_PHRED"]] = result2[["DISTANCE","CADD_PHRED"]].replace({"-": 0})
result2["DISTANCE"] = pd.to_numeric(result2["DISTANCE"])
result2["CADD_PHRED"] = pd.to_numeric(result2["CADD_PHRED"])
CADD_PHRED = result2.groupby('#Uploaded_variation')['CADD_PHRED'].max().to_frame()
result2 = result2.drop(columns={"CADD_PHRED"}).drop_duplicates()
CADD_PHRED["#Uploaded_variation"] = CADD_PHRED.index
CADD_PHRED = CADD_PHRED.reset_index(drop=True)
result2 = pd.merge(result2, CADD_PHRED, on=["#Uploaded_variation"], how= "left")
result2 = result2.drop_duplicates(subset=['#Uploaded_variation', 'Gene'])

N_features = result2.groupby('#Uploaded_variation')['Gene'].count().to_frame()
N_features["#Uploaded_variation"] = N_features.index
N_features = N_features.reset_index(drop=True)
result3 = result2.drop(columns={"Gene"})
result3 = pd.merge(result3, N_features, on=["#Uploaded_variation"], how= "left")
result3["IMPACT"] = result3["IMPACT"].replace({"MODERATE": 5, "MODIFIER":1.5,"HIGH":10, "LOW":2})
result3["Score_fix"] = result3["IMPACT"] * (result3.Score +1) * (1-(result3.DISTANCE/5000)*0.9)
result3["Score2_fix"] = result3["IMPACT"] * (result3.Score2 +1) * (1-(result3.DISTANCE/5000)*0.9)
result3 = result3.drop(columns={"DISTANCE", "BIOTYPE"}).drop_duplicates()
result3 = result3.drop(columns={"IMPACT"}).drop_duplicates()
Score_fix = result3.groupby('#Uploaded_variation')['Score_fix'].sum().to_frame()
Score2_fix = result3.groupby('#Uploaded_variation')['Score2_fix'].sum().to_frame()
result3 = result3.drop(columns={"Score_fix", "Score2_fix"}).drop_duplicates()
Score_fix["#Uploaded_variation"] = Score_fix.index
Score_fix = Score_fix.reset_index(drop=True)
Score2_fix["#Uploaded_variation"] = Score2_fix.index
Score2_fix = Score2_fix.reset_index(drop=True)
result3 = pd.merge(result3, Score_fix, on=["#Uploaded_variation"], how= "left")
result3 = pd.merge(result3, Score2_fix, on=["#Uploaded_variation"], how= "left")
result3 = result3.drop(columns={"Score", "Score2"}).drop_duplicates()

result4 = result3.drop_duplicates()
result4["Score_fix"] = result4.CADD_PHRED * result4.Score_fix
result4["Score2_fix"] = result4.CADD_PHRED * result4.Score2_fix
X_train = result4[["Score_fix", "Score2_fix"]].values
min_max_scaler = preprocessing.MinMaxScaler()
result4[["Score_fix", "Score2_fix"]] = min_max_scaler.fit_transform(X_train)
result4["Score_fix"] = (2 * result4.Score_fix) + 1 
result4["Score2_fix"] = (2 * result4.Score2_fix) + 1 

info_snp = pd.read_csv("info_snp_tab", sep="\t", header= 0).drop(columns={"best_grid_nosp_mean_K","beta_inf_mean_K","final_beta_auto_mean_K"})
info_snp = pd.merge(info_snp, result4, left_on=["rsid"], right_on=["#Uploaded_variation"], how= "left").drop(columns={"#Uploaded_variation","IMPACT","Gene","CADD_PHRED"}).fillna(value={"Score_fix":1,"Score2_fix":1})
info_snp.to_csv('info_snp_tab2', index=False, sep= "\t")
r('info_snp2 <- read.delim("info_snp_tab2", sep = "\t")')
r('info_snp2$Score_fix <- info_snp2$Score_fix * info_snp2$best_grid_sp_mean_K')
r('info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K')
r('pred_sp <- big_prodVec(G, info_snp2$Score_fix, ind.row = ind.final.val, ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)')
r('pred_sp2 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = ind.final.val, ind.col = df_beta[["_NUM_ID_"]], ncores = NCORES)')
r('AUC_pred_sp <- AUCBoot(pred_sp, y[ind.final.val], seed = 1)')
r('AUC_pred_sp_baseline <- AUCBoot(best_auc_pred, y[ind.final.val], seed = 1)')
r('AUC_pred_sp2 <- AUCBoot(pred_sp2, y[ind.final.val], seed = 1)')

datos_gencode_not_found_int = gencode[gencode['gene_id'].isin(id_not_found)]
datos_gencode_not_found_int["gene_type"].value_counts()
gencode_info_first_interacion_lvl_score = gencode[gencode['gene_id'].isin(first_interacion_lvl_score["Gene"])]
gencode_info_first_interacion_lvl_score["gene_type"].value_counts()

datos_gencode_not_found_int[datos_gencode_not_found_int["gene_type"] == "protein_coding"]
gencode_info_first_interacion_lvl_score[gencode_info_first_interacion_lvl_score["gene_type"] == "protein_coding"]
info.loc[info["element1"] == "ATF3"]

list(set(id_not_found) - set(gencode['gene_id']))



