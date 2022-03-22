#!/usr/bin/env python3
import os
from numpy.lib.function_base import append
import pandas as pd
import sys
from pandas.core.reshape.merge import merge
from sklearn import preprocessing
from rpy2.robjects import r
#######################################################################################################################################
r('library(bigsnpr)')
r('load("../../1er_ObjetivoV2/2_entorno_validacion_cruzada.RData")')
#r('best_auc_pred <- pred_sp')
r('write.table(info_snp,file = "info_snp_tab", quote = F, col.names=T, row.names= F, sep = "\t")')
#######################################################################################################################################
gencode = pd.read_csv("gencode_genes.v38.annotation.tab", sep="\t", header=0).drop(['level','hgnc_id','ccdsid', 'exon_id', 'transcript_name', 'transcript_support_level', 'havana_transcript', 'transcript_id', 'transcript_type', 'frame', 'score', 'havana_gene', 'ont', 'exon_number', 'protein_id'], axis=1)
def Priorizacion_data():
    genesymbol_to_geneid = gencode[["gene_id", "gene_type", "gene_name"]]
    DE_data_mRNA = pd.read_csv("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_mRNAs-gene_level/all_samples_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",").rename(columns={"Unnamed: 0":"X"})
    DE_data_mRNA = DE_data_mRNA.merge(genesymbol_to_geneid, left_on= "X", right_on = "gene_id").drop_duplicates()
    DE_data_mRNA[["S_DE_all","Dir_all"]] = 0
    DE_data_mRNA.loc[(DE_data_mRNA.padj < 0.05) & (DE_data_mRNA.gene_type == "lncRNA") & ((DE_data_mRNA.log2FoldChange >= 1) | (DE_data_mRNA.log2FoldChange <= -1)), "S_DE_all"] = 1
    DE_data_mRNA.loc[(DE_data_mRNA.padj < 0.05) & (DE_data_mRNA.gene_type != "lncRNA") & ((DE_data_mRNA.log2FoldChange >= 2) | (DE_data_mRNA.log2FoldChange <= -2)), "S_DE_all"] = 1
    DE_data_mRNA.loc[(DE_data_mRNA.S_DE_all == 1) & (DE_data_mRNA.log2FoldChange < 0), "Dir_all"] = -1
    DE_data_mRNA.loc[(DE_data_mRNA.S_DE_all == 1) & (DE_data_mRNA.log2FoldChange > 0), "Dir_all"] = 1
    DE_data_mRNA_final = DE_data_mRNA[["X","gene_type","S_DE_all","Dir_all"]]
    subtipos = list(["Basal", "Her2", "LumB", "LumA"])
    for i in subtipos:
        print(i)
        DE_data_mRNA = pd.read_csv("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_mRNAs-gene_level/" + i + "_samples_all_controls_mRNAs_Normal_vs_Tumoral_DE.tab", sep = ",").rename(columns={"Unnamed: 0":"X"})
        DE_data_mRNA = DE_data_mRNA.merge(genesymbol_to_geneid, left_on= "X", right_on = "gene_id").drop_duplicates()    
        DE_data_mRNA[["S_DE","Dir"]] = 0
        DE_data_mRNA.loc[(DE_data_mRNA.padj < 0.05) & (DE_data_mRNA.gene_type == "lncRNA") & ((DE_data_mRNA.log2FoldChange >= 1) | (DE_data_mRNA.log2FoldChange <= -1)),"S_DE"] = 1
        DE_data_mRNA.loc[(DE_data_mRNA.padj < 0.05) & (DE_data_mRNA.gene_type != "lncRNA") & ((DE_data_mRNA.log2FoldChange >= 2) | (DE_data_mRNA.log2FoldChange <= -2)), "S_DE"] = 1
        DE_data_mRNA.loc[(DE_data_mRNA.S_DE == 1) & (DE_data_mRNA.log2FoldChange < 0), "Dir"] = -1
        DE_data_mRNA.loc[(DE_data_mRNA.S_DE == 1) & (DE_data_mRNA.log2FoldChange > 0), "Dir"] = 1
        DE_data_mRNA = DE_data_mRNA[["X","gene_type","S_DE","Dir"]].rename(columns={"S_DE":"S_DE_" + i, "Dir":"Dir_"+i})    
        DE_data_mRNA_final = DE_data_mRNA_final.merge(DE_data_mRNA, on=["X","gene_type"], how= "outer").drop_duplicates().fillna(value={"S_DE_all":0, "S_DE_" + i:0,"Dir_all":0,"Dir_"+i:0})
    DE_data_mRNA_final.loc[DE_data_mRNA_final.S_DE_Basal == 1, "S_DE_Basal"] = 2
    DE_data_mRNA_final["DE_Score"] = DE_data_mRNA_final[DE_data_mRNA_final.columns[DE_data_mRNA_final.columns.str.contains("S_DE")]].sum(axis=1)
    DE_data_mRNA_final["Dir_Score"] = DE_data_mRNA_final[DE_data_mRNA_final.columns[DE_data_mRNA_final.columns.str.contains("Dir_")]].sum(axis=1)
    DE_data_miRNA = pd.read_csv("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_miRNAs-gene_level/all_samples_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",").rename(columns={"Unnamed: 0":"X"}).drop_duplicates()
    DE_data_miRNA[["S_DE_all","Dir_all"]] = 0
    DE_data_miRNA.loc[(DE_data_miRNA.padj < 0.05) & ((DE_data_miRNA.log2FoldChange >= 1) | (DE_data_miRNA.log2FoldChange <= -1)), "S_DE_all"] = 1
    DE_data_miRNA.loc[(DE_data_miRNA.S_DE_all == 1) & (DE_data_miRNA.log2FoldChange < 0), "Dir_all"] = -1
    DE_data_miRNA.loc[(DE_data_miRNA.S_DE_all == 1) & (DE_data_miRNA.log2FoldChange > 0), "Dir_all"] = 1
    DE_data_miRNA_final = DE_data_miRNA[["X","S_DE_all","Dir_all"]]
    for i in subtipos:
        print(i)
        DE_data_miRNA = pd.read_csv("analisis_muestras_y_DE/Corregido_resultados_expresion_diferencial_miRNAs-gene_level/" + i + "_samples_all_controls_miRNAs_Normal_vs_Tumoral_DE.tab", sep = ",").rename(columns={"Unnamed: 0":"X"}).drop_duplicates()     
        DE_data_miRNA[["S_DE","Dir"]] = 0
        DE_data_miRNA.loc[(DE_data_miRNA.padj < 0.05) & ((DE_data_miRNA.log2FoldChange >= 1) | (DE_data_miRNA.log2FoldChange <= -1)),"S_DE"] = 1
        DE_data_miRNA.loc[(DE_data_miRNA.S_DE == 1) & (DE_data_miRNA.log2FoldChange < 0), "Dir"] = -1
        DE_data_miRNA.loc[(DE_data_miRNA.S_DE == 1) & (DE_data_miRNA.log2FoldChange > 0), "Dir"] = 1    
        DE_data_miRNA = DE_data_miRNA[["X","S_DE","Dir"]].rename(columns={"S_DE":"S_DE_" + i, "Dir":"Dir_"+i})      
        DE_data_miRNA_final = DE_data_miRNA_final.merge(DE_data_miRNA, on=["X"], how= "outer").drop_duplicates().fillna(value={"S_DE_all":0, "S_DE_" + i:0,"Dir_all":0,"Dir_"+i:0})
    DE_data_miRNA_final.loc[DE_data_miRNA_final.S_DE_Basal == 1, "S_DE_Basal"] = 2
    DE_data_miRNA_final["DE_Score"] = DE_data_miRNA_final[DE_data_miRNA_final.columns[DE_data_miRNA_final.columns.str.contains("S_DE")]].sum(axis=1)
    DE_data_miRNA_final["Dir_Score"] = DE_data_miRNA_final[DE_data_miRNA_final.columns[DE_data_miRNA_final.columns.str.contains("Dir_")]].sum(axis=1)
    DE_data_miRNA_final["gene_type"] = "miRNA"
    DE_data_mRNA_final = DE_data_mRNA_final.loc[DE_data_mRNA_final.gene_type != "miRNA"] ############################### eliminar luego
    DE_data_miRNA_final["X2"] = DE_data_miRNA_final["X"].replace('hsa-mir-', 'MIR', regex=True).str.upper().replace('HSA-LET-', 'MIRLET', regex=True)
    not_found_miRNAs = DE_data_miRNA_final.loc[~DE_data_miRNA_final["X2"].isin(gencode["gene_name"])]
    not_found_miRNAs["X2"] = not_found_miRNAs["X2"].replace("-", "", regex=True)
    not_found_miRNAs2 = not_found_miRNAs.loc[~not_found_miRNAs["X2"].isin(gencode["gene_name"])] # En su mayoria fueron eliminados o no tienen ensembl ID (16 en total)
    DE_data_miRNA_final = DE_data_miRNA_final.loc[~DE_data_miRNA_final["X"].isin(not_found_miRNAs["X"])]
    not_found_miRNAs = not_found_miRNAs.loc[~not_found_miRNAs["X"].isin(not_found_miRNAs2["X"])]
    DE_data_miRNA_final = DE_data_miRNA_final.append(not_found_miRNAs, ignore_index=True).drop_duplicates()
    DE_data_miRNA_final["X"] = DE_data_miRNA_final["X2"]
    DE_data_miRNA_final = DE_data_miRNA_final.drop(columns={"X2"}).drop_duplicates()
    DE_data_miRNA_final = DE_data_miRNA_final.merge(gencode, left_on= ["X","gene_type"], right_on=["gene_name","gene_type"]).drop(columns={"seqname","feature","start","end","gene_name","tag","source","strand"}).drop_duplicates()
    DE_data_miRNA_final["X"] = DE_data_miRNA_final["gene_id"]
    DE_data_miRNA_final = DE_data_miRNA_final.drop(columns={"gene_id"}).drop_duplicates()
    final_table = DE_data_mRNA_final.append(DE_data_miRNA_final, ignore_index=True)
    final_table["Consensus_DE_direction"] = 0
    final_table.loc[(final_table.DE_Score != 0) & (final_table.Dir_Score >= 0), "Consensus_DE_direction"]= "+"
    final_table.loc[(final_table.DE_Score != 0) & (final_table.Dir_Score < 0), "Consensus_DE_direction"]= "-"
    final_table.loc[final_table.Consensus_DE_direction == "+"]
    final_table = final_table[final_table.columns[~final_table.columns.str.contains("Dir_")]]
    ###########################################################################
    modulos_info = pd.read_csv("co-expression/normal_vs_tumoral/Tables/module.tsv", sep = "\t").rename(columns={"genes":"X", "modules": "Module"}).drop_duplicates()
    change_id = pd.read_csv("co-expression/2_gene_ID_to_gene_symbol_TCGA.tab", sep = "\t").drop_duplicates()
    modulos_info = modulos_info.merge(change_id, left_on= "X", right_on= "gene_name").drop_duplicates()[["gene_id","Module"]].rename(columns={"gene_id":"X"})
    final_table2 = final_table.merge(modulos_info, on="X", how="outer").fillna(value={"Module":"-"})
    final_table3 = final_table2
    final_table3.loc[(final_table3.gene_type =="miRNA") & (final_table3.DE_Score != 0)]
    final_table3.loc[final_table3.Module =="M3"]
    final_table2 = final_table2[final_table2.columns[~final_table2.columns.str.contains("S_DE_")]].rename(columns={"X":"Gene"}).drop(columns={"gene_type"})
    return final_table2

Tabla_de_priorizacion = Priorizacion_data()
#######################################################################################################################################
A = "element1"
B = "element2"
C = "int_type"
info = pd.read_csv("final_interaction.tab", sep="\t", header=0)
info = info[[A, B, C]]
###########################
def Fix_miRNAs_IDs(data): 
    info = data    
    info_miRNA = info.loc[info.int_type == "miRNA-mRNA"].copy()
    info_miRNA["original"] = info_miRNA.element1
    info_miRNA.element1 = info_miRNA.element1.str.upper().replace('-3P','', regex=True).replace('-5P','', regex=True).replace('HSA-MIR-', 'MIR', regex=True).replace('HSA-LET-', 'MIRLET', regex=True)
    info_miRNA.loc[info_miRNA.int_type == "miRNA-mRNA","element1"].loc[info_miRNA.loc[info_miRNA.int_type == "miRNA-mRNA","element1"].str.contains("LET")]
    not_found_miRNAs_int = info_miRNA.loc[info_miRNA.int_type == "miRNA-mRNA"].loc[~info_miRNA["element1"].isin(gencode["gene_name"])]
    not_found_miRNAs_int["element1"] = not_found_miRNAs_int["element1"].replace("-", "", regex=True)
    not_found_miRNAs_int2 = not_found_miRNAs_int.loc[~not_found_miRNAs_int["element1"].isin(gencode["gene_name"])]
    option_1 = not_found_miRNAs_int2.copy()
    option_2 = not_found_miRNAs_int2.copy()
    option_3 = not_found_miRNAs_int2.copy()
    option_1["X"] = option_1.element1 + "-1"
    option_2["X"]  = option_2.element1 + "-2"
    option_3["X"]  = option_3.element1 + "-3"
    option_1 = option_1.loc[option_1["X"].isin(gencode["gene_name"])]
    option_3 = option_3.loc[option_3["X"].isin(gencode["gene_name"])]
    option_2 = option_2.loc[option_2["X"].isin(gencode["gene_name"])]
    option_1 = option_1.append([option_2,option_3])
    not_found_miRNAs_int2 = not_found_miRNAs_int2.loc[~not_found_miRNAs_int2.element1.isin(option_1.element1)]
    option_4 = not_found_miRNAs_int2.copy()
    option_5 = not_found_miRNAs_int2.copy()
    option_6 = not_found_miRNAs_int2.copy()
    option_4["X"]  = option_4.element1 + "1"
    option_5["X"]  = option_5.element1 + "2"
    option_6["X"]  = option_6.element1 + "3"
    option_4 = option_4.loc[option_4["X"].isin(gencode["gene_name"])]
    option_5 = option_5.loc[option_5["X"].isin(gencode["gene_name"])]
    option_6 = option_6.loc[option_6["X"].isin(gencode["gene_name"])]
    option_4 = option_4.append([option_5,option_6])
    not_found_miRNAs_int2 = not_found_miRNAs_int2.loc[~not_found_miRNAs_int2.element1.isin(option_4.element1)]
    not_found_miRNAs_int2 = option_1.append(option_4)
    not_found_miRNAs_int2 = not_found_miRNAs_int2.merge(gencode, left_on="X", right_on="gene_name")[["element1","element2","int_type","gene_id","original"]]
    not_found_miRNAs_int = not_found_miRNAs_int.loc[~not_found_miRNAs_int.element1.isin(not_found_miRNAs_int2.element1)].merge(gencode, left_on="element1", right_on="gene_name")[["element1","element2","int_type","gene_id","original"]]
    not_found_miRNAs_int = not_found_miRNAs_int.append(not_found_miRNAs_int2)
    info_miRNA = info_miRNA.loc[~info_miRNA.original.isin(not_found_miRNAs_int.original)].merge(gencode, left_on="element1", right_on="gene_name")[["element1","element2","int_type","gene_id","original"]]
    info_miRNA = info_miRNA.append(not_found_miRNAs_int)
    info.loc[info.element1.isin(info_miRNA.original)] 
    info.loc[info.int_type == "miRNA-mRNA"]
    info = info.loc[~info.element1.isin(info_miRNA.original)] 
    info_miRNA.element1 = info_miRNA.gene_id
    info_miRNA = info_miRNA.drop(columns={"gene_id", "original"}).drop_duplicates()
    info = info.append(info_miRNA)
    info_miRNA["intermediaria"] = info_miRNA.element1
    info_miRNA.element1 = info_miRNA.element2
    info_miRNA.element2 = info_miRNA.intermediaria
    info_miRNA = info_miRNA.drop(columns=["intermediaria"]).drop_duplicates()
    info = info.append(info_miRNA)
    return info

info = Fix_miRNAs_IDs(info)
##############################
#     Valores de interaccion #
info[C] = info[C].replace(['ceRNA'],1)
info[C] = info[C].replace(['PPI'],2)
info[C] = info[C].replace(['TF-Target'],5)
info[C] = info[C].replace(['miRNA-mRNA'],3)
info[C] = info[C].replace(to_replace='\w+',value = 1, regex=True) # Redes co-expresion
##############################
info3 = pd.crosstab(info[A],info[B],values= info[C], aggfunc='sum').fillna(0) #,values= info["int_type"]).fillna(0)
info2 = info3.sum(axis=1)
#######################################################################################################################################
VEP_file = pd.read_csv("VEP_GRCh38_700K_snp.tsv", sep="\t", header= 0)[["#Uploaded_variation","Gene","IMPACT","BIOTYPE","DISTANCE","CADD_PHRED","Consequence"]]
variantes_usadas = r('variantes_geneticas')
VEP_file = VEP_file.loc[VEP_file["#Uploaded_variation"].isin(variantes_usadas)]
GTEx_expressed = pd.read_csv("GTEx_TFs_PPI/2_mRNAs_gene-level_TMM_GTEx_matrix_with_geneID.tab", sep="\t", header= 0,index_col = 0).index.to_frame().reset_index(drop=True).rename(columns={0:"Genes"})
afected_genes_5_kb = pd.DataFrame(VEP_file["Gene"].unique()) #### genes afectados segun 
afected_genes_5_kb = afected_genes_5_kb.loc[afected_genes_5_kb[0].isin(gencode.gene_id)].rename(columns={0:"ids"}) 
#######################################################################################################################################
def Find_interactions(Data):
    afected_genes_5_kb = Data    
    id_not_found = list() ## lista de genes a los cuales no se les encontro interaccion 
    first_interacion_lvl_score = pd.DataFrame(columns = ["Gene","Score"])       
    second_interacion_lvl_score = pd.DataFrame(columns = ["Gene","Score2"])  
    for afected_gene in afected_genes_5_kb.ids:    
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
        print(str(round(len(second_interacion_lvl_score)/len(afected_genes_5_kb)*100, 2))+"% Completado")
    return first_interacion_lvl_score, second_interacion_lvl_score, id_not_found

interaction_data = Find_interactions(afected_genes_5_kb)
first_interacion_lvl_score = interaction_data[0]
second_interacion_lvl_score = interaction_data[1]
id_not_found = interaction_data[2]
len(id_not_found)
len(afected_genes_5_kb) 
GTEx_expressed.loc[GTEx_expressed["Genes"].isin(id_not_found),]
#######################################################################################################################################
#                      Puntuacion a nivel de Gen afectado
Tabla_de_priorizacion.Consensus_DE_direction = Tabla_de_priorizacion.Consensus_DE_direction.replace(['-'],0).replace(['+'],1)
def Data_prep_1():
    result = pd.merge(first_interacion_lvl_score, second_interacion_lvl_score, on=["Gene"])
    result2 = pd.merge(VEP_file, result, on=["Gene"], how= "left")
    result2 = result2.drop(columns={"Consequence"}).drop_duplicates() ########################  .loc[result2.Consequence != "intergenic_variant"]
    result2 = pd.merge(result2, Tabla_de_priorizacion, on=["Gene"], how= "left")
    result2 = result2.fillna(value={"Score":0,"Score2":0,"Consensus_DE_direction":0,"DE_Score":0, "Module":"-"})
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
    return result3

result2 = Data_prep_1()

result3 = result2.copy()
result3["IMPACT"] = result3["IMPACT"].replace({"MODERATE": 5, "MODIFIER":1.5,"HIGH":10, "LOW":2})


result3["Score_fix"] = result3["IMPACT"] * (result3.Score +1) * (1 + result3.DE_Score * result3.Consensus_DE_direction) #* (1-(result3.DISTANCE/5000)*0.9) 
result3["Score2_fix"] = result3["IMPACT"] * (result3.Score2 +1) * (1 + result3.DE_Score * result3.Consensus_DE_direction)# * (1-(result3.DISTANCE/5000)*0.9) 
#######################################################################################################################################
#                      Puntuacion a nivel de SNP
def Data_prep_2(Data):
    result3 = Data
    result3 = result3.drop(columns={"DISTANCE", "BIOTYPE","DE_Score","Consensus_DE_direction","Module"}).drop_duplicates()
    result3 = result3.drop(columns={"IMPACT"}).drop_duplicates()
    X_train = result3[["Score_fix", "Score2_fix"]].values
    min_max_scaler = preprocessing.MinMaxScaler()
    result3[["Score_fix", "Score2_fix"]] = min_max_scaler.fit_transform(X_train)
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
    #result4["Score_fix"] = result4.Score_fix * result4.CADD_PHRED
    #result4["Score2_fix"] = result4.Score2_fix * result4.CADD_PHRED 
    #X_train = result4[["Score_fix", "Score2_fix"]].values
    #min_max_scaler = preprocessing.MinMaxScaler()
    #result4[["Score_fix", "Score2_fix"]] = min_max_scaler.fit_transform(X_train)
    return result4

result3 = Data_prep_2(result3)

result4 = result3.copy()
result4["Score_fix"] = (2 * result4.Score_fix) + 1 
result4["Score2_fix"] = (2 * result4.Score2_fix) + 1 

#######################################################################################################################################
def Evaluar_rendimiento():
    info_snp = pd.read_csv("info_snp_tab", sep="\t", header= 0).drop(columns={"beta_inf_mean_K","final_beta_auto_mean_K","chr","a0","a1","N","MAF","INFO","_NUM_ID_.ss"})
    info_snp = pd.merge(info_snp, result4, left_on=["rsid"], right_on=["#Uploaded_variation"], how= "inner").drop(columns={"#Uploaded_variation","Gene","CADD_PHRED"}).fillna(value={"Score_fix":1,"Score2_fix":1})
    info_snp.to_csv('../../1er_ObjetivoV2/info_snp_tab2', index=False, sep= "\t")
    #r('info_snp2 <- read.delim("info_snp_tab2", sep = "\t")')
    #r('info_snp2$Score_fix <- info_snp2$Score_fix * info_snp2$best_grid_nosp_mean_K')
    #r('info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_nosp_mean_K')
    #r('pred_sp4 <- big_prodVec(G, info_snp2$Score_fix, ind.row = ind.final.val, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)')
    #r('pred_sp <- big_prodVec(G, info_snp2$Score_fix, ind.row = poblacion_restante, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)')
    #r('pred_sp5 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = ind.final.val, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)')
    #r('pred_sp2 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = poblacion_restante, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)')
    #r('pred_sp6 <- big_prodVec(G, info_snp2$best_grid_nosp_mean_K, ind.row = ind.final.val, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)')
    #r('pred_sp3 <- big_prodVec(G, info_snp2$best_grid_nosp_mean_K, ind.row = poblacion_restante, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)')
    #AUC_pred_sp4 = r('AUC_pred_sp <- AUCBoot(pred_sp4, y[ind.final.val], seed = 1)')
    #AUC_pred_sp = r('AUC_pred_sp <- AUCBoot(pred_sp, y[poblacion_restante], seed = 1)')
    #AUC_pred_sp_baseline2 = r('AUC_pred_sp_baseline <- AUCBoot(pred_sp6, y[ind.final.val], seed = 1)')
    #AUC_pred_sp_baseline = r('AUC_pred_sp_baseline <- AUCBoot(pred_sp3, y[poblacion_restante], seed = 1)')
    #AUC_pred_sp5 = r('AUC_pred_sp2 <- AUCBoot(pred_sp5, y[ind.final.val], seed = 1)')
    #AUC_pred_sp2 = r('AUC_pred_sp2 <- AUCBoot(pred_sp2, y[poblacion_restante], seed = 1)')
    #d = {'Mean AUC': [AUC_pred_sp_baseline[0], AUC_pred_sp[0], AUC_pred_sp2[0]], 'AUC sd': [AUC_pred_sp_baseline[3], AUC_pred_sp[3], AUC_pred_sp2[3]], 'Mean AUC2':[AUC_pred_sp_baseline2[0],AUC_pred_sp4[0],AUC_pred_sp5[0]], 'AUC sd2':[AUC_pred_sp_baseline2[3],AUC_pred_sp4[3],AUC_pred_sp5[3]]}
    #rendimientos_finales = pd.DataFrame(data=d, index=["Baseline", "1era capa de interaccion", "2da capa de interaccion"])
    #return rendimientos_finales

Evaluar_rendimiento()
#rendimiento = Evaluar_rendimiento()
#rendimiento
#######################################################################################################################################

    








#datos_gencode_not_found_int = gencode[gencode['gene_id'].isin(id_not_found)].drop_duplicates()
#datos_gencode_not_found_int["gene_type"].value_counts()
#gencode.loc[gencode.gene_type == "miRNA"]
##gencode_info_first_interacion_lvl_score = gencode[gencode['gene_id'].isin(first_interacion_lvl_score["Gene"])]
#gencode_info_first_interacion_lvl_score["gene_type"].value_counts()
#datos_gencode_not_found_int[datos_gencode_not_found_int["gene_type"] == "protein_coding"]
#gencode_info_first_interacion_lvl_score[gencode_info_first_interacion_lvl_score["gene_type"] == "protein_coding"]
#info.loc[info["element1"] == "ATF3"]
#list(set(id_not_found) - set(gencode['gene_id']))