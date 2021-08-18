#!/usr/bin/env python3
import pandas as pd
from pyliftover import LiftOver

File = "p-Value_threshold_1_hapmap3_all_variant_effect_non_zero_GRCh37.txt"
Input_file = pd.read_csv(File, index_col=None, header=None, sep = " ")
lo = LiftOver('hg19', 'hg38')
Input_file[6] = ""
Input_file[7] = ""
#hg38 = []
id_not_found = list()
Asd = list()
for var in range(0, len(Input_file[1])):
    print(var)
    try:
        Asd = lo.convert_coordinate("chr" + Input_file[0][var].astype(str), Input_file[1][var])
        Asdf = lo.convert_coordinate("chr" + Input_file[0][var].astype(str), Input_file[2][var])            
        Input_file[6][var] = Asd[0][1]
        Input_file[7][var] = Asdf[0][1]
    except IndexError:        
        id_not_found.append(list([var, Input_file[5][var]]))
        #Input_file[6][var] = Asd[2]
        #Input_file[7][var] = Asdf[2]
        pass
Input_file.loc[Input_file[5] == "rs12728058"][6] = 555
Input_file.loc[Input_file[5] == "rs12728058"][7] = ""

Input_file.loc[id_not_found[5][0],]
Input_file.loc[id_not_found[4][0],6] = 142739784     
Input_file.loc[id_not_found[4][0],7] = 142739784     
Input_file.loc[id_not_found[4][0],]
#Input_file["hg38-position"] = hg38
#Input_file = Input_file.drop(columns = [1])
#Input_file["VEP"] = Input_file[0].astype(str) + " " + Input_file["hg38-position"].astype(str) + " " + (Input_file["hg38-position"]+Input_file["Reference Allele"].str.len()-1).astype(str) + " " + Input_file["Reference Allele"] + "/" + Input_file["Effect Allele"]
#Input_file = Input_file.sort_values(by=[0, "hg38-position"])
Input_file.to_csv(File + "_VEP_ready.tab", sep="\t",header=True, index=False)



id_not_found = list()
first_interacion_lvl_score = pd.DataFrame(columns = ["Gene","Score"])              
for afected_gene in afected_genes_5_kb_GTEx_expressed:
    print(afected_gene)
    try:
        efecto_ponderado = info2[afected_gene]
        interacion_lvl_score_loop = pd.DataFrame([[afected_gene, efecto_ponderado]],columns = ["Gene","Score"]) 
        first_interacion_lvl_score = first_interacion_lvl_score.append(interacion_lvl_score_loop, ignore_index=True)
    except KeyError:
        print("Key error: " + afected_gene + " not found")
        id_not_found.append(afected_gene)
        pass
    print(efecto_ponderado)