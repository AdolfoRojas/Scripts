#!/usr/bin/env python3
#from abc import ABCMeta
import os
import pandas as pd
import sys
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
info[C] = info[C].replace(['miRNA-mRNA'],0.5)
info[C] = info[C].replace(to_replace='\w+',value = 0.1, regex=True) # Redes co-expresion
###########################


#info = pd.concat([info, info2], ignore_index = True, axis = 0).drop_duplicates()
info3 = pd.crosstab(info[A],info[B],values= info[C], aggfunc='sum').fillna(0) #,values= info["int_type"]).fillna(0)
print(info3)

#info3.to_csv(out_file, sep="\t",header=True, index=True)

#info3 = pd.read_csv("matrix_final_interaction.tab", sep="\t", header=0)
info2 = info3.sum(axis=1)


VEP_file = pd.read_csv("/run/media/vinicius/run-projects/Adolfo/Resultados_Tesis/Objetivo_1/VEP_p-Value_threshold_1_hapmap3_all_variant_effect", 
sep="\t", header= None)
VEP_file = VEP_file.drop(columns= VEP_file.columns[7:13]).drop(columns = [1]) #  remove : cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation and Location
#Uploaded_variation		Allele	Gene	Feature	Feature_type	Consequence		Extra
VEP_file["protein_external_id"] = VEP_file[13].str.split('SYMBOL=', n=-1, expand = True)[1].str.split(';', n=-1, expand = True)[1]





snp_af_genes = "TP53", "BRCA1", "BRCA2", "A2ML1"
               
for afected_gene in snp_af_genes:
    efecto_ponderado = info2[afected_gene]
    asd = info3.loc[[afected_gene]][info3.loc[[afected_gene]].gt(0)].dropna(axis = 1).columns
    for interactor_1st_grade in asd:
        try:
            efecto_ponderado += info2[interactor_1st_grade] * info3.loc[afected_gene,interactor_1st_grade]     
            asdf = info3.loc[[interactor_1st_grade]][info3.loc[[interactor_1st_grade]].gt(0)].dropna(axis = 1).columns
            for interactor_2nd_grade in asdf:
                try:
                    efecto_ponderado += info2[interactor_2nd_grade] * info3.loc[afected_gene,interactor_1st_grade] * info3.loc[interactor_1st_grade,interactor_2nd_grade]
                except KeyError:
                    print("Key error in " + afected_gene + " 2nd grade interactor " + interactor_2nd_grade)
                    pass
        except KeyError:
            print("Key error in " + afected_gene + " 1st grade interactor " + interactor_1st_grade)
            pass
    print(efecto_ponderado)
    
#ZCCHC9 1 150 634 408
    