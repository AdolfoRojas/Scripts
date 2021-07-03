#!/usr/bin/env python3
import pandas as pd
def Replace_lncRNA_to_matrix_ID():
    print("Iniciando Carga")        
    raw_matrix = pd.read_csv("IDs_matrix.txt", sep="|",header=None)  # Requerido
    ref_lncRNAs = pd.read_csv("Human_lncRNAs.bed12", sep="\t",header=None) # Requerido
    ref_lncRNAs[3] = ref_lncRNAs[3].str.split(".", n=1, expand = True)[0]
    raw_matrix[0] = raw_matrix[0].str.split(".", n=1, expand = True)[0]
    print("Iniciando Reemplazo")
    
    count = 0
    Replaced_types = pd.DataFrame({"ID_transcript":[],"Type_replaced":[]})

    for lncRNA in range(0, len(raw_matrix[0])):
        count += 1
        print(str(round(count/len(raw_matrix[0])*100, 3)) + "% Completado")
        if raw_matrix[7][lncRNA] != "protein_coding":        
            query = ref_lncRNAs[ref_lncRNAs[3] == raw_matrix[0][lncRNA]]       
            if query.empty != True:
                in_loop_df = pd.DataFrame({"ID_transcript":[raw_matrix[0][lncRNA]],"Type_replaced":[raw_matrix[7][lncRNA]]})
                Replaced_types = Replaced_types.append(in_loop_df)                
                raw_matrix[7][lncRNA] = "lncRNA"
    raw_matrix = raw_matrix.rename(columns= {0:"Transcripts"})
    Summary_replaced = Replaced_types["Type_replaced"].value_counts()
    print(Summary_replaced)
    Replaced_types.to_csv("lncRNAs_reemplazados_changed_types.txt", sep="\t",header=True, index=False)
    raw_matrix.to_csv("lncRNAs_reemplazados.txt", sep="|",header=True, index=False)
Replace_lncRNA_to_matrix_ID()

def Join_Matrix_ID_to_matrix_and_sum():
    print("Iniciando")
    data = pd.read_csv("TCGA_BRCA_tpm.tsv", sep="\t",header=0).drop(columns=["Unnamed: 0"]) # Requerido
    print("Archivo de datos cargado")
    data = data + 0.0001
    ids = pd.read_csv("lncRNAs_reemplazados.txt", sep="\t",header=0).rename(columns={"Transcripts|1|2|3|4|5|6|7|8":"Transcripts"})    
    ids = ids.join(data)
    ids.to_csv("TCGA-BRCA_tpm_+0.0001.txt", sep="\t",header=True, index=False)
Join_Matrix_ID_to_matrix_and_sum()

def Generate_tables_subtype_BRCA():
    Count_file = pd.read_csv("TCGA-BRCA_tpm_+0.0001.txt", sep="\t",header=0)
    Subtypes_file = pd.read_csv("Samples_Subtype_BRCA.tsv", sep="\t",header=0) # Requerido
    for subtype in Subtypes_file["BRCA_Subtype_PAM50"].unique():         
        SubType_df = Subtypes_file[Subtypes_file["BRCA_Subtype_PAM50"] == subtype]
        SubType_cols = ["Transcripts"]
        for Sample in SubType_df["CGHubAnalysisID"].unique():    
            Count_file = Count_file.rename(columns={Count_file.columns[Count_file.columns.str.contains(Sample[4:11], case=True, regex=True)][0]:Sample})
            SubType_cols.append(Sample)
        SubType_counts = Count_file[SubType_cols]
        SubType_counts.to_csv(subtype + "_tpm_+0.0001_counts.tab", sep="\t",header=True, index=False)
        print(subtype)
        SubType_df = SubType_df[["patient","age_at_diagnosis","pathologic_stage","BRCA_Subtype_PAM50","race"]]
        SubType_df["age_at_diagnosis"] = SubType_df["age_at_diagnosis"]/365.25
        SubType_df = SubType_df.drop_duplicates()
        print("Individuals: "+ str(len(SubType_df["patient"])) +"\nAverage Age at diagnosis, years(std): " + str(round(SubType_df["age_at_diagnosis"].mean(), 1)) + "("+ str(round(SubType_df["age_at_diagnosis"].std(), 1)) + ")")
        print(pd.DataFrame({"race":SubType_df["race"].value_counts().keys().tolist(), "# Individuals":SubType_df["race"].value_counts().tolist(), "%":(round(SubType_df["race"].value_counts()/SubType_df["race"].value_counts().sum()*100, 2)).tolist()}))
Generate_tables_subtype_BRCA()