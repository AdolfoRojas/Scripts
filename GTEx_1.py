#!/usr/bin/env python3
import pandas as pd
input_file = input("Full file name")
female_only = input("Only female data? (True/False)")

def Prepare_info():
    sample_info = pd.read_csv("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep="\t", header=0)
    individual_info = pd.read_csv("https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",sep="\t", header=0).rename(columns = {"SUBJID":"Subject_ID", "DTHHRDY": "Death(Hardy_Scale)"})
    individual_info["SEX"] = individual_info["SEX"].replace([1,2],["Male", "Female"])
    individual_info["Death(Hardy_Scale)"] = individual_info["Death(Hardy_Scale)"].replace([0,1,2,3,4],["Ventilator Case", "Violent and fast death", "Fast death of natural causes","Intermediate death","Slow death"])
    Breast_samples = sample_info[sample_info["SMTSD"]=="Breast - Mammary Tissue"].reset_index()
    Breast_samples["Subject_ID"] = Breast_samples["SAMPID"].str.split("-", expand= True)[0] + "-" + Breast_samples["SAMPID"].str.split("-", expand= True)[1]
    Breast_samples_with_individual_info = individual_info.merge(Breast_samples, on = "Subject_ID")
    Breast_samples_with_individual_info.to_csv("GTEx_Breast_samples_info.tab", sep="\t",header=True, index=False)
#Prepare_info()

def Generate_tables_genre_Breast():
    Count_file = pd.read_csv(input_file, sep="\t", header=2, compression="gzip")
    Info_file = pd.read_csv("GTEx_Breast_samples_info.tab", sep="\t",header=0)
    ids_count_file = Count_file[["Name", "Description"]]
    for genre in Info_file["SEX"].unique():
        if female_only == True:
                if genre == "Male":
                        continue
        genre_df = Info_file[Info_file["SEX"] == genre]
        genre_cols = []
        for Sample in genre_df["SAMPID"].unique():
                if Count_file.columns[Count_file.columns.str.contains(Sample, case=True, regex=True)].empty != True:
                        genre_cols.append(Sample)
        genre_counts = Count_file[genre_cols]        
        genre_counts = ids_count_file.join(genre_counts)
        genre_counts.to_csv("1_" + genre + "_" + input_file + "_gene_level.tab", sep="\t",header=True, index=False)
        print(genre + " samples:")
        print(len(genre_counts.columns)-2)
Generate_tables_genre_Breast()