#!/usr/bin/env python3
import pandas as pd
from pyliftover import LiftOver
File = "mmc4"
Input_file = pd.read_excel("Tesis/"+File+".xlsx", index_col=None, header=0, sheet_name='Sheet 1')
lo = LiftOver('hg19', 'hg38')
hg38 = []
for var in range(0, len(Input_file["Positionb"])):
    Asd = lo.convert_coordinate("chr" + Input_file["Chromosome"][var].astype(str), Input_file["Positionb"][var])    
    hg38.append(Asd[0][1])
Input_file["hg38-position"] = hg38
Input_file = Input_file.drop(columns = ["Positionb"])
Input_file["VEP"] = Input_file["Chromosome"].astype(str) + " " + Input_file["hg38-position"].astype(str) + " " + (Input_file["hg38-position"]+Input_file["Reference Allele"].str.len()-1).astype(str) + " " + Input_file["Reference Allele"] + "/" + Input_file["Effect Allele"]
Input_file = Input_file.sort_values(by=["Chromosome", "hg38-position"])
Input_file.to_csv("Tesis/" + File + "_VEP_ready.tab", sep="\t",header=True, index=False)