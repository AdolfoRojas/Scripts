#!/usr/bin/env python3
import os
import pandas as pd
A = "tf"
B = "target"
#info  = pd.read_csv("prueba2.tsv", sep="\t", header=0)
#info2  = pd.read_csv("prueba2.tsv", sep="\t", header=0)
info = pd.read_csv("5_Interacciones_TF-Target_tejido_mamario.tsv", sep="\t", header=0)
info2 = pd.read_csv("5_Interacciones_TF-Target_tejido_mamario.tsv", sep="\t", header=0)
info = info[[A, B]]
info2 = info2[[A, B]]
A = "Element_A"
B = "Element_B"
info.columns = [A, B]
info2.columns = [B, A]

info = pd.concat([info, info2], ignore_index = True, axis = 0).drop_duplicates()
info3 = pd.crosstab(info[A],info[B]).fillna(0)

os.system("ls -lh")
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.system("ls -lh")