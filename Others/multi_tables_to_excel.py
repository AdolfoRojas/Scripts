#!/usr/bin/env python3
import pandas as pd
import openpyxl
import xlsxwriter
import os 
from os import walk

Output_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
tab_files = Output_dir_files.loc[Output_dir_files.file.str.contains("_variants.tab")].copy()

# We'll define an Excel writer object and the target file
Excelwriter = pd.ExcelWriter("variantes_exomas.xlsx",engine="xlsxwriter")

#We now loop process the list of dataframes
for i in tab_files.file:
    df = pd.read_csv(i, sep = "\t")
    df.to_excel(Excelwriter, sheet_name=i.replace(".tab", ""),index=False)

#And finally save the file
Excelwriter.save()
os.system("scp variantes_exomas.xlsx adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Otros_proyectos/bizamaExoma/")