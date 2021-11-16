
import pyupset as pyu
from os import walk
import os
import argparse
import pandas as pd

Output_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
table_files = Output_dir_files.loc[Output_dir_files.file.str.contains("sumary_table.tsv")].copy().reset_index(drop=True)
tabla_final = pd.DataFrame()
contador = 0
for tabla in table_files.file:#[range(0,2)]:
    contador+=1
    print(tabla)
    muestra = tabla.replace(".vcf.sumary_table.tsv","").replace("output_","")
    archivo = pd.read_csv(tabla,sep = "\t")
    archivo["Identificador"] = archivo.CHROM.map(str) + archivo.POS.map(str) + archivo.REF.map(str) + archivo.ALT.map(str) + archivo.TYPE.map(str)
    #archivo = archivo.loc[archivo.NCALLED == 1].assign(Muestra=muestra)
    archivo = archivo[["Identificador","NCALLED"]].rename(columns={"NCALLED":muestra}).set_index('Identificador')
    #archivo = archivo
    if contador == 1:
        tabla_final = tabla_final.append(archivo)
    else:
        tabla_final = pd.concat([tabla_final, archivo], axis=1)
    print(tabla_final)

tabla_final.to_csv("upset_plot_df.tsv", sep = "\t")
doublequote=True
pyu.plot(tabla_final)
pyu.plot(tabla_final, unique_keys = ['Identificador'])

R
library(UpSetR)
datos <- read.delim("upset_plot_df.tsv", sep = "\t")
datos <- datos[!duplicated(datos$Identificador),]
rownames(datos) <- datos$Identificador
datos$Identificador <- NULL
#jpeg(filename="upset_plot.jpg", width = 1920, height = 1080,res=600, bg = "white")

plot_upset <- upset(datos,order.by = "freq",nsets = length(colnames(datos)),scale.intersections= "log2",text.scale = 0.75,set_size.show = FALSE, number.angles = 20)#, scale.sets = "log2"

pdf(file="upset_plot.pdf")
plot_upset
dev.off()
system("scp upset_plot.* adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/")
