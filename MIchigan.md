## Descargar datos desde el servidor

## Descomprimir archivos ZIP
~~~
cd /media/storage/Adolfo/otros/3er_Objetivo/Michigan2
python3
from os import walk
import pandas as pd
import os
import time
Sample_dir_files = pd.DataFrame(next(walk("./"), (None, None, []))[2], columns={"file"})
Sample_files = Sample_dir_files.loc[Sample_dir_files.file.str.contains(".zip")].copy()

for i in Sample_files.file:
    os.system("unzip -P 'ChTXb2Rpi7V(dX' " + i + " &")
    time.sleep(30)

quit()
~~~

# Usar Dosage convertor 

~~~
$Docker
cd working_directory/otros/3er_Objetivo/Michigan2/
for i in {1..22}; do echo "--vcfDose chr$i.dose.vcf.gz --info chr$i.info.gz --prefix CIMBA_chr$i --type plink --format 1" >> commands; done
cat commands | xargs -P6 -n10 DosageConvertor & # P4 indica 4 trabajos a la vez n12 la cantidad de argumentos del comando
rm commands
~~~
# Covertir a BED con plink2
~~~ 
for i in {1..22}; do gunzip CIMBA_chr$i.plink.dosage.gz; done
for i in {1..22}; do ./plink2 --import-dosage CIMBA_chr$i.plink.dosage --map CIMBA_chr$i.plink.map --fam CIMBA_chr$i.plink.fam --make-bed --out CIMBA_chr$i; done
for i in {2..22}; do echo "CIMBA_chr$i.bed CIMBA_chr$i.bim CIMBA_chr$i.fam" >> unir_cromosomas.txt; done
plink --noweb --make-bed --bfile CIMBA_chr1 --merge-list unir_cromosomas.txt --out data_plink/CIMBA_Imputed
R
datos <- read.delim("data_plink/CIMBA_Imputed.fam",sep=" ", header= F)
datos <- datos[c("V1","V2")]
datos$V3 <- sapply(strsplit(as.character(datos$V2),'_'), "[", 2)
datos$V4 <- sapply(strsplit(as.character(datos$V2),'_'), "[", 2)
write.table(datos, "data_plink/id_fam_correction", sep = " ", row.names=F, col.names=F, quote=F)
quit()
n
plink --bfile data_plink/CIMBA_Imputed --update-ids data_plink/id_fam_correction --make-bed --out data_plink/CIMBA_Ids_fix
plink --bfile data_plink/CIMBA_Ids_fix --make-pheno ../cases_pheno_CIMBA.txt '*' --update-sex ../sex_CIMBA.txt --make-bed --out data_plink/CIMBA_pheno_added
cp data_plink/CIMBA_pheno_added ../
cd ../
~~~