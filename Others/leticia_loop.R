#!/usr/bin/env Rscript
muestras_dir <- "/media/storage2/leticia/tutorial/"

samples_scRNA <- c("99491", "99492","99493","99494","99495","99496") 

ref_data <- "/media/storage2/software/refdata-gex-mm10-2020-A"






for (muestra in samples_scRNA){
    print(paste("Ejecutando Cellranger en muestra: ", muestra, sep = ""))
    print(paste("cellranger count --id=run_count_", 
    "2191_", ##### este numero tambien debe ser cambiado para que corresponda con el ID correcto --id=run_count_2191_99491
    muestra,  
    " --fastqs=", muestras_dir, muestra, " --sample=", muestra, " --transcriptome=", ref_data, sep = ""))
    print("Terminado")
}
