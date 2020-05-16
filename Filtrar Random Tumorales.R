Clasificados <- read.delim("New_Dis2.tsv")
Tumoral_samples <- read.delim("Filtered_Dis2.tsv")
Random <- Tumoral_samples[grep("RANDOM", Tumoral_samples$Selection), ]
Projects <- Random$Bioproject
df3 <- as.data.frame(table(Projects, dnn = list("Projects")), responseName = "Experiments")

#####   Resultados   ####
PRJNA396019 -> single cell
PRJEB9083
PRJEB15096
PRJEB13586
PRJEB21619 * Enrichment
PRJNA358250



View()