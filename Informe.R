##### Informe Datos SRA

# Nº de individuos | Link proyecto | Nacionalidades


Informe <- read.table("Informe.tsv", sep = "\t", header = T)
# del archivco completo se necesitan la 1 (Runs), 11(Experiments), 22(BioSamples) y 25 (Bioprojects)

library(dplyr)

SampleVSRuns = distinct(Informe[, c(5,1)])
DFSampleVSRuns <- as.data.frame(table(SampleVSRuns$X.4, dnn = list("Sample")), responseName = "Runs")

SampleVSExperiments = distinct(Informe[, c(5,2)])
DFSampleVSExperiments <- as.data.frame(table(SampleVSExperiments$X.4, dnn = list("Sample")), responseName = "Experiments")

SampleVSProjects = distinct(Informe[, c(5,4)])
DFSampleVSProjects <- as.data.frame(table(SampleVSProjects$X.4, dnn = list("Sample")), responseName = "Projects")

M1 <- merge(DFSampleVSRuns,DFSampleVSExperiments, by.x = "Sample", by.y = "Sample")
M2 <- merge(M1,DFSampleVSProjects, by.x = "Sample", by.y = "Sample")

InformeDF <- Reduce(merge, list(DFSampleVSRuns, DFSampleVSExperiments, DFSampleVSProjects))
 write.table(InformeDF,file = "InformeDF", sep = "\t")

 