##### Informe Datos SRA

# Nº de individuos | Link proyecto | Nacionalidades

setwd("/home/adolfo/Escritorio/LIB/Scripts")


Informe <- read.table("informe3.tsv", sep = "\t", header = T)
# del archivco completo se necesitan la 1 (Runs), 11(Experiments), 22(BioSamples) y 25 (Bioprojects)

library(dplyr)

SampleVSRuns = distinct(Informe[, c(25,1)])
DFSampleVSRuns <- as.data.frame(table(SampleVSRuns$BioSample, dnn = list("Sample")), responseName = "Runs")

SampleVSExperiments = distinct(Informe[, c(25,11)])
DFSampleVSExperiments <- as.data.frame(table(SampleVSExperiments$BioSample, dnn = list("Sample")), responseName = "Experiments")

SampleVSProjects = distinct(Informe[, c(25,22)])
DFSampleVSProjects <- as.data.frame(table(SampleVSProjects$BioSample, dnn = list("Sample")), responseName = "Projects")

M1 <- merge(DFSampleVSRuns,DFSampleVSExperiments, by.x = "Sample", by.y = "Sample")
M2 <- merge(M1,DFSampleVSProjects, by.x = "Sample", by.y = "Sample")
InformeDF <- Reduce(merge, list(DFSampleVSRuns, DFSampleVSExperiments, DFSampleVSProjects))
