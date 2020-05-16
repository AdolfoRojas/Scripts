#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ddCt")
library(Biobase)
library(lattice)
library(RColorBrewer)
library(ddCt)
datadir <- function(x) system.file("extdata", x, package="ddCt")
savedir <- function(x) file.path(tempdir(), x)
file.names <- datadir("Experiment3.txt") # SDMFrame Data object which holds a data set containing columns with the following names: 'Ct','Sample','Detector','Platename'. The column 'Ct' must contain numeric values.
#info <- datadir("sampleData.txt")
warningFile <- savedir("warnings.txt")
a<-read.delim(datadir("Experiment3.txt"))
### 3.2 Reference sample and housekeeping gene
name.reference.sample <- c("B")  ## Calibradores, no tratados
name.reference.gene <- c("u6") ## housekeeping gene
CtData <- QuantStudioFrame(file.names)
#sampleInformation <- read.AnnotatedDataFrame(info,header=TRUE, row.names=NULL)
result <- ddCtExpression(CtData, calibrationSample=name.reference.sample, housekeepingGene=name.reference.gene)
CtErr(result)
br <- errBarchart(result)
print(br)
?errBarchart()
