setwd("e:/")
file <- "Nuevo.tsv"
rs <- read.table(file, header=T, sep = "\t")
rs$Chromosome <- paste("chr",rs$Chromosome, sep = "")
rs$PositionA <- rs$Positionb
rs$PositionB <- rs$PositionA + nchar(as.character(rs$Reference.Allele))

LiftOver_file <- paste(rs$Chromosome, rs$PositionA, rs$PositionB, sep = " ")
write.table(LiftOver_file, sep = "\t",
              file = "e:/LiftOver",
              row.names = F, quote = F, col.names = F)

Lifted_coordinates <- read.table("hglft_genome_1a3ec_78b7b0.bed", header = F, sep = "\t")
rs$PositionA <- Lifted_coordinates$V2
rs$PositionB <- Lifted_coordinates$V3-1

rs$VEP <- paste(sub("chr", "", rs$Chromosome), rs$PositionA, rs$PositionB, paste(rs$Reference.Allele, rs$Effect.Allele, sep = "/"), sep = " ")


rs <- rs[c(1,2,10,11,4:9,12)]

write.table(rs$VEP, sep = "\t",
            file = "e:/file_GRCh38_3820SNP",
            row.names = F, quote = F, col.names = F)
