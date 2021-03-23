library(VennDiagram)

file <- "../VEP_p-Value_threshold_1_hapmap3_all_variant_effect_non_zero"
#rs_inicial <- read.table(file, header=F, sep = "\t", col.names = c("Uploaded_variation", "Location", "Allele", "Consequence", "IMPACT", "SYMBOL", "Feature_type", "Features(transcripts)", "Feature", "BIOTYPE", "Existing_variation"), skip = 1)


rs_inicial <- read.table(file, header=T, sep = "\t")


nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE)))



rs <-rs_inicial[c(1,2,4,5,6,8,10,13)] 


grid.newpage()
draw.quad.venn(area1 = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE))), area2 = nrow(subset(rs_inicial, grepl("protein_coding",rs_inicial$BIOTYPE))),
               area3 = nrow(subset(rs_inicial, grepl("RegulatoryFeature",rs_inicial$Feature_type))), area4 = nrow(subset(rs_inicial, grepl("intergenic_variant",rs_inicial$Consequence))), 
               n12 = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("protein_coding",rs_inicial$BIOTYPE))),
               n23 = nrow(subset(rs_inicial, grepl("protein_coding",rs_inicial$BIOTYPE) & grepl("RegulatoryFeature",rs_inicial$Feature_type))), n13 = nrow(subset(rs_inicial,grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("RegulatoryFeature",rs_inicial$Feature_type))),
               n14 = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("intergenic_variant",rs_inicial$Consequence))),
               n24 = nrow(subset(rs_inicial, grepl("protein_coding",rs_inicial$BIOTYPE) & grepl("intergenic_variant",rs_inicial$Consequence))),
               n34 = nrow(subset(rs_inicial, grepl("RegulatoryFeature",rs_inicial$Feature_type) & grepl("intergenic_variant",rs_inicial$Consequence))),
               n124 = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("protein_coding",rs_inicial$BIOTYPE) & grepl("intergenic_variant",rs_inicial$Consequence))),
               n123 = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("protein_coding",rs_inicial$BIOTYPE) & grepl("RegulatoryFeature",rs_inicial$Feature_type))),
               n134  = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("RegulatoryFeature",rs_inicial$Feature_type) & grepl("intergenic_variant",rs_inicial$Consequence))),
               n234 = nrow(subset(rs_inicial, grepl("protein_coding",rs_inicial$BIOTYPE) & grepl("RegulatoryFeature",rs_inicial$Feature_type) & grepl("intergenic_variant",rs_inicial$Consequence))),
               n1234 = nrow(subset(rs_inicial, grepl("lncRNA",rs_inicial$BIOTYPE) & grepl("protein_coding",rs_inicial$BIOTYPE) & grepl("RegulatoryFeature",rs_inicial$Feature_type) & grepl("intergenic_variant",rs_inicial$Consequence))),
               category = c("lncRNAs", "Genes codificantes", "Caracter�sticas regulatorias", "Interg�nicas"), lty = "blank", 
               fill = c("skyblue", "pink1", "mediumorchid", "red"))


d<-rs_inicial[subset(rs_inicial, !grepl("lncRNA",rs_inicial$BIOTYPE) & !grepl("protein_coding",rs_inicial$BIOTYPE) & !grepl("RegulatoryFeature",rs_inicial$Feature_type)),]

df1<-dplyr::filter(rs_inicial, !grepl("lncRNA",rs_inicial$BIOTYPE))
df2<-dplyr::filter(df1, !grepl("protein_coding",df1$BIOTYPE))
df3<-dplyr::filter(df2, !grepl("RegulatoryFeature",df2$Feature_type))
df4<-dplyr::filter(df3, !grepl("intergenic_variant",df3$Consequence))
# Fuera de las categorias en el diagrama
nrow(subset(rs_inicial, !grepl("lncRNA",rs_inicial$BIOTYPE) & !grepl("protein_coding",rs_inicial$BIOTYPE) & !grepl("RegulatoryFeature",rs_inicial$Feature_type) & !grepl("intergenic_variant",rs_inicial$Consequence)))
