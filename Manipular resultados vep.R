#library(VennDiagram)

file <- "VEP_p-Value_threshold_0.01_hapmap3_all_variant_effect"

rs_inicial <- read.table(file, header=F, sep = "\t")
names(rs_inicial) <- c("rsid", "Location", "Allele_alt", "Affected_gene", "Feature", "Feature_type", "Variant_type","cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Other_characteristics")

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(rs_inicial[rs_inicial$Variant_type == "intergenic_variant",]$rsid,
  rs_inicial[rs_inicial$subtype == "Her2" & rs_inicial$log2FoldChange > 0,]$rsid,
  rs_inicial[rs_inicial$subtype == "Luminal A" & rs_inicial$log2FoldChange > 0,]$rsid, 
  rs_inicial[rs_inicial$subtype == "Luminal B" & rs_inicial$log2FoldChange > 0,]$rsid),
  category.names = c("Intergenic variants","Her2", "Luminal A" , "Luminal B"),
  filename = 'Venn_DE_Protein_coding_up.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1920 , 
  resolution = 600,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans")





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
