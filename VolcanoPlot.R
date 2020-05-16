library(EnhancedVolcano)
library(dplyr)
library(stringr)

file <- "z_TCGA_BRCA_DiffExp.tab.txt"
#Volcano_file <- read.table(file, header=F, sep = "\t", col.names = c("Uploaded_variation", "Location", "Allele", "Consequence", "IMPACT", "SYMBOL", "Feature_type", "Features(transcripts)", "Feature", "BIOTYPE", "Existing_variation"), skip = 1)


Volcano_file <- read.table(file, header=T, sep = "\t")
ids_to_volcano <- read.table("ids_diff_exp.txt", header=T, sep = "|")

Volcano_file$gene <- ids_to_volcano$Gene


pvalue_filtered<- na.omit(Volcano_file[Volcano_file$pvalue <= 0.05,])


downregulated_transcripts <- na.omit(pvalue_filtered[pvalue_filtered$log2FoldChange <= -1.5,])

upregulated_transcripts <-na.omit(pvalue_filtered[pvalue_filtered$log2FoldChange >= 1.5,])


make_ring



write.table(downregulated_transcripts, sep = "\t",
            file = "C:/Users/adolf/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/adolf/Scripts/downregulated_transcripts", 
            row.names = F, quote = F, col.names = F)
write.table(upregulated_transcripts, sep = "\t",
            file = "C:/Users/adolf/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/adolf/Scripts/upregulated_transcripts", 
            row.names = F, quote = F, col.names = F)


EnhancedVolcano(Volcano_file,
                lab = Volcano_file$gene,
                title = 'Expresión diferencial Cáncer de mama',
                titleLabSize = 35,
                caption = paste0(nrow(upregulated_transcripts)+nrow(downregulated_transcripts), ' Transcritos diferencialmente expresados'),
                captionLabSize = 20,
                axisLabSize = 23,
                subtitle = NULL,
                #legendVisible = F,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-12,27),
                xlab = bquote(~Log[2]~ 'fold change'),
                labCol = 'black',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 3,
                shape = c(4, 4, 4, 18),
                col = c("grey","grey","grey","red3"),
                colAlpha = 1)


# load library
library(ggplot2)

# Create test data.
data <- data.frame(
  category=c("Upregulated", "Downregulated"),
  count=c(5594, 1879)
)

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n Genes: ", data$count)

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_manual(values = c("navyblue", "red4")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

