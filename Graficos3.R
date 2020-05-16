Breast_Cancer_data <- read.delim(file = "Disease.tsv", header = T, sep = "\t")
# Load ggplot2
library(ggplot2)
library(dplyr)
library(plyr)
Selection <- Breast_Cancer_data$Selection
Selection_exa <- Breast_Cancer_data$Selection
Selection <- mapvalues(Selection, from = c("Hybrid Selection","CAGE","PCR","other","RT-PCR","unspecified","Inverse rRNA","Oligo-dT","size fractionation","PolyA"), to =c(rep("Others", 10)) )
Selection <- mapvalues(Selection, from = "RANDOM", to = "RANDOM PCR")
#Selection_exa <- mapvalues(Selection_exa, from = c("Hybrid Selection","unspecified"), to = c("other","other"))
# Create Data

df <- as.data.frame(table(Selection, dnn = list("Selection.Method")), responseName = "Experiments")

df_exa <- as.data.frame(table(Selection_exa, dnn = list("Others")), responseName = "Experiments")

exa<-df_exa[c(1,3:8,11:13),]


# Compute the position of labels 
df <- df %>% 
  arrange(desc(Selection.Method)) %>%
  mutate(prop = df$Experiments/sum(df$Experiments)*100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
exa <- exa %>% 
  arrange(desc(Others)) %>%
  mutate(prop = exa$Experiments/sum(exa$Experiments)*100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
#lapply(df$prop, round(digits = 3))
#sapply(df$prop, mean))

# Basic piechart
library(cowplot)

ggplot(df, aes(x="", y=prop, fill=Selection.Method)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="left") +
  ggtitle(label = "Selection methods for library preparation") +
  theme(plot.title = element_text(color="black", size=25, face="bold.italic")) +
  geom_text(aes(y = ypos, label = paste(round(df$prop, digits = 1), "%",sep = "")), color = "white", size=6)


plot2<-ggplot(exa, aes(x="", y=prop, fill=Others)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="right") +
  #ggtitle(label = "Selection methods for library preparation") +
 # theme(plot.title = element_text(color="black", size=25, face="bold.italic")) +
  geom_text(aes(y = ypos, label = paste(round(exa$prop, digits = 1), "%",sep = "")), color = "white", size=5)
plot_grid(plot1, plot2, labels = "AUTO")
