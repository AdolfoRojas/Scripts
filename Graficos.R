  Breast_Cancer_data <- read.delim(file = "Disease.tsv", header = T, sep = "\t")
# Load ggplot2
library(ggplot2)
library(dplyr)
Selection <- Breast_Cancer_data$Selection
Selection <- mapvalues(Selection, from = c("Hybrid Selection","CAGE","PCR","other","RT-PCR","unspecified","Inverse rRNA","Oligo-dT","size fractionation","RANDOM","PolyA"), to =c(rep("Others", 11)) )

# Create Data
df <- as.data.frame(table(Selection, dnn = list("Selection.Method")), responseName = "Experiments")

# Compute the position of labels 
df <- df %>% 
  arrange(desc(Selection.Method)) %>%
  mutate(prop = df$Experiments/sum(df$Experiments)*100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# Compute the position of labels 2
df <- df %>% 
  arrange(desc(Selection.Method)) %>%
  mutate(prop = df$Experiments/sum(df$Experiments)*100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

#lapply(df$prop, round(digits = 3))
#sapply(df$prop, mean))

# Basic piechart
par(mfrow=c(1,2))
ggplot(df, aes(x="", y=prop, fill=Selection.Method)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="left") +
  ggtitle(label = "Selection methods for library preparation") +
  theme(plot.title = element_text(color="black", size=25, face="bold.italic")) +
  geom_text(aes(y = ypos, label = paste(round(df$prop, digits = 1), "%",sep = "")), color = "white", size=6)

