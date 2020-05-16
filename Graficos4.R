Breast_Cancer_data <- read.delim(file = "New_Dis2.tsv", header = T, sep = "\t")
# Load ggplot2
library(ggplot2)
library(dplyr)
library(plyr)
Clasification <- Breast_Cancer_data$Clasificacion
Clasification <- mapvalues(Clasification, from = c("Healthy control","miRNA","Other Tissues related to breast cancer","Breast cancer serum","Sin Clasificar","Treatments","Study of Immunity cells","Circulating Tumor cells"), to =c(rep("Others", 8)) )
# Create Data
df2 <- as.data.frame(table(Clasification, dnn = list("Clasification")), responseName = "Experiments")

# Compute the position of labels
df2 <- df2 %>% 
  arrange(desc(Clasification)) %>%
  mutate(prop = df2$Experiments/sum(df2$Experiments)*100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

#lapply(df$prop, round(digits = 3))
#sapply(df$prop, mean))

# Basic piechart
par(mfrow=c(1,2))
theme(plot.title = element_text(color="black", size=25, face="bold.italic"))

ggplot(df2, aes(x="", y=prop, fill=Clasification)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="right") +
  ggtitle(label = "Type of Sample") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(y = ypos, label = paste(round(df2$prop, digits = 1), "%",sep = "")), color = "white", size=4)

