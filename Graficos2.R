Breast_Cancer_data <- read.delim(file = "New_Dis2.tsv", header = T, sep = "\t")
# Load ggplot2
library(ggplot2)
library(dplyr)
Clasificacion <- Breast_Cancer_data$Clasificacion

# Create Data
df2 <- as.data.frame(table(Clasificacion, dnn = list("Clasificacion")), responseName = "Experiments")

# Compute the position of labels
df2 <- df2 %>% 
  arrange(desc(Clasificacion)) %>%
  mutate(prop = df$Experiments/sum(df$Experiments)*100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

#lapply(df$prop, round(digits = 3))
#sapply(df$prop, mean))

# Basic piechart
par(mfrow=c(1,2))
theme(plot.title = element_text(color="black", size=20, face="bold.italic"))

ggplot(df2, aes(x="", y=prop, fill=Clasificacion)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="right") +
  ggtitle(label = "Type of Sample") +
  geom_text(aes(y = ypos, label = paste(round(df$prop, digits = 1), "%",sep = "")), color = "white", size=4)

