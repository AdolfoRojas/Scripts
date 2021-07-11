#!/usr/bin/env Rscript
library(tidyr)
Repeated_Measures_ANOVA <- function(df_long, titulo) {
    library(ggplot2)
    library(ggpubr)
    library(rstatix)
    df_long$ID <- as.factor(df_long$ID)
    df_long$Condition <- as.factor(df_long$Condition)
    ##Summary statistics
    print(df_long %>% group_by(Condition) %>% get_summary_stats(Value, type = "mean_sd"))

    # Visualization
    bxp <- ggboxplot(df_long, x = "Condition", y = "Value") #+ geom_violin()  # , add = "point"
    # Check assumptions
    # Outliers

    outliers <- df_long %>% group_by(Condition) %>% identify_outliers(Value)
    if (dim(outliers)[1] > 0) {
        print("Se han detectado Outliers")
        print(outliers)
    }
    # Normality assumption
    normalidad <- df_long %>% group_by(Condition) %>% shapiro_test(Value)
    for (condition in normalidad$Condition) {
        if (normalidad[normalidad$Condition == condition,]$p < 0.05){
            print(paste("Los datos de la condicion ", condition, " no cumplen con el supuesto de normalidad", sep = ""))
        }
    }
    ggqqplot(df_long, "Value", facet.by = "Condition")
    # Assumption of sphericity
    res.aov <- anova_test(data = df_long, dv = Value, wid = ID, within = Condition)
    if (res.aov[[2]][3] < 0.05) {
       print("No se cumple el supuesto de esfericidad, se aplicara la correccion de Greenhouse-Geisser")
    }
    #get_anova_table(res.aov)

    # pairwise comparisons
    pwc <- df_long %>% pairwise_t_test(Value ~ Condition, paired = TRUE, p.adjust.method = "bonferroni")
    pwc

    # Visualization: box plots with p-values
    pwc <- pwc %>% add_xy_position(x = "Condition")
    grafico <- bxp + stat_pvalue_manual(pwc, linetype = 2, step.increase = 0.025) + 
    labs(subtitle = get_test_label(res.aov, detailed = FALSE), caption = get_pwc_label(pwc)) + labs(title = titulo)
    return(list(res.aov,grafico))}




asd <- read.delim("../Trabajillos/Editado_TMZ_v2.csv", sep = ";")

# valoración indirecta de la función ventricular derecha es la medida de la excursión sistólica del anillo tricúspide, 
# que conocemos por su acrónimo en inglés (TAPSE)

tapse_data <- asd[grepl("tapse", names(asd)) | grepl("ID", names(asd))]
tapse_data <- tapse_data[order(tapse_data$ID),]
tapse_data2 <- tapse_data[complete.cases(tapse_data),]
df_long <- pivot_longer(tapse_data2, cols=names(tapse_data2)[names(tapse_data2)!="ID"], names_to = "Condition", values_to = "Value")
df_long[df_long$Condition == "tapse_bl_tmz",]$Condition <- "TMZ Wash-out"
df_long[df_long$Condition == "tapse_3_tmz",]$Condition <- "TMZ 3 Months"
df_long[df_long$Condition == "tapse_bl_pl",]$Condition <- "Baseline"
df_long[df_long$Condition == "tapse_3_pl",]$Condition <- "Placebo 3 Months"
tapse_data_analisis <- Repeated_Measures_ANOVA(df_long, "TAPSE")
#write.csv(bnp_long, file = "df_long.csv", append = FALSE, quote = FALSE, sep = "\t", row.names= FALSE)

# Cuantificacion de peptido natriuretico BNP

bnp_data <- asd[grepl("bnp", names(asd))  | grepl("ID", names(asd))]
bnp_data <- bnp_data[order(bnp_data$ID),]
bnp_data2 <- bnp_data[complete.cases(bnp_data),]
bnp_long <- pivot_longer(bnp_data2, cols=names(bnp_data2)[names(bnp_data2)!="ID"], names_to = "Condition", values_to = "Value")
bnp_data_analysis <- Repeated_Measures_ANOVA(bnp_long, "BNP levels")

bnp_long <- bnp_long[bnp_long$ID != 4 & bnp_long$ID != 25,]
bnp_data_analysis <- Repeated_Measures_ANOVA(bnp_long, "BNP levels")

# Test de caminata

walk_data <- asd[grepl("6mwt", names(asd))  | grepl("ID", names(asd))]
walk_data <- walk_data[order(walk_data$ID),]
walk_data2 <- walk_data[complete.cases(walk_data) & walk_data$ID != "2",] # el individuo numero 2 no presenta informacion en una categoria haciendo alusion al valor "no"
walk_data2$X6mwt_bl_tmz <- as.numeric(walk_data2$X6mwt_bl_tmz)
walk_long <- pivot_longer(walk_data2, cols=names(walk_data2)[names(walk_data2)!="ID"], names_to = "Condition", values_to = "Value")
walk_data_analysis <- Repeated_Measures_ANOVA(walk_long, "Walk test")

# BORG 

borg_data <- asd[grepl("borg", names(asd))  | grepl("ID", names(asd))]
borg_data <- borg_data[order(borg_data$ID),]
borg_data2 <- borg_data[complete.cases(borg_data),]
borg_long <- pivot_longer(borg_data2, cols=names(borg_data2)[names(borg_data2)!="ID"], names_to = "Condition", values_to = "Value")
borg_data_analysis <- Repeated_Measures_ANOVA(borg_long, "Borg dyspnea index")

###############################################
# Datos de practica



data("selfesteem", package = "datarium")
colnames(selfesteem)[colnames(selfesteem) == "id"] <- "ID"
selfesteem <- pivot_longer(selfesteem, cols=names(selfesteem)[names(selfesteem)!="ID"], names_to = "Condition", values_to = "Value")
selfesteem_analysis <- Repeated_Measures_ANOVA(selfesteem, "selfesteem")

asd <- read.delim("../Trabajillos/ejemplo.csv", sep = ";")
ejemplo <- asd[complete.cases(asd),]
ejemplo2 <- pivot_longer(ejemplo, cols=names(ejemplo)[names(ejemplo)!="ID"], names_to = "Condition", values_to = "Value")
ejemplo_analisis <- Repeated_Measures_ANOVA(ejemplo2, "ejemplo2")

anovaModelRM = aov(lm(Value ~ Condition + Error(ID), ejemplo2))
shapiro.test(anovaModelRM$ID$residuals)


ggsave("tapse_data_analisis.png", plot = tapse_data_analisis[[2]], width = 9, height = 8.5, dpi = 300, units = "in")

asd <- read.csv("https://raw.githubusercontent.com/marsja/jupyter/master/Python_ANOVA/rmAOV1way.csv")
head(asd)
names(asd)[1] <- "ID"
names(asd)[2] <- "Value"
names(asd)[3] <- "Condition"
df_long <- asd
ejemplo_analisis <- Repeated_Measures_ANOVA(df_long, "ejemplo2")
