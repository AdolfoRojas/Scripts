#!/usr/bin/env Rscript

setwd("Otros_proyectos/Cohorte_TMZ/")
tmz_cohort_data <- read.delim("Editado_TMZ_v2.csv", sep = ";")
tmz_cohort_data$rut <- NULL # Dato no necesario 
tmz_cohort_data$fecha_ingreso <- NULL # Dato no necesario
tmz_cohort_data$fecha_nacimiento <- NULL # Dato no necesario
tmz_cohort_data$ano_inicio <- NULL # Dato no necesario
tmz_cohort_data$ano_diagnostico <- NULL # Dato no necesario

num <- Filter(is.numeric, tmz_cohort_data)
non_num <- tmz_cohort_data[,!(names(tmz_cohort_data) %in% names(num))]
tmz_cohort_data$emr_3_tmz <- as.numeric(as.character(tmz_cohort_data$emr_3_tmz))
tmz_cohort_data$peso_tmz3 <- as.numeric(as.character(tmz_cohort_data$peso_tmz3))
tmz_cohort_data$pad_pl2 <- as.numeric(as.character(tmz_cohort_data$pad_pl2))
tmz_cohort_data$X6mwt_bl_tmz <- as.numeric(as.character(tmz_cohort_data$X6mwt_bl_tmz))

num <- Filter(is.numeric, tmz_cohort_data)
non_num <- tmz_cohort_data[,!(names(tmz_cohort_data) %in% names(num))]

tmz_cohort_data_NA <- tmz_cohort_data[, sapply(tmz_cohort_data, function(x) sum(is.na(x))) > length(tmz_cohort_data$ID)*0.3]
tmz_cohort_data <- tmz_cohort_data[, sapply(tmz_cohort_data, function(x) sum(is.na(x))) <= length(tmz_cohort_data$ID)*0.3]
################ Imputacion ###################################################################
library(VIM)
tmz_cohort_data <- kNN(tmz_cohort_data, imp_var = FALSE)
nombres <- names(tmz_cohort_data)[grepl("_bl_pl$", names(tmz_cohort_data)) | grepl("_3_pl$", names(tmz_cohort_data)) | grepl("_3_tmz$", names(tmz_cohort_data)) | grepl("_bl_tmz$", names(tmz_cohort_data))]
library(stringr)
asdf <- str_replace_all(nombres, "_bl_pl", "")
asdf <- str_replace_all(asdf, "_3_pl", "")
asdf <- str_replace_all(asdf, "_3_tmz", "")
asdf <- str_replace_all(asdf, "_bl_tmz", "")

determinaciones <- list()
for (i in unique(asdf)) {
   if ( table(unlist(asdf))[[i]]==4) {
      determinaciones <- append(determinaciones, i)}}
determinaciones <- unlist(determinaciones)
determinaciones <- determinaciones[determinaciones!="uso_ox"]
###############################################
# T-test
Paired_T_test_or_Wilcoxon_TMZ <- function(df_long, titulo) {    
    library(ggplot2)
    library(ggpubr)
    library(rstatix)
    library(tidyr)
    df_long <- df_long[order(df_long$ID),]

    diferencia <- df_long[grepl("_3_", names(df_long))]
    names(diferencia) <- c("antes","despues")
    diferencia$dif <- diferencia$despues-diferencia$antes
    diferencia %>% identify_outliers(dif) ##############################################################    
    normalidad <- diferencia %>% shapiro_test(dif)

    df_long <- pivot_longer(df_long, cols=names(df_long)[names(df_long)!="ID"], names_to = "Condition", values_to = "Value")    
    df_long[grepl("_bl_pl",df_long$Condition),]$Condition <- "Baseline" 
    df_long[grepl("_3_tmz",df_long$Condition),]$Condition <- "TMZ 3 Months"    
    df_long[grepl("_3_pl",df_long$Condition),]$Condition <- "Placebo 3 Months" 
    
    general_stat <- df_long %>% group_by(Condition) %>% get_summary_stats(Value, type = "mean_sd")
    g_s_3_tmz <- general_stat[general_stat$Condition == "TMZ 3 Months",]   
    g_s_3_tmz <- paste(signif(g_s_3_tmz[[4]],2),"±",signif(g_s_3_tmz[[5]],2), sep = "")
    g_s_3_pl <- general_stat[general_stat$Condition == "Placebo 3 Months",]  
    g_s_3_pl <- paste(signif(g_s_3_pl[[4]],2),"±",signif(g_s_3_pl[[5]],2), sep = "")
    g_s_bl_pl <- general_stat[general_stat$Condition == "Baseline",] 
    g_s_bl_pl <- paste(signif(g_s_bl_pl[[4]],2),"±",signif(g_s_bl_pl[[5]],2), sep = "")
   
    if (normalidad$p < 0.05){                
        before <- subset(df_long,  Condition == "Placebo 3 Months" , Value, drop = TRUE)
        # subset weight data after treatment
        after <- subset(df_long,  Condition == "TMZ 3 Months", Value, drop = TRUE)
        res <- wilcox.test(before, after, paired = TRUE)        
        return(list(g_s_bl_pl,g_s_3_pl,g_s_3_tmz,"Wilcoxon",res$p.value))} 
    df_long <- df_long[!grepl("Baseline",df_long$Condition) &  !grepl("_bl_tmz",df_long$Condition),]
    stat.test <- df_long  %>% t_test(Value ~ Condition, paired = TRUE, detailed = TRUE) %>% add_significance()
    magnitud_efecto <- df_long  %>% cohens_d(Value ~ Condition, paired = TRUE)
    if (stat.test$p <0.05) {    
    # Create a box plot
    bxp <- ggpaired(df_long, x = "Condition", y = "Value", order = c("Placebo 3 Months", "TMZ 3 Months"), ylab = "Value", xlab = "Condition")+ labs(title = gsub("_", " ", titulo))

    # Add p-value and significance levels
    stat.test <- stat.test %>% add_xy_position(x = "Condition")
    bxp <- bxp + stat_pvalue_manual(stat.test, tip.length = 0, linetype = 2) + labs(subtitle = get_test_label(stat.test, detailed= TRUE))
    
    ggsave(paste("T-tests/", titulo, ".png", sep = ""), plot = bxp, width = 9, height = 8.5, dpi = 300, units = "in")
    return(list(g_s_bl_pl,g_s_3_pl,g_s_3_tmz,"t-test",stat.test$p, magnitud_efecto))}
    return(list(g_s_bl_pl,g_s_3_pl,g_s_3_tmz, "t-test", stat.test$p, magnitud_efecto))}

info_det_tmz <- read.delim("info_determinaciones.tsv", sep = ";")
Tabla_estadistica <- NULL
for (det in determinaciones) {
    det_data <- tmz_cohort_data[grepl(paste("^", det, "_", sep = ""), names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]
    det_data <- det_data[grepl("_bl_pl",names(det_data)) | grepl("_3_pl",names(det_data))  | grepl("_bl_tmz",names(det_data))  | grepl("_3_tmz",names(det_data)) | grepl("ID", names(det_data))]
    analisis <- Paired_T_test_or_Wilcoxon_TMZ(det_data, det)
    det_unit <- info_det_tmz[grepl(paste(det,"_bl_pl", sep=""), info_det_tmz$Determinacion),]$Unidad
    det_name <- str_replace(info_det_tmz[grepl(paste(det,"_bl_pl", sep=""), info_det_tmz$Determinacion),]$Nombre,"Baseline ", "")
    if (det =="X6mwt") {
       det_name <- str_replace(info_det_tmz[grepl("6mwt_bl_pl", info_det_tmz$Determinacion),]$Nombre,"Baseline ", "")
       det_unit <- info_det_tmz[grepl("6mwt_bl_pl", info_det_tmz$Determinacion),]$Unidad}
    det_unit <- as.character(det_unit)
    print(det_unit)       
    Tabla_estadistica <- rbind(Tabla_estadistica, as.data.frame(cbind(det_name,analisis[[1]],analisis[[2]],analisis[[3]],det_unit, analisis[[4]], round(analisis[[5]],3))))    
    } 
names(Tabla_estadistica) <- c("Determinacion","Baseline","Placebo 3 Months", "TMZ 3 Months","Unit", "Test", "p-Value") 
Tabla_estadistica
library(kableExtra)
Tabla_estadistica %>%
  kbl(caption = "Statistical Analysis of the Trimetazidine Cohort") %>%
  kable_classic(full_width = F, html_font = "Cambria")

Tabla_estadistica2 <- NULL
for (det in determinaciones) {
    det_data <- tmz_cohort_data[grepl(paste("^", det, "_", sep = ""), names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]
    det_data <- det_data[grepl("_bl_pl",names(det_data)) | grepl("_3_pl",names(det_data))  | grepl("_bl_tmz",names(det_data))  | grepl("_3_tmz",names(det_data)) | grepl("ID", names(det_data))]
    det_data[grepl("_3_tmz", names(det_data))] <- det_data[grepl("_3_tmz", names(det_data))] - det_data[grepl("_bl_pl", names(det_data))]
    det_data[grepl("_3_pl", names(det_data))] <- det_data[grepl("_3_pl", names(det_data))] - det_data[grepl("_bl_pl", names(det_data))]
    det_data[grepl("_bl_tmz", names(det_data))] <- det_data[grepl("_bl_tmz", names(det_data))] - det_data[grepl("_bl_pl", names(det_data))]
    #det_data[grepl("_bl_pl", names(det_data))] <- 0
    

    analisis <- Paired_T_test_or_Wilcoxon_TMZ(det_data, det)
    det_unit <- info_det_tmz[grepl(paste(det,"_bl_pl", sep=""), info_det_tmz$Determinacion),]$Unidad
    det_name <- str_replace(info_det_tmz[grepl(paste(det,"_bl_pl", sep=""), info_det_tmz$Determinacion),]$Nombre,"Baseline ", "")
    if (det =="X6mwt") {
       det_name <- str_replace(info_det_tmz[grepl("6mwt_bl_pl", info_det_tmz$Determinacion),]$Nombre,"Baseline ", "")
       det_unit <- info_det_tmz[grepl("6mwt_bl_pl", info_det_tmz$Determinacion),]$Unidad}
    det_unit <- as.character(det_unit)
    print(det_unit)       
    Tabla_estadistica2 <- rbind(Tabla_estadistica2, as.data.frame(cbind(det_name,analisis[[1]],analisis[[2]],analisis[[3]],det_unit, analisis[[4]], round(analisis[[5]],3))))    
    } 
names(Tabla_estadistica2) <- c("Determinacion","Baseline","Placebo 3 Months", "TMZ 3 Months","Unit", "Test", "p-Value") 
Tabla_estadistica2
Tabla_estadistica2 %>%
  kbl(caption = "Statistical Analysis of the Trimetazidine Cohort") %>%
  kable_classic(full_width = F, html_font = "Cambria")

############################################################
# ML preparacion

tmz_cohort_data2 <- tmz_cohort_data[!grepl("_tmz",names(tmz_cohort_data)) & !grepl("_3_pl",names(tmz_cohort_data))] # Descartar mediciones de TMZ
names(tmz_cohort_data2)[!grepl("_bl_pl", names(tmz_cohort_data2))]
tmz_cohort_data2 <- tmz_cohort_data2[!grepl("_pl",names(tmz_cohort_data2)) | grepl("_bl_pl",names(tmz_cohort_data2))] # Descartar mediciones placebo remanentes
tmz_cohort_data2 <- tmz_cohort_data2[!grepl("estudio",names(tmz_cohort_data2))]
names(tmz_cohort_data2)[grepl("_bl", names(tmz_cohort_data2)) & !grepl("_bl_pl", names(tmz_cohort_data2))]
tmz_cohort_data2 <- tmz_cohort_data2[names(tmz_cohort_data2) != "psap_bl" & names(tmz_cohort_data2) != "papm_bl"] # descartar ya que se encuentra ademas una columna con la misma medicion basal en forma _bl_pl

diferencia_6mwt_TMZ_treat <- tmz_cohort_data$X6mwt_3_tmz - tmz_cohort_data$X6mwt_bl_pl
tmz_cohort_data2$TMZ_6mwt_treat_eff <- tmz_cohort_data$X6mwt_3_tmz - tmz_cohort_data2$X6mwt_bl_pl
write.table(tmz_cohort_data2, sep = "\t", file = "ML_TZM.tab", row.names = F, quote = F, col.names = T)













###########################################################

df_long <- tmz_cohort_data[grepl(paste("^", "X6mwt", "_", sep = ""), names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]

















###############################################
## Repeated Measures ANOVA #
Repeated_Measures_ANOVA <- function(df_long, titulo) {
    library(ggplot2)
    library(ggpubr)
    library(rstatix)
    library(tidyr)

    df_long <- df_long[order(df_long$ID),]
    df_long <- pivot_longer(df_long, cols=names(df_long)[names(df_long)!="ID"], names_to = "Condition", values_to = "Value")
    df_long[grepl("tmz",df_long$Condition) & grepl("bl",df_long$Condition),]$Condition <- "TMZ Wash-out"
    df_long[grepl("tmz",df_long$Condition) & grepl("3",df_long$Condition),]$Condition <- "TMZ 3 Months"
    df_long[grepl("pl",df_long$Condition) & grepl("bl",df_long$Condition),]$Condition <- "Baseline"
    df_long[grepl("pl",df_long$Condition) & grepl("3",df_long$Condition),]$Condition <- "Placebo 3 Months"

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
            print("no se puede continuar analisis")
            return()}}
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
    labs(subtitle = get_test_label(res.aov, detailed = FALSE), caption = get_pwc_label(pwc)) + labs(title = gsub("_", " ", titulo))
    ggsave(paste(titulo, ".png", sep = ""), plot = grafico, width = 9, height = 8.5, dpi = 300, units = "in")
    #system(paste("scp ", titulo, ".png adolfo@200.89.65.156:/run/media/vinicius/run-projects/Adolfo/Otros_proyectos/Cohorte_TMZ/", titulo, ".png", sep =""))
    return(list(res.aov,grafico))}

tapse_data <- tmz_cohort_data[grepl("tapse", names(tmz_cohort_data)) | grepl("ID", names(tmz_cohort_data))]
tapse_data_analisis <- Repeated_Measures_ANOVA(tapse_data, "TAPSE")
# Test de caminata
walk_data <- tmz_cohort_data[grepl("6mwt", names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]
walk_data_analysis <- Repeated_Measures_ANOVA(walk_data, "Walk_test")
# Cuantificacion de peptido natriuretico BNP
bnp_data <- tmz_cohort_data[grepl("bnp", names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]
bnp_data_analysis <- Repeated_Measures_ANOVA(bnp_data, "BNP levels")
ggplot(data = bnp_long, mapping = aes(x = Condition, y = Value, colour = Condition)) + geom_boxplot() + theme_bw() + theme(legend.position = "none")
friedman.test(bnp_long$Value,bnp_long$Condition, bnp_long$ID)
pairwise.wilcox.test(bnp_long$Value, bnp_long$Condition, paired = TRUE, p.adjust.method = "holm")
bnp_long <- bnp_long[bnp_long$ID != 4 & bnp_long$ID != 25,]
bnp_data_analysis <- Repeated_Measures_ANOVA(bnp_long, "BNP_levels")
# BORG 
borg_data <- tmz_cohort_data[grepl("borg", names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]
borg_data_analysis <- Repeated_Measures_ANOVA(borg_data, "Borg_dyspnea_index")

info_det_tmz <- read.delim("info_determinaciones.tsv", sep = ";")
Tabla_estadistica <- NULL
for (det in determinaciones) {
    det_data <- tmz_cohort_data[grepl(paste("^", det, "_", sep = ""), names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data))]
    det_data <- det_data[grepl("_bl_pl",names(det_data)) | grepl("_3_pl",names(det_data))  | grepl("_bl_tmz",names(det_data))  | grepl("_3_tmz",names(det_data)) | grepl("ID", names(det_data))]
    analisis <- Repeated_Measures_ANOVA(det_data, det)
    det_unit <- info_det_tmz[grepl(paste(det,"_bl_pl", sep=""), info_det_tmz$Determinacion),]$Unidad
    det_name <- str_replace(info_det_tmz[grepl(paste(det,"_bl_pl", sep=""), info_det_tmz$Determinacion),]$Nombre,"Baseline ", "")
    if (det =="X6mwt") {
       det_name <- str_replace(info_det_tmz[grepl("6mwt_bl_pl", info_det_tmz$Determinacion),]$Nombre,"Baseline ", "")
       det_unit <- info_det_tmz[grepl("6mwt_bl_pl", info_det_tmz$Determinacion),]$Unidad}
    det_unit <- as.character(det_unit)
    print(det_unit)       
    #Tabla_estadistica <- rbind(Tabla_estadistica, as.data.frame(cbind(det_name,analisis[[1]],analisis[[2]],analisis[[3]],det_unit, analisis[[4]], round(analisis[[5]],3))))    
    } 
names(Tabla_estadistica) <- c("Determinacion","Baseline","Placebo 3 Months", "TMZ 3 Months","Unit", "Test", "p-Value") 
Tabla_estadistica




###############################################
## linear mixed model
library(nlme)
library(lme4)
#tapse_data_mlm <- tmz_cohort_data[grepl("tapse", names(tmz_cohort_data)) | grepl("ID", names(tmz_cohort_data)) | grepl("droga.estudio.1", names(tmz_cohort_data))]
df_long <- tmz_cohort_data[grepl("X6mwt", names(tmz_cohort_data))  | grepl("ID", names(tmz_cohort_data)) | grepl("droga.estudio.1", names(tmz_cohort_data))]
df_long <- df_long[order(df_long$ID),]
df_long <- df_long[,!grepl("_bl_tmz", names(df_long))]
#df_long <- df_long[,!grepl("_bl_pl", names(df_long))]

df_long <- pivot_longer(df_long, cols=names(df_long)[grepl("_3_tmz", names(df_long)) | grepl("_3_pl", names(df_long))], names_to = "Condition", values_to = "Value")
df_long$treatment <- "a"
df_long$period <- "b"
#df_long[grepl("tmz",df_long$Condition) & grepl("3",df_long$Condition),]$Condition <- "TMZ 3 Months"    
#df_long[grepl("pl",df_long$Condition) & grepl("3",df_long$Condition),]$Condition <- "Placebo 3 Months" 
    df_long[grepl("_3_tmz",df_long$Condition),]$Condition <- "TMZ 3 Months"    
    df_long[grepl("_3_pl",df_long$Condition),]$Condition <- "Placebo 3 Months"
    names(df_long)[grepl("_bl_pl",names(df_long))] <- "Baseline"  
df_long[df_long$droga.estudio.1 == "z" & df_long$Condition == "TMZ 3 Months",]$treatment <- "1"
df_long[df_long$droga.estudio.1 == "z" & df_long$Condition == "TMZ 3 Months",]$period <- "0"
df_long[df_long$droga.estudio.1 == "y" & df_long$Condition == "TMZ 3 Months",]$period <- "1"
df_long[df_long$droga.estudio.1 == "y" & df_long$Condition == "TMZ 3 Months",]$treatment <- "1"

df_long[df_long$droga.estudio.1 == "z" & df_long$Condition == "Placebo 3 Months",]$treatment <- "0"
df_long[df_long$droga.estudio.1 == "z" & df_long$Condition == "Placebo 3 Months",]$period <- "1"
df_long[df_long$droga.estudio.1 == "y" & df_long$Condition == "Placebo 3 Months",]$period <- "0"
df_long[df_long$droga.estudio.1 == "y" & df_long$Condition == "Placebo 3 Months",]$treatment <- "0"
df_long$treatment <- as.numeric(as.character(df_long$treatment))
df_long$period <- as.numeric(as.character(df_long$period))
mixed_model <-lme(Value~Condition+period, random= ~1|ID,data=df_long) 
summary(mixed_model)
pred <- ggpredict(mixed_model, c("Baseline", "Condition"))
pred <- ggpredict(mixed_model, "Condition")
plot(pred)
################################################
