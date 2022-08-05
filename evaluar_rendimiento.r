#!/usr/bin/env Rscript
R
evaluar_rendimientos <- function(){
    info_snp2 <- read.delim("info_snp_tab_INT_model.tsv", sep = "\t")
    info_snp2 <- info_snp2[c("rsid","Score_fix","Score2_fix")]
    info_snp2 <- merge(info_snp2,info_snp[info_snp$best_grid_sp_mean_K != 0,],by="rsid")
    info_snp2$Score_fix <- info_snp2$Score_fix * info_snp2$best_grid_sp_mean_K
    info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K
    print(dim(info_snp2))
    pred_basal_pob1 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = sub_pob1, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_pob1 <- AUCBoot(pred_basal_pob1, y[sub_pob1], seed = 1)
    pred_basal_pob2 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = sub_pob2, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_pob2 <- AUCBoot(pred_basal_pob2, y[sub_pob2], seed = 1)
    pred_basal_pob3 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = sub_pob3, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_pob3 <- AUCBoot(pred_basal_pob3, y[sub_pob3], seed = 1)
    
    pred_mod2_pob1 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = sub_pob1, ind.col = info_snp2$"_NUM_ID_")    
    AUC_pred_mod2_pob1 <- AUCBoot(pred_mod2_pob1, y[sub_pob1], seed = 1)
    pred_mod2_pob2 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = sub_pob2, ind.col = info_snp2$"_NUM_ID_")    
    AUC_pred_mod2_pob2 <- AUCBoot(pred_mod2_pob2, y[sub_pob2], seed = 1)
    pred_mod2_pob3<- big_prodVec(G, info_snp2$Score2_fix, ind.row = sub_pob3, ind.col = info_snp2$"_NUM_ID_")    
    AUC_pred_mod2_pob3 <- AUCBoot(pred_mod2_pob3, y[sub_pob3], seed = 1)
    rendimientos <- as.data.frame(rbind(AUC_pred_baseline_pob1,AUC_pred_baseline_pob2,AUC_pred_baseline_pob3,AUC_pred_mod2_pob1,AUC_pred_mod2_pob2,AUC_pred_mod2_pob3))
    print(rendimientos)
    print((AUC_pred_baseline_pob1[1]+AUC_pred_baseline_pob2[1]+AUC_pred_baseline_pob3[1])/3)
    print((AUC_pred_mod2_pob1[1]+AUC_pred_mod2_pob2[1]+AUC_pred_mod2_pob3[1])/3)}

evaluar_rendimientos_DE <- function(){
    info_snp2 <- read.delim("info_snp_tab_DE_model.tsv", sep = "\t")
    info_snp2 <- info_snp2[c("rsid","Score_fix","Score2_fix")]
    info_snp2 <- merge(info_snp2,info_snp[info_snp$best_grid_sp_mean_K != 0,],by="rsid")
    info_snp2$Score_fix <- info_snp2$Score_fix * info_snp2$best_grid_sp_mean_K
    info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K
    print(dim(info_snp2))
    pred_basal_pob1 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = sub_pob1, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_pob1 <- AUCBoot(pred_basal_pob1, y[sub_pob1], seed = 1)
    pred_basal_pob2 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = sub_pob2, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_pob2 <- AUCBoot(pred_basal_pob2, y[sub_pob2], seed = 1)
    pred_basal_pob3 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = sub_pob3, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_pob3 <- AUCBoot(pred_basal_pob3, y[sub_pob3], seed = 1)
    
    pred_mod2_pob1 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = sub_pob1, ind.col = info_snp2$"_NUM_ID_")    
    AUC_pred_mod2_pob1 <- AUCBoot(pred_mod2_pob1, y[sub_pob1], seed = 1)
    pred_mod2_pob2 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = sub_pob2, ind.col = info_snp2$"_NUM_ID_")    
    AUC_pred_mod2_pob2 <- AUCBoot(pred_mod2_pob2, y[sub_pob2], seed = 1)
    pred_mod2_pob3<- big_prodVec(G, info_snp2$Score2_fix, ind.row = sub_pob3, ind.col = info_snp2$"_NUM_ID_")    
    AUC_pred_mod2_pob3 <- AUCBoot(pred_mod2_pob3, y[sub_pob3], seed = 1)
    rendimientos <- as.data.frame(rbind(AUC_pred_baseline_pob1,AUC_pred_baseline_pob2,AUC_pred_baseline_pob3,AUC_pred_mod2_pob1,AUC_pred_mod2_pob2,AUC_pred_mod2_pob3))
    print(rendimientos)
    print((AUC_pred_baseline_pob1[1]+AUC_pred_baseline_pob2[1]+AUC_pred_baseline_pob3[1])/3)
    print((AUC_pred_mod2_pob1[1]+AUC_pred_mod2_pob2[1]+AUC_pred_mod2_pob3[1])/3)}

comparar_roc_curves <- function(){
    library(pROC)
    pred_basal_758_2 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = ind.final.val, ind.col = info_snp$'_NUM_ID_')
    AUC_pred_baseline_758_2 <- AUCBoot(pred_basal_758_2, y[ind.final.val], seed = 1)
    info_snp2 <- read.delim("info_snp_tab_INT_model.tsv", sep = "\t")
    info_snp2 <- info_snp2[c("rsid","Score_fix","Score2_fix")]
    info_snp2 <- merge(info_snp2,info_snp[info_snp$best_grid_sp_mean_K != 0,],by="rsid")
    info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K
    pred_INT_model <- big_prodVec(G, info_snp2$Score2_fix, ind.row = ind.final.val, ind.col = info_snp2$"_NUM_ID_")
    AUC_INT_model <- AUCBoot(pred_INT_model, y[ind.final.val], seed = 1)
    
    info_snp2 <- read.delim("info_snp_tab_DE_model.tsv", sep = "\t")
    info_snp2 <- info_snp2[c("rsid","Score_fix","Score2_fix")]
    info_snp2 <- merge(info_snp2,info_snp[info_snp$best_grid_sp_mean_K != 0,],by="rsid")
    info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K
    pred_DE_model <- big_prodVec(G, info_snp2$Score2_fix, ind.row = ind.final.val, ind.col = info_snp2$"_NUM_ID_")
    AUC_DE_model <- AUCBoot(pred_DE_model, y[ind.final.val], seed = 1)
    print(AUC_pred_baseline_758_2)
    print(AUC_INT_model)
    print(AUC_DE_model)
    roc1 <- roc(y[ind.final.val], pred_INT_model)
    roc2 <- roc(y[ind.final.val], pred_basal_758_2)
    roc.test(roc1, roc2)
    #tiff("roc_comparacion_modelos_generados.tiff",width = 1080, height = 1080, pointsize = 27)
    rocobj1 <- plot.roc(y[ind.final.val], pred_INT_model,
                    main="Statistical comparison",
                    percent=TRUE,
                    col="#1c61b6")
    rocobj2 <- lines.roc(y[ind.final.val], pred_basal_758_2, 
                     percent=TRUE, 
                     col="#008600")
    testobj <- roc.test(rocobj1, rocobj2, method="bootstrap")
    print(testobj)
    text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
    legend("bottomright", legend=c("INT Model", "PRS Basal"), col=c("#1c61b6", "#008600"), lwd=2)
    #dev.off()
    roc1 <- roc(y[ind.final.val], pred_DE_model)
    roc2 <- roc(y[ind.final.val], pred_basal_758_2)
    roc.test(roc1, roc2)
    #tiff("roc_comparacion_modelos_generados2.tiff",width = 1080, height = 1080, pointsize = 27)
    rocobj1 <- plot.roc(y[ind.final.val], pred_DE_model,
                    main="Statistical comparison",
                    percent=TRUE,
                    col="#1c61b6")
    rocobj2 <- lines.roc(y[ind.final.val], pred_basal_758_2, 
                     percent=TRUE, 
                     col="#008600")
    testobj <- roc.test(rocobj1, rocobj2, method="bootstrap")
    print(testobj)
    text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
    legend("bottomright", legend=c("DE Model", "PRS Basal"), col=c("#1c61b6", "#008600"), lwd=2)
    #dev.off()
    }

library(bigsnpr)
load("CIMBA_GRID-SP.RData")
jk <- y
load("../1er_Objetivo/2_entorno_validacion_cruzada.RData")
pred_basal_758 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = ind.final.val, ind.col = info_snp$'_NUM_ID_')
pred_basal_1500 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = poblacion_restante, ind.col = info_snp$'_NUM_ID_')
pred_basal_CIMBA <- big_prodVec(imputed, info_CIMBA$best_grid_sp_mean_K, ind.col = info_CIMBA$'_NUM_ID_')
AUC_pred_baseline_1500 <- AUCBoot(pred_basal_1500, y[poblacion_restante], seed = 1) #0.60762861
AUC_pred_baseline_758 <- AUCBoot(pred_basal_758, y[ind.final.val], seed = 1) #0.60928497
AUC_pred_baseline_CIMBA  <- AUCBoot(pred_basal_CIMBA, jk, seed = 1)

set.seed(1)
#ind.final.val <- sample(nrow(G), 758)
sub_pob1 <- sample(ind.final.val, 258)
sub_pob2 <- sample(setdiff(ind.final.val, sub_pob1), 250)
sub_pob3 <- setdiff(setdiff(ind.final.val, sub_pob1),sub_pob2)

evaluar_rendimientos()
evaluar_rendimientos_DE()

evaluar_rendimientos2 <- function(){
    info_snp2 <- read.delim("info_snp_tab_INT_model.tsv", sep = "\t")
    print(dim(info_snp2))
    info_snp2 <- info_snp2[c("rsid","Score_fix","Score2_fix")]
    info_snp_CIMBA <- merge(info_snp2,info_CIMBA,by="rsid")
    info_snp_CIMBA$Score_fix <- info_snp_CIMBA$Score_fix * info_snp_CIMBA$best_grid_sp_mean_K
    info_snp_CIMBA$Score2_fix <- info_snp_CIMBA$Score2_fix * info_snp_CIMBA$best_grid_sp_mean_K
    info_snp2 <- merge(info_snp2,info_snp[info_snp$best_grid_sp_mean_K != 0,],by="rsid")
    info_snp2$Score_fix <- info_snp2$Score_fix * info_snp2$best_grid_sp_mean_K
    info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K   
    print(dim(info_snp2))
    print(dim(info_snp_CIMBA))
    
    pred_basal_758 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = ind.final.val, ind.col = info_snp$'_NUM_ID_')
    pred_basal_1500 <- big_prodVec(G, info_snp$best_grid_sp_mean_K, ind.row = poblacion_restante, ind.col = info_snp$'_NUM_ID_')
    pred_basal_CIMBA <- big_prodVec(imputed, info_snp_CIMBA$best_grid_sp_mean_K, ind.col = info_snp_CIMBA$'_NUM_ID_', )
    pred_mod2_1500 <- big_prodVec(G, info_snp2$Score_fix, ind.row = poblacion_restante, ind.col = info_snp2$"_NUM_ID_")
    pred_mod2_758 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = ind.final.val, ind.col = info_snp2$"_NUM_ID_")    
    pred_CIMBA2 <- big_prodVec(imputed, info_snp_CIMBA$Score2_fix, ind.col = info_snp_CIMBA$"_NUM_ID_")

    AUC_pred_baseline_1500 <- AUCBoot(pred_basal_1500, y[poblacion_restante], seed = 1) #0.60762861
    AUC_pred_baseline_758 <- AUCBoot(pred_basal_758, y[ind.final.val], seed = 1) #0.60928497
    AUC_pred_baseline_CIMBA  <- AUCBoot(pred_basal_CIMBA, jk, seed = 1) #0.5362011
    AUC_pred_mod2_1500 <- AUCBoot(pred_mod2_1500, y[poblacion_restante], seed = 1)
    AUC_pred_mod2_758 <- AUCBoot(pred_mod2_758, y[ind.final.val], seed = 1)
    AUC_pred_CIMBA2 <- AUCBoot(pred_CIMBA2, jk, seed = 1)
    rendimientos <- as.data.frame(rbind(AUC_pred_baseline_758,AUC_pred_baseline_1500,AUC_pred_mod2_1500,AUC_pred_mod2_758,AUC_pred_baseline_CIMBA,AUC_pred_CIMBA2))
    print(rendimientos)}

evaluar_rendimientos2()

pred_mod1_1500 <- big_prodVec(G, info_snp2$Score_fix, ind.row = poblacion_restante, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)
pred_mod2_1500 <- big_prodVec(G, info_snp2$Score2_fix, ind.row = poblacion_restante, ind.col = info_snp2$X_NUM_ID_, ncores = NCORES)
AUC_pred_mod1_1500 <- AUCBoot(pred_mod1_1500, y[poblacion_restante], seed = 1)
AUC_pred_mod2_1500 <- AUCBoot(pred_mod2_1500, y[poblacion_restante], seed = 1)

library(pROC)
roc1 <- roc(y[poblacion_restante], pred_sp)
roc2 <- roc(y[poblacion_restante], pred_nosp)
roc.test(roc1, roc2)

tiff("roc_comparacion.tiff",width = 1080, height = 1080, pointsize = 27)
rocobj1 <- plot.roc(y[poblacion_restante], pred_nosp,
                    main="Statistical comparison",
                    percent=TRUE,
                    col="#1c61b6")
rocobj2 <- lines.roc(y[poblacion_restante], pred_sp, 
                     percent=TRUE, 
                     col="#008600")
testobj <- roc.test(rocobj1, rocobj2)
text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
legend("bottomright", legend=c("GRID-NOSP", "GRID-SP"), col=c("#1c61b6", "#008600"), lwd=2)
dev.off()

system("scp -P 1313 roc_comparacion.tiff adolfo@200.89.65.156:/media/run-projects/Adolfo/Resultados_Tesis/Objetivo_1")

comparar_roc_curves()

info_snp2 <- read.delim("info_snp_tab_DE_model.tsv", sep = "\t")
#info_snp2 <- read.delim("info_snp_tab_INT_model.tsv", sep = "\t")
info_snp2 <- info_snp2[c("rsid","Score_fix","Score2_fix")]
info_snp2 <- merge(info_snp2,info_snp[info_snp$best_grid_sp_mean_K != 0,],by="rsid")
info_snp2$Score2_fix <- info_snp2$Score2_fix * info_snp2$best_grid_sp_mean_K
pred_INT_model <- big_prodVec(G, info_snp2$Score2_fix, ind.row = ind.final.val, ind.col = info_snp2$"_NUM_ID_")
Prediccion <- as.data.frame(cbind(obj.bigSNP$fam$sample.ID[ind.final.val],y[ind.final.val],pred_INT_model))
Prediccion$V2 <- as.character(Prediccion$V2)
Prediccion[Prediccion$V2 == 0,]$V2 <- "Control"
Prediccion[Prediccion$V2 == 1,]$V2 <- "Case"
Prediccion$V2 <- as.factor(Prediccion$V2)
Prediccion$pred_INT_model <- as.character(Prediccion$pred_INT_model)
Prediccion$pred_INT_model <- as.numeric(Prediccion$pred_INT_model)
less_risk <- head(Prediccion[order(Prediccion$pred_INT_model),],length(ind.final.val)*0.95)
low_risk <- head(Prediccion[order(Prediccion$pred_INT_model),],length(ind.final.val)*0.05)
odds_less_risk <- length(less_risk[less_risk$V2 == "Case",]$V2)/length(less_risk[less_risk$V2 == "Control",]$V2)
greater_risk <- tail(Prediccion[order(Prediccion$pred_INT_model),],length(ind.final.val)*0.05)
odds_greater_risk <- length(greater_risk[greater_risk$V2 == "Case",]$V2)/length(greater_risk[greater_risk$V2 == "Control",]$V2)
odds_ratio <- odds_greater_risk/odds_less_risk
library(ggpubr)
ggdensity(Prediccion, x = "pred_INT_model", add = "mean", rug = TRUE, color = "V2", fill = "V2", palette = c("#00AFBB", "#E7B800")) +
    ggtitle("Distribution of test samples PRS", subtitle = paste(best_auc_name, " = ", round(best_auc_value,4), sep = ""))+ 
    theme(legend.position = c(0.85,0.75), plot.title = element_text(hjust = 0, size=20, vjust= 0), legend.title = element_text(size = 12), legend.text = element_text(size = 9), legend.direction = "vertical") +
    ylab("Density")+ xlab("PRS") + labs(fill = "Sample type", color ="Sample type")#+ geom_vline(xintercept = greater_risk$pred_INT_model[1], na.rm = FALSE, show.legend = NA)

ggsave(paste("Density_plot_",best_auc_name2,"_p-Value_threshold_", p_value,".png", sep = ""), dpi = 600, units = "in") # width = 16, height = 9,