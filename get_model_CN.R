library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(pscl)
library(betareg)

args=commandArgs(TRUE)
args=c("SBS96A_InClu","Indel_signatures_raw.tsv","CNA_by_gene_all_samples_restrained.tsv","/g/strcombio/fsupek_cancer3/SV_clusters_project/Pipeline_inputs/Hartwig_PCAWG_TCGA_MMRF_CPTAC_OVCARE_MUTes.tsv","tes","test")

signature <- as.character(args[1])
exposures <- as.data.frame(fread(args[2]))
cna <- as.data.frame(fread(args[3]))
metadata <- as.data.frame(fread(args[4])); metadata$LizaCancerType=gsub("_MSI","",metadata$LizaCancerType)
model_type <- as.character(args[5])
covariates <- as.logical(args[6])

exposures<-exposures[,c("sample",signature)] 
cna %>% pivot_longer(cols=-c(chromosome,start,end,gene),names_to="sample",values_to="Freq") %>% dplyr::select(sample, gene, Freq) %>% 
  dplyr::filter(sample %in% exposures$sample) %>% left_join(exposures) %>% left_join(metadata[,c("sample","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType","primaryTumorLocation")]) -> cna
colnames(cna)<-c("sample","gene","CN","Exposures","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType","primaryTumorLocation")
  
if (model_type=="GLMglog2"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN > 0, df$CN, 0)
      if(nrow(df[which(df$CN > 0),])>=500){
          df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
          df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"
          model <- glm(log2(Exposures + 1) ~ CN + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)
          beta <- coef(model)["CN"]
          se <- summary(model)$coefficients["CN", "Std. Error"]
          p_value <- summary(model)$coefficients["CN", "Pr(>|t|)"]
          results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_amp.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN < 0, df$CN, 0)
      if(nrow(df[which(df$CN < 0),])>=500){
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"
        model <- glm(log2(Exposures + 1) ~ CN + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)
        beta <- coef(model)["CN"]
        se <- summary(model)$coefficients["CN", "Std. Error"]
        p_value <- summary(model)$coefficients["CN", "Pr(>|t|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_del.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
  
if (model_type=="beta"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN > 0, df$CN, 0)
      if(nrow(df[which(df$CN > 0),])>=500){
          df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
          df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"
          df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
          df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
          model <- betareg(Exposures ~ CN + LizaCancerType + gender + purity + tmbStatus + msStatus, data = df)
          beta <- coef(model)["CN"]
          se <- summary(model)$coefficients$mean["CN", "Std. Error"]
          p_value <- summary(model)$coefficients$mean["CN", "Pr(>|z|)"]
          results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_amp.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN < 0, df$CN, 0)
      if(nrow(df[which(df$CN < 0),])>=500){
          df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
          df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"
          df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
          df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
          model <- betareg(Exposures ~ CN + LizaCancerType + gender + purity + tmbStatus + msStatus, data = df)
          beta <- coef(model)["CN"]
          se <- summary(model)$coefficients$mean["CN", "Std. Error"]
          p_value <- summary(model)$coefficients$mean["CN", "Pr(>|z|)"]
          results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_del.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
