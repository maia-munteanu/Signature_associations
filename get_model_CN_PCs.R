library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(pscl)
library(betareg)

args=commandArgs(TRUE)
signature <- as.character(args[1])
exposures <- as.data.frame(fread(args[2]))
cna_pcs <- read.table(args[3],row.names = 1); cna_pcs$sample=rownames(cna_pcs)
metadata <- as.data.frame(fread(args[4])); metadata$LizaCancerType=gsub("_MSI","",metadata$LizaCancerType)
model_type <- as.character(args[5])
covariates <- as.logical(args[6])

exposures<-exposures[,c("sample",signature)] 
cna_pcs %>% pivot_longer(cols=c(-sample)) %>% left_join(metadata[,c("sample","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType")]) %>% dplyr::filter(sample %in% exposures$sample)  %>% left_join(exposures) -> cna_pcs
colnames(cna_pcs)<-c("sample","PC","PC_value","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType","Exposures")

if (covariates){
  if (model_type=="GLMnb"){
      results=data.frame(Signature = c(), PC = c(), Beta = c(), SE = c(), P_Value = c())
      for (PC in unique(cna_pcs$PC)){
        df<-cna_pcs[which(cna_pcs$PC == PC),]
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
        model <- glm.nb(Exposures ~ PC_value + LizaCancerType + gender + purity + tmbStatus + msStatus, data = df)
        beta <- coef(model)["PC_value"]
        se <- summary(model)$coefficients["PC_value", "Std. Error"]
        p_value <- summary(model)$coefficients["PC_value", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, PC = PC, Beta = beta, SE = se, P_Value = p_value))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  
  if (model_type=="beta"){
      results=data.frame(Signature = c(), PC = c(), Beta = c(), SE = c(), P_Value = c())
      for (PC in unique(cna_pcs$PC)){
        df<-cna_pcs[which(cna_pcs$PC == PC),]
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
        df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
        df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
        model <- betareg(Exposures ~ PC_value + LizaCancerType + gender + purity + tmbStatus + msStatus, data = df)
        beta <- coef(model)["PC_value"]
        se <- summary(model)$coefficients$mean["PC_value", "Std. Error"]
        p_value <- summary(model)$coefficients$mean["PC_value", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, PC = PC, Beta = beta, SE = se, P_Value = p_value))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
}else{
  if (model_type=="GLMnb"){
    results=data.frame(Signature = c(), PC = c(), Beta = c(), SE = c(), P_Value = c())
    for (PC in unique(cna_pcs$PC)){
      df<-cna_pcs[which(cna_pcs$PC == PC),]
      df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
      model <- glm.nb(Exposures ~ PC_value , data = df)
      beta <- coef(model)["PC_value"]
      se <- summary(model)$coefficients["PC_value", "Std. Error"]
      p_value <- summary(model)$coefficients["PC_value", "Pr(>|z|)"]
      results<-rbind(results,data.frame(Signature = signature, PC = PC, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  
  if (model_type=="beta"){
    results=data.frame(Signature = c(), PC = c(), Beta = c(), SE = c(), P_Value = c())
    for (PC in unique(cna_pcs$PC)){
      df<-cna_pcs[which(cna_pcs$PC == PC),]
      df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
      df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
      df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
      model <- betareg(Exposures ~ PC_value, data = df)
      beta <- coef(model)["PC_value"]
      se <- summary(model)$coefficients$mean["PC_value", "Std. Error"]
      p_value <- summary(model)$coefficients$mean["PC_value", "Pr(>|z|)"]
      results<-rbind(results,data.frame(Signature = signature, PC = PC, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
}

