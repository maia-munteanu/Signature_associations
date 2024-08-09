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
cna <- as.data.frame(fread(args[3]))
metadata <- as.data.frame(fread(args[4])); metadata$LizaCancerType=gsub("_MSI","",metadata$LizaCancerType)
model_type <- as.character(args[5])
covariates <- as.logical(args[6])
output <- as.character(args[7])

exposures<-exposures[,c("sample",signature)] 
cna %>% pivot_longer(cols=-c(chromosome,start,end,gene),names_to="sample",values_to="Freq") %>% dplyr::select(sample, gene, Freq) %>% 
  dplyr::filter(sample %in% exposures$sample) %>% left_join(exposures) %>% left_join(metadata[,c("sample","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType")]) -> cna
colnames(cna)<-c("sample","gene","CN","Exposures","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType")

if (covariates){
  if (model_type=="GLMnb"){
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
      for (gene in unique(cna$gene)){
          df<-cna[which(cna$gene == gene),]
          df$CN <- ifelse(df$CN > 0, df$CN, 0)
          df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
          model <- glm.nb(Exposures ~ CN + LizaCancerType + gender + purity + ploidy + tmbStatus + msStatus, data = df)
          beta <- coef(model)["CN"]
          se <- summary(model)$coefficients["CN", "Std. Error"]
          p_value <- summary(model)$coefficients["CN", "Pr(>|z|)"]
          results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      write.table(results, file = paste0(signature, "_amp.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
      
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
      for (gene in unique(cna$gene)){
        df<-cna[which(cna$gene == gene),]
        df$CN <- ifelse(df$CN < 0, df$CN, 0)
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
        model <- glm.nb(Exposures ~ CN + LizaCancerType + gender + purity + ploidy + tmbStatus + msStatus, data = df)
        beta <- coef(model)["CN"]
        se <- summary(model)$coefficients["CN", "Std. Error"]
        p_value <- summary(model)$coefficients["CN", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      write.table(results, file = paste0(signature, "_del.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  
  if (model_type=="beta"){
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
      for (gene in unique(cna$gene)){
          df<-cna[which(cna$gene == gene),]
          df$CN <- ifelse(df$CN > 0, df$CN, 0)
          df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
          df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
          df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
          model <- betareg(Exposures ~ CN + LizaCancerType + gender + purity + ploidy + tmbStatus + msStatus, data = df)
          beta <- coef(model)["CN"]
          se <- summary(model)$coefficients$mean["CN", "Std. Error"]
          p_value <- summary(model)$coefficients$mean["CN", "Pr(>|z|)"]
          results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      write.table(results, file = paste0(signature, "_amp.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
      
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
      for (gene in unique(cna$gene)){
        df<-cna[which(cna$gene == gene),]
        df$CN <- ifelse(df$CN < 0, df$CN, 0)
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
        df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
        df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
        model <- betareg(Exposures ~ CN + LizaCancerType + gender + purity + ploidy + tmbStatus + msStatus, data = df)
        beta <- coef(model)["CN"]
        se <- summary(model)$coefficients$mean["CN", "Std. Error"]
        p_value <- summary(model)$coefficients$mean["CN", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      write.table(results, file = paste0(signature, "_del.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
}else{
  if (model_type=="GLMnb"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN > 0, df$CN, 0)
      df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
      model <- glm.nb(Exposures ~ CN, data = df)
      beta <- coef(model)["CN"]
      se <- summary(model)$coefficients["CN", "Std. Error"]
      p_value <- summary(model)$coefficients["CN", "Pr(>|z|)"]
      results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_amp.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN < 0, df$CN, 0)
      df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
      model <- glm.nb(Exposures ~ CN , data = df)
      beta <- coef(model)["CN"]
      se <- summary(model)$coefficients["CN", "Std. Error"]
      p_value <- summary(model)$coefficients["CN", "Pr(>|z|)"]
      results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_del.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  
  if (model_type=="beta"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN > 0, df$CN, 0)
      df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
      df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
      df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
      model <- betareg(Exposures ~ CN , data = df)
      beta <- coef(model)["CN"]
      se <- summary(model)$coefficients$mean["CN", "Std. Error"]
      p_value <- summary(model)$coefficients$mean["CN", "Pr(>|z|)"]
      results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_amp.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(cna$gene)){
      df<-cna[which(cna$gene == gene),]
      df$CN <- ifelse(df$CN < 0, df$CN, 0)
      df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"
      df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
      df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
      model <- betareg(Exposures ~ CN, data = df)
      beta <- coef(model)["CN"]
      se <- summary(model)$coefficients$mean["CN", "Std. Error"]
      p_value <- summary(model)$coefficients$mean["CN", "Pr(>|z|)"]
      results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, "_del.tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
}
