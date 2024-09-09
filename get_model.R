library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(pscl)
library(betareg)
library(arm)
library(statmod)
library(mgcv)
library(tweedie)
#library(tweedie, lib.loc="/g/strcombio/fsupek_home/mmunteanu/.conda/envs/SigProfilerAssignment/lib/R/library/")

args=commandArgs(TRUE)
signature <- as.character(args[1])
input <- as.character(args[2])
model_type <- as.character(args[3])
metadata <- as.data.frame(fread(args[4])); metadata$LizaCancerType=gsub("_MSI","",metadata$LizaCancerType)
covariates <- as.logical(args[5])

germline=as.data.frame(fread(input))
germline %>% pivot_longer(cols = -c(sample, Gene.refGene, Freq, Clustered, Unclustered), names_to = "Signature", values_to = "Exposures") %>% filter(Signature==signature) -> germline
germline %>% left_join(metadata[,c("sample","gender","purity","ploidy","msStatus","tmbStatus","primaryTumorLocation","LizaCancerType")]) -> germline

# bd.xi <- tweedie.profile(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, data=df,do.plot=TRUE, xi.vec=seq(1.1,1.9,by=0.1))
# bd.xi$xi.max
# bd.m <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, data=df,family=tweedie(link.power=0, var.power=bd.xi$xi.max))
# qqnorm( resid(bd.m), las=1 ); qqline( resid(bd.m) )

if (covariates) {
    if (model_type=="GLMglog2") {
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
            df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
            model <- glm(log2(Exposures + 1) ~ Mutation_Score + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="GLMglog2_logSum") {
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
            df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
            if (grepl("Clu",signature)){model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)} 
            else {model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)}
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="bGLMglog2_logSum") {
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
            df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
            if (grepl("Clu",signature)){model <- bayesglm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)} 
            else {model <- bayesglm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, family = gaussian(), data = df)}
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
     if (model_type=="Tweedielog2_logSum") {
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
            df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
            if (grepl("Clu",signature)){model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = 1.5, link.power = 0),  data = df)} 
            else {model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = 1.5, link.power = 0),  data = df)}
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
     } 
    if (model_type=="pTweedie_logoSum_pvar") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        if (grepl("Clu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = p, link.power = 0),  data = df)} 
        if (grepl("Uclu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = p, link.power = 0),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pTweedie_logoSum_pvar_inter") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        if (grepl("Clu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = p, link.power = 0),  data = df)} 
        if (grepl("Uclu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = p, link.power = 0),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }   
    if (model_type=="pTweedie_logoSum_pvar_2") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        if (grepl("Clu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + purity, family = tweedie(var.power = p, link.power = 0),  data = df)} 
        if (grepl("Uclu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation + msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation + msStatus + purity,family = tweedie(var.power = p, link.power = 0),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pTweedie_logoSum_pvar_2_inter") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        if (grepl("Clu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation*msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation*msStatus + purity, family = tweedie(var.power = p, link.power = 0),  data = df)} 
        if (grepl("Uclu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation*msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation*msStatus + purity,family = tweedie(var.power = p, link.power = 0),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pTweedie_logoSum_pvar_3") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        if (grepl("Clu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Clustered)) + LizaCancerType + msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + LizaCancerType + msStatus + purity, family = tweedie(var.power = p, link.power = 0),  data = df)} 
        if (grepl("Uclu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + LizaCancerType + msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + LizaCancerType + msStatus + purity,family = tweedie(var.power = p, link.power = 0),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pTweedie_logoSum_pvar_3_inter") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        if (grepl("Clu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Clustered)) + LizaCancerType*msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + LizaCancerType*msStatus + purity, family = tweedie(var.power = p, link.power = 0),  data = df)} 
        if (grepl("Uclu",signature)){
          p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + LizaCancerType*msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + LizaCancerType*msStatus + purity,family = tweedie(var.power = p, link.power = 0),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pTweedie_logoSumUnclustered_pvar_3_inter") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Power = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType); df$LizaCancerType = relevel(df$LizaCancerType, ref = "Other")
        p=as.numeric(sub(".*p=([0-9.]+).*", "\\1",gam(Exposures ~ Mutation_Score + offset(log(Unclustered)) + LizaCancerType*msStatus + purity, data = df, family = tw(link="log"), method = "ML")$family$family))
        model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + LizaCancerType*msStatus + purity,family = tweedie(var.power = p, link.power = 0),  data = df)
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        AIC <- AICtweedie(model, dispersion=NULL, k = 2, verbose=TRUE)
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Power = p, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pGLM_logoSum") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType)
        if (grepl("Clu",signature)){
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)} 
        if (grepl("Uclu",signature)){
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC(model)))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="pGLM_logoSum_inter") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType)
        if (grepl("Clu",signature)){
          model <- glm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)} 
        if (grepl("Uclu",signature)){
          model <- glm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC(model)))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
   if (model_type=="bpGLM_logoSum") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType)
        if (grepl("Clu",signature)){
          model <- bayesglm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)} 
        if (grepl("Uclu",signature)){
          model <- bayesglm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC(model)))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
   } 
    if (model_type=="bpGLM_logoSum_inter") {
      results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c(), Deviance=c(), Df=c(), nullDeviance=c(), nullDf=c(), AIC=c())
      for (gene in unique(germline$Gene.refGene)){
        print(gene)
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
        df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType)
        if (grepl("Clu",signature)){
          model <- bayesglm(Exposures ~ Mutation_Score + offset(log(Clustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)} 
        if (grepl("Uclu",signature)){
          model <- bayesglm(Exposures ~ Mutation_Score + offset(log(Unclustered)) + primaryTumorLocation*msStatus + tmbStatus + purity + ploidy + gender,family = poisson(),  data = df)}
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value, Deviance=model$deviance, Df=model$df.residual, nullDeviance=model$null.deviance, nullDf=model$df.null, AIC=AIC(model)))}
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
      results$PseudoR2 <- 1 - (results$Deviance / results$nullDeviance)
      write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="beta"){
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation); df$primaryTumorLocation = relevel(df$primaryTumorLocation, ref = "Other")
            df$LizaCancerType[df$LizaCancerType %in%  names(table(df$LizaCancerType)[table(df$LizaCancerType) < 10])] <- "Other"; df$LizaCancerType=factor(df$LizaCancerType)
            df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
            df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
            model <- betareg(Exposures ~ Mutation_Score + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender, data = df)
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients$mean["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients$mean["Mutation_Score", "Pr(>|z|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }
}
