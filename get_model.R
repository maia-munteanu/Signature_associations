library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(pscl)
library(betareg)
library(arm)
library(statmod)

args=commandArgs(TRUE)
signature <- as.character(args[1])
input <- as.character(args[2])
model_type <- as.character(args[3])
metadata <- as.data.frame(fread(args[4])); metadata$LizaCancerType=gsub("_MSI","",metadata$LizaCancerType)
covariates <- as.logical(args[5])

germline=as.data.frame(fread(input))
germline %>% pivot_longer(cols = -c(sample, Gene.refGene, Freq, Clustered, Unclustered), names_to = "Signature", values_to = "Exposures") %>% filter(Signature==signature) -> germline
germline %>% left_join(metadata[,c("sample","gender","purity","ploidy","msStatus","tmbStatus","primaryTumorLocation")]) -> germline

if (covariates) {
    if (model_type=="GLMglog2") {
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation)
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
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation)
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
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation)
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
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation)
            if (grepl("Clu",signature)){model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = 1.5, link.power = 1),  data = df)} 
            else {model <- model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered) + primaryTumorLocation + msStatus + tmbStatus + purity + ploidy + gender,family = tweedie(var.power = 1.5, link.power = 1),  data = df)}
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="beta"){
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$primaryTumorLocation[df$primaryTumorLocation %in%  names(table(df$primaryTumorLocation)[table(df$primaryTumorLocation) < 10])] <- "Other"; df$primaryTumorLocation=factor(df$primaryTumorLocation)
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
} else {
    if (model_type=="GLMglog2") {
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            model <- glm(log2(Exposures + 1) ~ Mutation_Score, family = gaussian(), data = df)
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
            if (grepl("Clu",signature)){model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered), family = gaussian(), data = df)} else {model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered), family = gaussian(), data = df)}
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
            if (grepl("Clu",signature)){model <- bayesglm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered), family = gaussian(), data = df)} else {model <- bayesglm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered), family = gaussian(), data = df)}
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
            if (grepl("Clu",signature)){model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Clustered),family = tweedie(var.power = 1.5, link.power = 1),  data = df)} else {model <- model <- glm(log2(Exposures + 1) ~ Mutation_Score + log(Unclustered),family = tweedie(var.power = 1.5, link.power = 1),  data = df)}
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    } 
    if (model_type=="beta"){
        results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
        for (gene in unique(germline$Gene.refGene)){
            print(gene)
            df<-germline[which(germline$Gene.refGene == gene),]
            df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
            df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
            df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
            model <- betareg(Exposures ~ Mutation_Score, data = df)
            beta <- coef(model)["Mutation_Score"]
            se <- summary(model)$coefficients$mean["Mutation_Score", "Std. Error"]
            p_value <- summary(model)$coefficients$mean["Mutation_Score", "Pr(>|z|)"]
            results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }
}
