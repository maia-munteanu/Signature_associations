library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(pscl)
library(betareg)

args=commandArgs(TRUE)
signature <- as.character(args[1])
input <- as.character(args[2])
model_type <- as.character(args[3])
metadata <- as.data.frame(fread(args[4])); metadata$LizaCancerType=gsub("_MSI","",metadata$LizaCancerType)

germline=as.data.frame(fread(input))
germline %>% pivot_longer(cols = -c(sample, Gene.refGene, Freq), names_to = "Signature", values_to = "Exposures") %>% filter(Signature==signature) -> germline
germline %>% left_join(metadata[,c("sample","gender","purity","ploidy","msStatus","tmbStatus","LizaCancerType")]) -> germline

if (model_type=="GLMnb") {
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(germline$Gene.refGene)){
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        model <- glm.nb(Exposures ~ Mutation_Score + gender + purity + ploidy + msStatus + tmbStatus + LizaCancerType, data = df)
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
} 

if (model_type=="LM"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(germline$Gene.refGene)){
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        model <- lm(Exposures ~ Mutation_Score + gender + purity + ploidy + msStatus + tmbStatus + LizaCancerType, data = df)
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
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        df$Exposures[df$Exposures == 0] <- df$Exposures[df$Exposures == 0] + 0.00001
        df$Exposures[df$Exposures == 1] <- df$Exposures[df$Exposures == 1] - 0.00001
        model <- betareg(Exposures ~ Mutation_Score + gender + purity + ploidy + msStatus + tmbStatus + LizaCancerType, data = df)
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients$mean["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients$mean["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

