library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(pscl)

args=commandArgs(TRUE)
signature <- as.character(args[1])
input <- as.character(args[2])
output_info <- as.character(args[3])

germline=as.data.frame(fread(input))
germline %>% pivot_longer(cols = -c(sample, Gene.refGene, Freq), names_to = "Signature", values_to = "Exposures") %>% filter(Signature==signature) -> germline

if (model=="GLMnb") {
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(germline$Gene.refGene)){
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        model <- glm.nb(Exposures ~ Mutation_Score, data = df)
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
        write.table(results, file = paste0(signature, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
} 

if (model=="LM"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(germline$Gene.refGene)){
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        model <- lm(Exposures ~ Mutation_Score, data = df)
        beta <- coef(model)["Mutation_Score"]
        se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
}

if (model=="ZInb"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(germline$Gene.refGene)){
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
        model <- zeroinfl(Exposures ~ Mutation_Score | Mutation_Score, data = data, dist = "negbin")
        beta <- coef(model)["count_Mutation_Score"]
        se <- summary(model)$coefficients$count["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients$count["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
}

if (model=="Hnb"){
    results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
    for (gene in unique(germline$Gene.refGene)){
        df<-germline[which(germline$Gene.refGene == gene),]
        df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)

        model <- hurdle(Exposures ~ Mutation_Score | Mutation_Score, data = data, dist = "negbin")
        beta <- coef(model)["count_Mutation_Score"]
        se <- summary(model)$coefficients$count["Mutation_Score", "Std. Error"]
        p_value <- summary(model)$coefficients$count["Mutation_Score", "Pr(>|z|)"]
        results<-rbind(results,data.frame(Signature = signature, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
        results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
}
