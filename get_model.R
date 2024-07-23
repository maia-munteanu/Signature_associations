library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)

args=commandArgs(TRUE)
signature <- as.character(args[1])
input <- as.character(args[2])
output_info <- as.character(args[3])

germline=as.data.frame(fread(input))
germline %>% pivot_longer(cols = -c(sample, Gene.refGene, Freq), names_to = "Signatures", values_to = "Exposures") -> germline
germline[which(germline$Signatures==signature),] -> germline

results=data.frame(Signature = c(), Gene = c(), Beta = c(), SE = c(), P_Value = c())
for (gene in unique(germline$Gene.refGene)){
    df<-germline[which(germline$Gene.refGene == gene & germline$Signatures == sig),]
    df$Mutation_Score <- ifelse(df$Freq > 0, 1, 0)
    model <- glm.nb(Exposures ~ Mutation_Score, data = df)
    beta <- coef(model)["Mutation_Score"]
    se <- summary(model)$coefficients["Mutation_Score", "Std. Error"]
    p_value <- summary(model)$coefficients["Mutation_Score", "Pr(>|t|)"]
    results<-rbind(results,data.frame(Signature = sig, Gene = gene, Beta = beta, SE = se, P_Value = p_value))}
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")
    write.table(results, file = paste0("GLMnb_", output_info,"_", sig, ".tsv"),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
