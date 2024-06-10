# The purpose of this code is to compare the GSEA clustersf for the PREVENT analysis

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(DOSE)
library(ggplot2)
library(ggh4x)

# Generate output directory
opPath <- "Outputs/013_GSEA_Cluster_Comparison_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Set up output directory architecture for Rds files
rdsPath <- "Outputs/013_GSEA_Cluster_Comparison_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Set up output directory architecture for results files
resultsPath <- "Outputs/013_GSEA_Cluster_Comparison_Outputs/Results"
if (!dir.exists(resultsPath)) {
  dir.create(resultsPath)
}

# Set up output directory architecture for custom figures
figPath <- "Outputs/013_GSEA_Cluster_Comparison_Outputs/CustomFigures"
if (!dir.exists(figPath)) {
  dir.create(figPath)
}



# Read in the input lists
controlVnormalGSEA <- read.csv("Outputs/007_TControl_vs_NControl_Outputs/GSEA/TControl_vs_NControl_GESA_Input_List.csv")
controlVnormalGSEA <- controlVnormalGSEA[controlVnormalGSEA$Symbols != "2510039O18Rikl",]
controlVnormalGSEA <- controlVnormalGSEA[controlVnormalGSEA$Symbols != "6430548M08Rikl",]
controlVnormalGSEA <- controlVnormalGSEA[!grepl("^LOC\\d+$", controlVnormalGSEA$Symbols),]
controlVnormalGSEA <- controlVnormalGSEA[!grepl("^RGD\\d+$", controlVnormalGSEA$Symbol),]
colnames(controlVnormalGSEA) <- c("Entrez", "Log2FC", "entrez", "Ensembl", "Symbols")
one <- controlVnormalGSEA$Log2FC
names(one) <- controlVnormalGSEA$Entrez

# Read in the naproxen vs control list
napVcontrol <- read.csv("Outputs/008_TNaproxen_vs_TControl_Outputs/GSEA/TNaproxen_vs_TControl_GSEA_Input_List.csv")
napVcontrol<- napVcontrol[napVcontrol$Symbols != "2510039O18Rikl",]
napVcontrol <- napVcontrol[napVcontrol$Symbols != "6430548M08Rikl",]
napVcontrol <- napVcontrol[!grepl("^LOC\\d+$", napVcontrol$Symbols),]
napVcontrol <- napVcontrol[!grepl("^RGD\\d+$", napVcontrol$Symbol),]
colnames(napVcontrol) <- c("Entrez", "Log2FC", "entrez", "Ensembl", "Symbols")
two <- napVcontrol$Log2FC
names(two) <- napVcontrol$Entrez

# Read in the epa vs control list
epaVcontrol <- read.csv("Outputs/009_TEPA_vs_TControl_Outputs/GSEA/TEPA_vs_TControl_GSEA_Input_List.csv")
epaVcontrol <- epaVcontrol[epaVcontrol$Symbols != "2510039O18Rikl",]
epaVcontrol <- epaVcontrol[epaVcontrol$Symbols != "6430548M08Rikl",]
epaVcontrol <- epaVcontrol[!grepl("^LOC\\d+$", epaVcontrol$Symbols),]
epaVcontrol <- epaVcontrol[!grepl("^RGD\\d+$", epaVcontrol$Symbol),]
colnames(epaVcontrol) <- c("Entrez", "Log2FC", "entrez", "Ensembl", "Symbols")
three <- epaVcontrol$Log2FC
names(three) <- epaVcontrol$Entrez

# Read in the combo vs control
comboVcontrol <- read.csv("Outputs/010_TCombo_vs_TControl_Outputs/GSEA/TCombo_vs_TControl_GSEA_Input_List.csv")
comboVcontrol <- comboVcontrol[comboVcontrol$Symbols != "2510039O18Rikl",]
comboVcontrol <- comboVcontrol[comboVcontrol$Symbols != "6430548M08Rikl",]
comboVcontrol <- comboVcontrol[!grepl("^LOC\\d+$", comboVcontrol$Symbols),]
comboVcontrol <- comboVcontrol[!grepl("^RGD\\d+$", comboVcontrol$Symbol),]
colnames(comboVcontrol) <- c("Entrez", "Log2FC", "entrez", "Ensembl", "Symbols")
four <- comboVcontrol$Log2FC
names(four) <- comboVcontrol$Entrez

# Initialize a list to hold all the GSEA inputs
comparisonLists <- list(Control_T = one, Naproxen = two, EPA = three, Combo = four)
str(comparisonLists)

# Run compareCluster function for KEGG
test.out <- compareCluster(geneClusters=comparisonLists,  fun = "gseKEGG", organism = "rno")
test.out <- setReadable(test.out, "org.Rn.eg.db", "ENTREZID")
test.out@compareClusterResult$Group <- ifelse(test.out@compareClusterResult$enrichmentScore < 0, "Suppressed", "Enriched")
result <- test.out@compareClusterResult
result$Group <- ifelse(result$enrichmentScore < 0, "Suppressed", "Enriched")
resultPlot <- dotplot(test.out, showCategory = 15, split = "Cluster")  + facet_nested(~Cluster + Group, scales = "free") +
  theme(axis.text.x = element_blank())
ggsave("Outputs/013_GSEA_Cluster_Comparison_Outputs/CustomFigures/PREVENT_GSEA_Cluster_Comparison_Dotplot.tiff", resultPlot, width = 20, height = 12, dpi = 300)

# Save Rds and results
saveRDS(test.out, file = "Outputs/013_GSEA_Cluster_Comparison_Outputs/Rds_Files/PREVENT_GSEA_Cluster_Comparison.Rds")
write.csv(result, file = "Outputs/013_GSEA_Cluster_Comparison_Outputs/Results/PREVENT_GSEA_Cluster_Comparison_Results.csv")



