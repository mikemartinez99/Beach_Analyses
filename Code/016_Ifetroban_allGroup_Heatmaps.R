# The purpose of this code is to generate Heatmaps for Ifetroban (for all groups)
# based on what was showin in script 014.

# Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(circlize)
library(magick)
library(ComplexHeatmap)
library(RColorBrewer)
library(org.Rn.eg.db)
library(AnnotationDbi)

# Clear environment
rm(list = ls())

# Initilize a function to generate new output directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# Generate an output directory
OP <- newDir("Outputs/016_Ifetroban_allGroup_Heatmaps_Outputs")

# Custom Heatmap directory
figDir <- newDir("Outputs/016_Ifetroban_allGroup_Heatmaps_Outputs/Custom_Figures")

################################################################################

#----- Read in the tumor counts data (containing all groups)
counts <- read.csv("Data_Files/Ifetroban_Study/IfetrobanStudy_MasterCounts.csv")
counts <- counts %>%
  column_to_rownames(var = "X")

# Read in the metadata
meta <- read.csv("Data_Files/Ifetroban_Study/IfetrobanStudy_MasterMeta.csv")

# Adjust Sample column in meta to match the counts
meta$Sample <- gsub("-", ".", meta$Sample)
rownames(meta) <- meta$Sample

#----- Generate TPM values for the counts

# Read in the gene length file
lengths <- read.csv("Data_Files/Rattus_norvegicus_GTF/Rattus_norvegicus_exon_lengths_in_BP.csv")

# Get gene symbol and format the gene column to match the counts
lengths$Symbols <- mapIds(org.Rn.eg.db,
                          column = "SYMBOL",
                          keys = lengths$Gene,
                          keytype = "ENSEMBL")
rownames(lengths) <- paste(lengths$Gene, lengths$Symbols, sep = " - ")

# Convert lengths to KB
lengths$exon.lengths <- lengths[,2] / 1000

# Check that the genes in the counts are in the same order as the lengths
all(rownames(counts) %in% rownames(lengths))
all(rownames(counts) == rownames(lengths))

# Function to generate TPM values
countToTPM <- function(count_df, length_df) {
  rate <- count_df/length_df
  tpm <- c()
  for (i in 1:nrow(counts)) {
    tpm[i] <- rate[i]/sum(rate) * 1e6
  }
  return(tpm)
}

# Run the function on the counts to generate TPMs
TPMs <- apply(counts, 2, countToTPM, lengths$exon.lengths)
TPMs <- as.data.frame(TPMs)

# Add a small constant
TPMs <- TPMs + 0.01

# Log transform the tpm values
logTPMs <- as.data.frame(t(apply(TPMs,1, log2)))

# Set rownames
rownames(TPMs) <- rownames(counts)
rownames(logTPMs) <- rownames(TPMs)

# Get gene symbols
TPMs$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(TPMs))

# Save as a csv
write.csv(TPMs, file = "Data_Files/Ifetroban_Study/AllGroup_TPMs.csv")
write.csv(logTPMs, file = "Data_Files/Ifetroban_Study/AllGroup_Log2TPMs.csv")

################################################################################


#----- Fibroblast Heatmap

# Get symbols for logTPMs
logTPMs$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(logTPMs)))

# Read in the fibroblast heatmap genes
fibroblasts <- read.csv("Data_Files/Ifetroban_Study/Custom_Heatmap_Gene/Ifetroban_Fibroblasts.csv")
fibroblasts$Specific.Cell <- factor(fibroblasts$Specific.Cell, levels = c("Inflammatory CAF",
                                                                          "Myofibroblast",
                                                                          "Crypt bottom fibroblast",
                                                                          "TGF-Beta signaling"))

# Set rownames
fibroblasts <- fibroblasts %>%
  column_to_rownames(var = "Genes")

# Get a list of genes we want to keep
fibroGenes <- rownames(fibroblasts)

# Isolate these genes from the logTPM dataframe
fibroCounts <- logTPMs[logTPMs$Symbols %in% fibroGenes,]
fibroCounts$X <- NULL
rownames(fibroCounts) <- fibroCounts$Symbols
fibroCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
fibroCounts <- fibroCounts[fibroGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate min/max scaling on transposed data
fibroMat <- t(apply(fibroCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
a <- 5 # Untreated (control)
b <- 5 # Naproxen
c <- 5 # EPA
d <- 5 # Ifetroban

A <- "Untreated"
B <- "Naproxen"
C <- "EPA"
D <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
names(sample_colors) <- c(A, B, C, D)

# Set heatmap splitting patterns for sample groups
hmSplit <- rep(1:4, c(a, b, c, d))

# Set heatmap splitting pattern for gene groups
rowSplit <- rep(1:4, c(4,4,4,4))

# Define the number of slices in the heatmap
slices <- a + b + c + d

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

# Create a heatmap annotation
fibroAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), 
                               fontisze = 14, 
                               fontface = "bold", 
                               fontfamily = "Times"), 
                     labels = c(A, B, C, D),
                     labels_gp = gpar(col = "black", 
                                      fontsize = 16, 
                                      fontface = 2, 
                                      fontfamily = "Times")),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

# Specify the brewer pal
cellTypeCols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C")
cellTypeLevels <- levels(fibroblasts$Specific.Cell)
names(cellTypeCols) <- cellTypeLevels

# Create the row annotation
fibroGroupAnno <- rowAnnotation(
  `Cell Type` = fibroblasts$Specific.Cell,
  col = list(`Cell Type` = cellTypeCols),
  show_annotation_name = FALSE
)

#Heatmap for scaled data
fibroScaled <- Heatmap(fibroMat,
                       show_column_names = FALSE,
                       name = "Z Score",
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       top_annotation = fibroAnno,
                       left_annotation = fibroGroupAnno,
                       column_split = hmSplit,
                       column_title = NULL,
                       row_title = NULL,
                       row_split = rowSplit,
                       row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
png("Outputs/016_Ifetroban_allGroup_Heatmaps_Outputs/Custom_Figures/Ifetroban_All_Groups_Fibroblast_All_Samples_HM.png", width = 8, height = 10, units = "in", res = 300)
draw(fibroScaled)
dev.off()

################################################################################

# Average the expresion between the samples per group: Control
control <- logTPMs[,1:5]
controlMean <- as.data.frame(rowMeans(control))
colnames(controlMean) <- c("Untreated")

# Naproxen
naproxen <- logTPMs[,6:10]
naproxenMean <- as.data.frame(rowMeans(naproxen))
colnames(naproxenMean) <- c("Naproxen")

# EPA
EPA <- logTPMs[,11:15]
EPAMean <- as.data.frame(rowMeans(EPA))
colnames(EPAMean) <- c("EPA")

# Ifetroban
ifetroban <- logTPMs[,16:20]
ifetrobanMean <- as.data.frame(rowMeans(ifetroban))
colnames(ifetrobanMean) <- c("Ifetroban")

# Cbind the 4 mean dataframes
meanData <- cbind(cbind(cbind(controlMean, naproxenMean), EPAMean), ifetrobanMean)

# Get gene symbols for meanData
meanData$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(meanData)))

# Isolate these genes from the logTPM dataframe
fibroCounts <- meanData[meanData$Symbols %in% fibroGenes,]
fibroCounts$X <- NULL
rownames(fibroCounts) <- fibroCounts$Symbols
fibroCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
fibroCounts <- fibroCounts[fibroGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate min/max scaling on transposed data
fibroMat <- t(apply(fibroCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
a <- 1 # Untreated (control)
b <- 1 # Naproxen
c <- 1 # EPA
d <- 1 # Ifetroban

A <- "Untreated"
B <- "Naproxen"
C <- "EPA"
D <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
names(sample_colors) <- c(A, B, C, D)

# Set heatmap splitting patterns for sample groups
hmSplit <- rep(1:4, c(a, b, c, d))

# Set heatmap splitting pattern for gene groups
rowSplit <- rep(1:4, c(4,4,4,4))

# Define the number of slices in the heatmap
slices <- a + b + c + d

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

# Create a heatmap annotation
fibroAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), 
                               fontisze = 14, 
                               fontface = "bold", 
                               fontfamily = "Times"), 
                     labels = c(A, B, C, D),
                     labels_gp = gpar(col = "black", 
                                      fontsize = 16, 
                                      fontface = 2, 
                                      fontfamily = "Times")),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

# Specify the brewer pal
cellTypeCols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C")
cellTypeLevels <- levels(fibroblasts$Specific.Cell)
names(cellTypeCols) <- cellTypeLevels

# Create the row annotation
fibroGroupAnno <- rowAnnotation(
  `Cell Type` = fibroblasts$Specific.Cell,
  col = list(`Cell Type` = cellTypeCols),
  show_annotation_name = FALSE
)

#Heatmap for scaled data
fibroScaled <- Heatmap(fibroMat,
                       show_column_names = FALSE,
                       name = "Z Score",
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       top_annotation = fibroAnno,
                       left_annotation = fibroGroupAnno,
                       column_split = hmSplit,
                       column_title = NULL,
                       row_title = NULL,
                       row_split = rowSplit,
                       row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
png("Outputs/016_Ifetroban_allGroup_Heatmaps_Outputs/Custom_Figures/Ifetroban_All_Groups_Fibroblast_Averaged_Samples_HM.png", width = 8, height = 10, units = "in", res = 300)
draw(fibroScaled)
dev.off()

################################################################################

# Read in the fibroblast heatmap genes
immune <- read.csv("Data_Files/Ifetroban_Study/Custom_Heatmap_Gene/Ifetroban_Immune.csv")
immune$Specific.Cell <- factor(immune$Specific.Cell, levels = c("T-cell activation",
                                                                "Immune recruitment",
                                                                "T-cell exhaustion",
                                                                "Immune suppression"))

# Set rownames
immune <- immune %>%
  column_to_rownames(var = "Genes")

# Get a list of genes we want to keep
immuneGenes <- rownames(immune)

# Isolate these genes from the logTPM dataframe
immuneCounts <- logTPMs[logTPMs$Symbols %in% immuneGenes,]
immuneCounts$X <- NULL
rownames(immuneCounts) <- immuneCounts$Symbols
immuneCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
immuneCounts <- immuneCounts[immuneGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate min/max scaling on transposed data
immuneMat <- t(apply(immuneCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
a <- 5 # Untreated (control)
b <- 5 # Naproxen
c <- 5 # EPA
d <- 5 # Ifetroban

A <- "Untreated"
B <- "Naproxen"
C <- "EPA"
D <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
names(sample_colors) <- c(A, B, C, D)

# Set heatmap splitting patterns for sample groups
hmSplit <- rep(1:4, c(a, b, c, d))

# Set heatmap splitting pattern for gene groups
rowSplit <- rep(1:4, c(5,4,3,5))

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

#Create a heatmap annotation
immuneAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), 
                               fontisze = 14, 
                               fontface = "bold", 
                               fontfamily = "Times"), 
                     labels = c(A, B, C, D),
                     labels_gp = gpar(col = "black", 
                                      fontsize = 16, 
                                      fontface = 2, 
                                      fontfamily = "Times")),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

# Specify the brewer pal
cellTypeCols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C")
cellTypeLevels <- levels(immune$Specific.Cell)
names(cellTypeCols) <- cellTypeLevels

# Create the row annotation
immuneGroupAnno <- rowAnnotation(
  `Cell Type` = immune$Specific.Cell,
  col = list(`Cell Type` = cellTypeCols),
  show_annotation_name = FALSE
)


#Heatmap for scaled data
immuneScaled <- Heatmap(immuneMat,
                        show_column_names = FALSE,
                        name = "Z Score",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        top_annotation = immuneAnno,
                        left_annotation = immuneGroupAnno,
                        column_split = hmSplit,
                        column_title = NULL,
                        row_title = NULL,
                        row_split = rowSplit,
                        row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
png("Outputs/016_Ifetroban_allGroup_Heatmaps_Outputs/Custom_Figures/Ifetroban_All_Groups_Immune_All_Samples_HM.png", width = 8, height = 10, units = "in", res = 300)
draw(immuneScaled)
dev.off()

################################################################################

# Isolate these genes from the logTPM dataframe
immuneCounts <- meanData[meanData$Symbols %in% immuneGenes,]
immuneCounts$X <- NULL
rownames(immuneCounts) <- immuneCounts$Symbols
immuneCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
immuneCounts <- immuneCounts[immuneGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate min/max scaling on transposed data
immuneMat <- t(apply(immuneCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
a <- 1 # Untreated (control)
b <- 1 # Naproxen
c <- 1 # EPA
d <- 1 # Ifetroban

A <- "Untreated"
B <- "Naproxen"
C <- "EPA"
D <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
names(sample_colors) <- c(A, B, C, D)

# Set heatmap splitting patterns for sample groups
hmSplit <- rep(1:4, c(a, b, c, d))

# Set heatmap splitting pattern for gene groups
rowSplit <- rep(1:4, c(5,4,3,5))

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

#Create a heatmap annotation
immuneAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"), 
                               fontisze = 14, 
                               fontface = "bold", 
                               fontfamily = "Times"), 
                     labels = c(A, B, C, D),
                     labels_gp = gpar(col = "black", 
                                      fontsize = 16, 
                                      fontface = 2, 
                                      fontfamily = "Times")),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

# Specify the brewer pal
cellTypeCols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C")
cellTypeLevels <- levels(immune$Specific.Cell)
names(cellTypeCols) <- cellTypeLevels

# Create the row annotation
immuneGroupAnno <- rowAnnotation(
  `Cell Type` = immune$Specific.Cell,
  col = list(`Cell Type` = cellTypeCols),
  show_annotation_name = FALSE
)

#Heatmap for scaled data
immuneScaled <- Heatmap(immuneMat,
                        show_column_names = FALSE,
                        name = "Z Score",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        top_annotation = immuneAnno,
                        left_annotation = immuneGroupAnno,
                        column_split = hmSplit,
                        column_title = NULL,
                        row_title = NULL,
                        row_split = rowSplit,
                        row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
png("Outputs/016_Ifetroban_allGroup_Heatmaps_Outputs/Custom_Figures/Ifetroban_All_Groups_Immune_Averaged_Samples_HM.png", width = 8, height = 10, units = "in", res = 300)
draw(immuneScaled)
dev.off()

#







