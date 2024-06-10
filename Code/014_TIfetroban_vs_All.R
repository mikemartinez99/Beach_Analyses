# Clear environment
rm(list = ls())

# Load libraries
library(dplyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(grDevices)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(ggrepel)
library(ComplexHeatmap)
library(Magick)
library(circlize)

# Generate output directory
opPath <- "Outputs/014_TIfetroban_vs_All_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Set up output directory architecture for Rds files
rdsPath <- "Outputs/014_TIfetroban_vs_All_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Set up output directory architecture for results files
resultsPath <- "Outputs/014_TIfetroban_vs_All_Outputs/DESeq2"
if (!dir.exists(resultsPath)) {
  dir.create(resultsPath)
}

# Set up output directory architecture for custom figures
figPath <- "Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures"
if (!dir.exists(figPath)) {
  dir.create(figPath)
}

# Set up output directory architecture Ifetroban data
dataPath <- "Data_Files/Ifetroban_Study"
if (!dir.exists(dataPath)) {
  dir.create(dataPath)
}

################################################################################

# Path to tumor counts and metadata
tumorCounts <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Counts.csv"
tumorMetadata <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Metadata.csv"

# Read in the counts
tumor <- read.csv(tumorCounts, header = TRUE, sep = ",", row.names = 1)
colnames(tumor) <- sub("\\.", "-", colnames(tumor))
meta <- read.csv(tumorMetadata, header = TRUE, sep = ",", row.names = 1)

# Filter out samples so each group has 5, 
# also remove the entire combination group for this particular comparison
remove <- c("T2-2443", "T3-2538", 
            "T6-2444", "T37", 
            "T39", "T41",
            "T27-2846", "T28-3102", "T29-4247", "T42", "T43", "T44", "T45")

# Remove the samples specified
counts <- tumor[,!(colnames(tumor) %in% remove)]
meta <- meta[!rownames(meta) %in% remove,]

# Remove TP252 and TP252_Naproxen
meta <- meta[meta$Group != "TP252",]
meta <- meta[meta$Group != "TP252_Naproxen",]

# Filter counts respectively
counts <- counts[,colnames(counts) %in% rownames(meta)]

# Save these counts and meta in the data directory
write.csv(counts, file = "Data_Files/Ifetroban_Study/IfetrobanStudy_MasterCounts.csv")
write.csv(meta, file = "Data_Files/Ifetroban_Study/IfetrobanStudy_MasterMeta.csv")

################################################################################

# Check that everything is the same between the two files
colnames(counts) == rownames(meta)

groups <- c("Control", "Ifetroban")
meta <- meta[meta$Group %in% groups,]
counts <- counts[,colnames(counts) %in% rownames(meta)]

# Prepare the DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Group)

# Set a reference level (in this case, it doesn't really matter)
dds$Group <- relevel(dds$Group, ref = "Control")

# Pre-filter for low counts
smallestGroup <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroup
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Save Rds
saveRDS(dds, file = "Outputs/014_TIfetroban_vs_All_Outputs/Rds_Files/TIfetroban_vs_TControl_dds.Rds")

################################################################################

# Extract results
results <- as.data.frame(results(dds))
results$Symbols <- gsub("^[^-]+-(.*)$", "\\1", rownames(dds))

# Extract the normalized counts
counts <- as.data.frame(counts(dds, normalized = TRUE))

# Merge results and counts
results <- merge(results, counts, by = 0)

# Remove unknown genes for each group
results <- results[!grepl("^ LOC\\d+$", results$Symbols),]
results <- results[!grepl("^ RGD\\d+$", results$Symbols),]
results <- results[!grepl("^ NA$", results$Symbols),]

# Save results as csv
write.csv(results, file = "Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_dds.csv")

# Extract the top results
results005 <- results[abs(results$log2FoldChange) > 2 & results$padj < 0.05,]
results005 <- results005[order(results005$log2FoldChange, decreasing = TRUE),]
results005 <- na.omit(results005)
numRes <- nrow(results005)
numRes <- numRes-19
results005top <- results005[c(1:20, numRes:nrow(results005)),]
rownames(results005top) <- trimws(results005top$Symbols)
results005top$Symbols <- NULL
results005top[,c(3,4,5)] <- NULL
results005top$Significance <- ifelse(results005top$padj <= 0.001, "p < 0.001", "p < 0.05")

# Save as csv
write.csv(results005top, file = "Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_Top_DEGs_Tabular.csv")

################################################################################

# All group PCA
# Variance-stabilization transformation
vsd <- vst(dds)

# Run PCA
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- factor(PCA$Group, levels = c("Control", "Ifetroban"))

# Assign colors
custom_colors <- c("Control" = "#66C2A5", "Ifetroban" = "#E78AC3")

# Plot the explortatory PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text_repel(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        legend.text = element_text(size = 24),
        title = element_text(size = 26))
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_PCA.tiff", plot, width = 8, height = 8, dpi = 100)

################################################################################

# Read in the ifetroban vs Control results
res <- read.csv("Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_dds.csv")

# Get the gene symbols
res$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", res$X))

# Set thresholds
significance_threshold <- 0.05
fc_threshold <- 1

# Create a grouping variable
res$Groups <- ifelse(res$padj < significance_threshold & res$log2FoldChange*-1 > fc_threshold, "Enriched in Control Tumor",
                     ifelse(res$padj < significance_threshold & res$log2FoldChange > fc_threshold, "Enriched in Ifetroban Tumor",
                            ifelse(res$padj < significance_threshold, "padj < 0.05", "ns")))
res$Groups <- factor(res$Groups, levels = c("Enriched in Control Tumor", "Enriched in Ifetroban Tumor",
                                            "padj < 0.05", "ns"))

# Plot custom volcano
Volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Groups)) +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(name = "",
                     values = c("Enriched in Control Tumor" = "steelblue", "Enriched in Ifetroban Tumor" = "firebrick2", 
                                "padj < 0.05" = "darkgrey", "ns" = "black")) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = fc_threshold*-1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change",
       y = "-log10(padj)") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24))
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/TIfetroban_vs_TControl_VolcanoPlot.tiff", Volcano, dpi = 300, width = 10, height = 10)

################################################################################
# Prepare Log10 TPM data

# Read in the gene lengths
lengths <- read.csv("Data_Files/Rattus_norvegicus_GTF/Rattus_norvegicus_exon_lengths_in_BP.csv")
rownames(lengths) <- lengths$X

# Read in the counts files
counts <- read.csv("Data_Files/Ifetroban_Study/IfetrobanStudy_MasterCounts.csv")
rownames(counts) <- counts$X
counts$X <- NULL
meta <- read.csv("Data_Files/Ifetroban_Study/IfetrobanStudy_MasterMeta.csv")
meta$Sample<- gsub("-", "\\.", meta$Sample)
rownames(meta) <- meta$Sample

# Filter for the correct groups
groups <- c("Control", "Ifetroban")
meta <- meta[meta$Group %in% groups,]
counts <- counts[,colnames(counts) %in% rownames(meta)]

# Get the ENSEMBL geneIDs
counts$ENSEMBL <- gsub("^(.*?)\\s-.*", "\\1", rownames(counts))
rownames(counts) <- counts$ENSEMBL
counts$ENSEMBL <- NULL

# Take the gene lengths for the genes in our dataset
genes <- rownames(counts)
lengths <- lengths[rownames(lengths) %in% genes,]

# Convert each gene length to KBP
kb <- lengths[,2] / 1000

# Divide each gene by it's length in KBP
rpk <- sweep(counts[,1:ncol(counts)], 1, kb, FUN="/")

# Sum up all of the RPKs for a sample
total_rpk <- colSums(rpk[,1:ncol(rpk)])
total_rpk_mill <- total_rpk/1000000

# Divide the RPK values by the total_rpk_mill (per million) scaling factor
tpm <- sweep(rpk[,1:ncol(rpk)], 2, total_rpk_mill, FUN="/")

# Add a small constant
tpm <- tpm + 0.01

# Log transform the tpm values
logTPM <- as.data.frame(t(apply(tpm,1, log2)))

# Get gene symbols
logTPM$Symbols <- mapIds(org.Rn.eg.db, key = rownames(logTPM),
                         column = "SYMBOL", keytype = "ENSEMBL",
                         multiVals = "first")

# Set up data directory to hold the log10 TPM counts
dataDir <- "Data_Files/Ifetroban_Study/Log10_TPM_Values"
if (dir.exists(dataDir)) {
  dir.create(dataDir)
}

# Save as csv
write.csv(logTPM, file = "TIfetroban_vs_TControl_Log10_TPM_Values.csv")

# Create a directory to hold stats information
statsDir <- "Outputs/014_TIfetroban_vs_All_Outputs/Stats"
if (!dir.exists(statsDir)) {
  dir.create(statsDir)
}


################################################################################

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

res <- read.csv("Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_dds.csv")
res$X <- NULL
res$Symbols <- trimws(res$Symbols)
res.filt <- res[res$Symbols %in% fibroGenes,]
rownames(res.filt) <- res.filt$Symbols
res.filt <- res.filt[fibroGenes,]
write.csv(res.filt, file = "Outputs/014_TIfetroban_vs_All_Outputs/Stats/Ifetroban_Fibroblast_Stats.csv")

# Isolate these genes from the logTPM dataframe
fibroCounts <- logTPM[logTPM$Symbols %in% fibroGenes,]
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
n <- 5
m <- 5
ref <- "Control"
treatment <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#E78AC3")
names(sample_colors) <- c(ref, treatment)

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:4, c(4,4,4,4))

#Define the number of slices in the heatmap
slices <- n+m

#Create a heatmap annotation
fibroAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#E78AC3"), fontisze = 14, fontface = "bold"), 
                     labels = c(ref, treatment),
                     labels_gp = gpar(col = "black", fontsize = 16, fontface = 2)),
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

#shades <- brewer.pal(name = "RdYlGn", n = 100)

#Heatmap for scaled data
fibroScaled <- Heatmap(fibroMat,
                       show_column_names = FALSE,
                       name = "Feature Scale",
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       top_annotation = fibroAnno,
                       left_annotation = fibroGroupAnno,
                       column_split = hmSplit,
                       column_title = NULL,
                       row_title = NULL,
                       row_split = rowSplit,
                       row_names_gp = gpar(fontsize = 14))
tiff("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Ifetroban_Fibroblast_Curated_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
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
immuneCounts <- logTPM[logTPM$Symbols %in% immuneGenes,]
rownames(immuneCounts) <- immuneCounts$Symbols
immuneCounts$Symbols <- NULL

# Get the stats associated with each curated gene
res <- read.csv("Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_dds.csv")
res$X <- NULL
res$Symbols <- trimws(res$Symbols)
res.filt <- res[res$Symbols %in% immuneGenes,]
rownames(res.filt) <- res.filt$Symbols
res.filt <- res.filt[immuneGenes,]
write.csv(res.filt, file = "Outputs/014_TIfetroban_vs_All_Outputs/Stats/Ifetroban_Immune_Stats.csv")

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
immuneCounts <- immuneCounts[immuneGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate min/max scaling on transposed data
immuneMat <- t(apply(immuneCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
n <- 5
m <- 5
ref <- "Control"
treatment <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#E78AC3")
names(sample_colors) <- c(ref, treatment)

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:4, c(5,4,4,5))

#Define the number of slices in the heatmap
slices <- n+m

#Create a heatmap annotation
immuneAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#E78AC3"), fontisze = 14, fontface = "bold"), 
                     labels = c(ref, treatment),
                     labels_gp = gpar(col = "black", fontsize = 16, fontface = 2)),
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

#shades <- brewer.pal(name = "RdYlGn", n = 100)

#Heatmap for scaled data
immuneScaled <- Heatmap(immuneMat,
                        show_column_names = FALSE,
                        name = "Feature Scale",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        top_annotation = immuneAnno,
                        left_annotation = immuneGroupAnno,
                        column_split = hmSplit,
                        column_title = NULL,
                        row_title = NULL,
                        row_split = rowSplit,
                        row_names_gp = gpar(fontsize = 14))
tiff("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Ifetroban_Immune_Curated_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
draw(immuneScaled)
dev.off()


















