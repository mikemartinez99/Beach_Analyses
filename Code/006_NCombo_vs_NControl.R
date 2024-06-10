# The purpose of this script is to analyze normal Pirc Combination EPA and Naproxen treated relative to normal control pirc

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(AnnotationDbi)

# Generate output directory
opPath <- "Outputs/006_NCombo_vs_NControl_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Set up output directory architecture for Rds files
rdsPath <- "Outputs/006_NCombo_vs_NControl_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Set up output diretory architecture for DESeq2 files
ddsPath <- "Outputs/006_NCombo_vs_NControl_Outputs/DESeq2"
if (!dir.exists(ddsPath)) {
  dir.create(ddsPath)
}

# Set up output diretory architecture for GSEA files
gseaPath <- "Outputs/006_NCombo_vs_NControl_Outputs/GSEA"
if (!dir.exists(gseaPath)) {
  dir.create(gseaPath)
}

################################################################################


# Read in the normal counts and associated metadata from `002_Exploratory_PCA.R`
normal <- read.csv("Data_Files/Filtered_Normal_Counts/Final_Filtered_Normal_Counts.csv", header = TRUE, sep = ",", row.names = 1)
colnames(normal) <- sub("\\.", "-", colnames(normal))
meta <- read.csv("Data_Files/Filtered_Normal_Counts/Final_Filtered_Normal_Metadata.csv", header = TRUE, sep = ",", row.names = 1)
meta$SampleID <- rownames(meta)

# Get a vector of SampleIDs for Control and WT samples
groups <- meta[meta$Group == "Combo" | meta$Group == "Control",]$SampleID

# Subset meta and counts based on the `groups` vector
counts <- normal[,colnames(normal) %in% groups]
meta <- meta[rownames(meta) %in% groups,]

# Prepare DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Group)


# Set a reference level. in this case, we want to compare control-normal relative to WT
dds$Group <- relevel(dds$Group, ref = "Control")

# Pre-filter for low counts
smallestGroup <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroup
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Sanity check that the comparison was run correctly. Output should be "Group_Control_vs_WT"
resultsNames(dds)

# Save the dds object
saveRDS(dds, file = "Outputs/006_NCombo_vs_NControl_Outputs/Rds_Files/NCombo_vs_NControl_dds.Rds")

# Variance stabilize transform the data
vsd <- vst(dds)

# Run PCA and return data for custom plotting
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- factor(PCA$Group, levels = c("Control", "EPA"))

# Plot PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = "N Combo vs  N Control") +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        legend.text = element_text(size = 24),
        title = element_text(size = 26))
ggsave("Outputs/006_NCombo_vs_NControl_Outputs/DESeq2/NCombo_vs_NControl_PCA.tiff", plot, width = 12, height = 8, dpi = 100)

################################################################################

# Let's get the results and the normalized counts
results <- as.data.frame(results(dds))
counts <- as.data.frame(counts(dds, normalized = TRUE))
results <- merge(results, counts, by = 0)
rownames(results) <- results$Row.names
results$Row.names <- NULL

# Write results to csv
write.csv(results, file = "Outputs/006_NCombo_vs_NControl_Outputs/DESeq2/NCombo_vs_NControl_DESeq2_Results.csv")

################################################################################
# GSEA

# Extract Ensembl IDs
results$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(results))

# Map Entrez IDs to Ensembl IDs
results$Entrez <- mapIds(org.Rn.eg.db, key = results$Ensembl,
                         column = "ENTREZID", keytype = "ENSEMBL",
                         multiVals = "first")

# Order the genes and remove NA values
results <- results[order(results$log2FoldChange, decreasing = TRUE),]
results <- na.omit(results)

# Get a vector of fold changes and name them by their Entrez ID
genes <- results$log2FoldChange
names(genes) <- results$Entrez

# Remove any duplicate Entrez IDs and for sanity, order in decreasing order again
unique_entrez_genes <- names(genes[!duplicated(names(genes))])
unique_genes <- genes[unique_entrez_genes]
unique_genes <- sort(unique_genes, decreasing = TRUE)

# Check how many genes there are
length(unique_genes)

# Save the named list as a GSEA input
GSEA_input <- as.data.frame(unique_genes)
GSEA_input$Entrez <- rownames(GSEA_input)
GSEA_input$Ensembl <- mapIds(org.Rn.eg.db, key = GSEA_input$Entrez,
                             column = "ENSEMBL", keytype = "ENTREZID",
                             multiVals = "first")
GSEA_input$Symbols <- mapIds(org.Rn.eg.db, key = GSEA_input$Ensembl,
                             column = "SYMBOL", keytype = "ENSEMBL",
                             multiVals = "first")
write.csv(GSEA_input, file = "Outputs/006_NCombo_vs_NControl_Outputs/GSEA/NCombo_vs_NControl_GESA_Input_List.csv")


# Set seed for reproducibility
set.seed(1234)

# Run GSEA for GO terms
GO <- gseGO(unique_genes,
            ont = "all",
            OrgDb = "org.Rn.eg.db",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            eps = 1e-300,
            verbose = TRUE,
            by = "fgsea",
            seed = TRUE)


# Save gseaResult object as a .Rds object
saveRDS(GO, file = "Outputs/006_NCombo_vs_NControl_Outputs/Rds_Files/NCombo_vs_NControl_WT_gseGO.Rds")

# Save results as a dataframe and write to csv
GO.df <- as.data.frame(setReadable(GO, "org.Rn.eg.db", "ENTREZID"))
write.csv(GO.df, file = "Outputs/006_NCombo_vs_NControl_Outputs/GSEA/NCombo_vs_NControl_gseGO_results.csv")

# Now do GSEA for KEGG terms
# Set seed for reproductibility
set.seed(1234)

# Run GSEA for KEGG terms
KEGG <- gseKEGG(unique_genes,
                organism = "rno",
                keyType = "kegg",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                eps = 1e-300,
                verbose = TRUE,
                by = "fgsea",
                seed = TRUE)

# Save gseaResult object as a .Rds file
saveRDS(KEGG, file = "Outputs/006_NCombo_vs_NControl_Outputs/Rds_Files/NCombo_vs_NControl_gseKEGG.Rds")

# Save results as a dataframe and write to csv
KEGG.df <- as.data.frame(setReadable(KEGG, "org.Rn.eg.db", "ENTREZID"))
write.csv(KEGG.df, file = "Outputs/006_NCombo_vs_NControl_Outputs/GSEA/NCombo_vs_NControlgseKEGG_results.csv")


