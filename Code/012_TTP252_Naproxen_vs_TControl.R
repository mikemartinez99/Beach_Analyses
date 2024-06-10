# The purpose of this script is to analyze naproxen pirc tumor relative to control pirc tumor

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
library(dendextend)

# Generate output directory
opPath <- "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Set up output directory architecture for Rds files
rdsPath <- "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Set up output diretory architecture for DESeq2 files
ddsPath <- "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/DESeq2"
if (!dir.exists(ddsPath)) {
  dir.create(ddsPath)
}

# Set up output diretory architecture for GSEA files
gseaPath <- "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/GSEA"
if (!dir.exists(gseaPath)) {
  dir.create(gseaPath)
}

# CSet up output diretory architecture for custom figures
figDir <- "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/CustomFigures"
if (!dir.exists(figDir)){
  dir.create(figDir)
}

################################################################################

# Path to counts and metadata
tumorCounts <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Counts.csv"
tumorMetadata <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Metadata.csv"

# Specify comarison name
comparison <- "TP252_Combo vs Control"
reference <- "Control"
treatment <- "TP252_Naproxen"

# Read in the counts
tumor <- read.csv(tumorCounts, header = TRUE, sep = ",", row.names = 1)
colnames(tumor) <- sub("\\.", "-", colnames(tumor))
meta <- read.csv(tumorMetadata, header = TRUE, sep = ",", row.names = 1)

# !!!!!
# Because some groups in the normal tissue dataset have high within group variability, running DESeq2 with all groups included might inflate the 
# per-gene dispersion estimate for other groups. Therefore, we will subset the data beforehand, then run DESeq2 on the subset of groups.

# Get a vector of SampleIDs for Control and WT samples
groups <- meta[meta$Group == reference | meta$Group == treatment,]$Sample

# Subset meta and counts based on the `groups` vector
counts <- tumor[,colnames(tumor) %in% groups]
meta <- meta[rownames(meta) %in% groups,]

# Prepare DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Group)

# Set a reference level. in this case, we want to compare control-normal relative to WT
dds$Group <- relevel(dds$Group, ref = reference)

# Pre-filter for low counts
smallestGroup <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroup
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Sanity check that the comparison was run correctly. Output should be "Group_Control_vs_WT"
resultsNames(dds)

# Save the DESeq2 object as an Rds file
saveRDS(dds, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/Rds_Files/TTP252_Naproxen_vs_TControl_dds.Rds")

################################################################################

# Get the results from DESeq2
results <- as.data.frame(results(dds))
counts <- as.data.frame(counts(dds, normalized = TRUE))
results <- merge(results, counts, by = 0)
results$Symbols <- gsub("^[^-]+-(.*)$", "\\1", results$Row.names)
rownames(results) <- results$Row.names
results$Row.names <- NULL

# Remove unknown genes from the results
results <- results[!grepl("^ LOC\\d+$", results$Symbols),]
results <- results[!grepl("^ RGD\\d+$", results$Symbols),]
results <- results[!grepl("^ NA$", results$Symbols),]

# Save results file
write.csv(results, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/DESeq2/TTP252_Naproxen_vs_TControl_dds_results.csv")

# Get the significant results
results005 <- results[abs(results$log2FoldChange) > 2 & results$padj < 0.05,]
results005 <- results005[order(results005$log2FoldChange, decreasing = TRUE),]
numRes <- nrow(results005)
numRes <- numRes-19
results005top <- results005[c(1:20, numRes:nrow(results005)),]
results005top <- na.omit(results005top)
rownames(results005top) <- trimws(results005top$Symbols)
results005top$Symbols <- NULL
results005top[,c(3,4,5)] <- NULL
results005top$Significance <- ifelse(results005top$padj <= 0.001, "p < 0.001", "p < 0.05")

# Save top DEGs file
write.csv(results005top, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/DESeq2/TTP252_Naproxen_vs_TControl_Top_DEGs_Tabular_Data.csv")

################################################################################

# Variance stabilize transform the data
vsd <- vst(dds)

# Run PCA and return data for custom plotting
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- factor(PCA$Group, levels = c(treatment, reference))

custom_colors <- c("TP252_Naproxen" = "#66C2A5", "Control" = "#8DA0CB")


# Plot PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  coord_fixed() +
  theme_bw() +
  labs(title = paste(treatment, "vs", reference, sep = " ")) +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26),
        legend.text = element_text(size = 24),
        title = element_text(size = 26))
ggsave("Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/DESeq2/TTP252_Naproxen_vs_TControl_PCA.tiff", plot, width = 8, height = 8, dpi = 100)

################################################################################

# Prepare data for GSEA
# Extract Ensembl IDs
results$Ensembl <- gsub("^(.*?) - .*", "\\1", rownames(results))

# Map Entrez IDs to Ensembl IDss
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
write.csv(GSEA_input, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/GSEA/TTP252_Naproxen_vs_TControl_GSEA_Input_List.csv")


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

# Simplify the  redundant GO terms
GO.simp <- clusterProfiler::simplify(GO, cutoff = 0.6, by = "p.adjust", select_fun = min)


# Save gseaResult object as a .Rds object
saveRDS(GO.simp, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/Rds_Files/TTP252_Naproxen_vs_TControl_Simplified_gseGO.rds")

# Save results as a dataframe and write to csv
GO.df <- as.data.frame(setReadable(GO.simp, "org.Rn.eg.db", "ENTREZID"))
write.csv(GO.df, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/GSEA/TTP252_Naproxen_vs_TControl_Simplified_gseGO_results.csv")

# Now do GSEA for KEGG terms
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
saveRDS(KEGG, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/Rds_Files/TTP252_Naproxen_vs_TControl_gseKegg.rds")

# Save results as a dataframe and write to csv
KEGG.df <- as.data.frame(setReadable(KEGG, "org.Rn.eg.db", "ENTREZID"))
write.csv(KEGG.df, file = "Outputs/012_TTP252_Naproxen_vs_TControl_Outputs/GSEA/TTP252_Naproxen_vs_TControl_gseKEGG_results.csv")














