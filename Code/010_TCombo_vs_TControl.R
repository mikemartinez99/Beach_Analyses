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
library(extrafont)
library(ComplexHeatmap)
library(magick)

# Generate output directory
opPath <- "Outputs/010_TCombo_vs_TControl_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Set up output directory architecture for Rds files
rdsPath <- "Outputs/010_TCombo_vs_TControl_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Set up output diretory architecture for DESeq2 files
ddsPath <- "Outputs/010_TCombo_vs_TControl_Outputs/DESeq2"
if (!dir.exists(ddsPath)) {
  dir.create(ddsPath)
}

# Set up output diretory architecture for GSEA files
gseaPath <- "Outputs/010_TCombo_vs_TControl_Outputs/GSEA"
if (!dir.exists(gseaPath)) {
  dir.create(gseaPath)
}

# CSet up output diretory architecture for custom figures
figDir <- "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures"
if (!dir.exists(figDir)){
  dir.create(figDir)
}

################################################################################

# Path to counts and metadata
tumorCounts <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Counts.csv"
tumorMetadata <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Metadata.csv"

# Specify comarison name
comparison <- "Combo vs Control"
reference <- "Control"
treatment <- "Combo"

# Read in the normal counts and associated metadata from `002_Exploratory_PCA.R`
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
saveRDS(dds, file = "Outputs/010_TCombo_vs_TControl_Outputs/Rds_Files/TCombo_vs_TControl_dds.Rds")

# Unhash this line if you are regenerating figures
# dds <- readRDS("Outputs/010_TCombo_vs_TControl_Outputs/Rds_Files/TCombo_vs_TControl_dds.Rds")

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
write.csv(results, file = "Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_dds_results.csv")

# Get the significant results
results005 <- results[abs(results$log2FoldChange) > 2 & results$padj < 0.05,]
results005 <- results005[order(results005$log2FoldChange, decreasing = TRUE),]
numRes <- nrow(results005)
numRes <- numRes-19
results005top <- results005[c(1:20, numRes:nrow(results005)),]
rownames(results005top) <- trimws(results005top$Symbols)
results005top$Symbols <- NULL
results005top[,c(3,4,5)] <- NULL
results005top$Significance <- ifelse(results005top$padj <= 0.001, "p < 0.001", "p < 0.05")

# Save top DEGs file
write.csv(results005top, file = "Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_Top_DEGs_Tabular_Data.csv")

################################################################################

# Variance stabilize transform the data
vsd <- vst(dds)

# Run PCA and return data for custom plotting
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- factor(PCA$Group, levels = c("Control", "Combo"))

custom_colors <- c("Control" = "#66C2A5", "Combo" = "#E78AC3")


# Plot PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(color = "",
       fill = "") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio = 1) +
  theme(legend.position = "bottom") +
  theme(text = element_text(family = "Times New Roman"),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26, face = "bold"),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 26, face = "bold"),
        legend.text = element_text(size = 24),
        title = element_text(size = 26))
ggsave("Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_PCA.tiff", plot, width = 10, height = 10, dpi = 300)

################################################################################

# Prepare data for GSEA
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
write.csv(GSEA_input, file = "Outputs/010_TCombo_vs_TControl_Outputs/GSEA/TCombo_vs_TControl_GSEA_Input_List.csv")


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
saveRDS(GO.simp, file = "Outputs/010_TCombo_vs_TControl_Outputs/Rds_Files/TCombo_vs_TControl_Simplified_gseGO.rds")

# Save results as a dataframe and write to csv
GO.df <- as.data.frame(setReadable(GO.simp, "org.Rn.eg.db", "ENTREZID"))
write.csv(GO.df, file = "Outputs/010_TCombo_vs_TControl_Outputs/GSEA/TCombo_vs_TControl_Simplified_gseGO_results.csv")

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
saveRDS(KEGG, file = "Outputs/010_TCombo_vs_TControl_Outputs/Rds_Files/TCombo_vs_TControl_gseKegg.rds")

# Save results as a dataframe and write to csv
KEGG.df <- as.data.frame(setReadable(KEGG, "org.Rn.eg.db", "ENTREZID"))
write.csv(KEGG.df, file = "Outputs/010_TCombo_vs_TControl_Outputs/GSEA/TCombo_vs_TControl_gseKEGG_results.csv")

################################################################################

# Generate curated GO dotplots

# Read in the curated GO and KEGG results
curGO <- read.csv("Data_Files/TCombo_vs_TControl_Curated/Curated_GSEA/TCombo_vs_TControl_Curated_GO.csv")
curKEGG <- read.csv("Data_Files/TCombo_vs_TControl_Curated/Curated_GSEA/TCombo_vs_TControl_Curated_KEGG.csv")

# Function to format the dataframe
formatGSEA <- function(x, reference, treatment) {
  formatted <- x %>%
    mutate(enrichment = ifelse(enrichmentScore > 0, paste("Upregulated in", treatment, sep = " "), paste("Downregulated in", treatment, sep = " "))) %>%
    mutate(GeneRatio = length(strsplit(as.character(core_enrichment), "/")) / setSize) %>%
    mutate(enrichment = factor(enrichment, levels = c(paste("Downregulated in", treatment, sep = " "), paste("Upregulated in", treatment, sep = " ")))) %>%
    arrange(enrichmentScore) 
  
  return(formatted)
}


# Format the dataframes and set factors
curGO <- formatGSEA(curGO, "Control", "Combination")
curGO$Description <- factor(curGO$Description, levels = curGO$Description)

curKEGG <- formatGSEA(curKEGG, "Control", "Combination")
curKEGG$Description <- factor(curKEGG$Description, levels = curKEGG$Description)

# Function to plot dotplot
plotDot <- function(df) {
  plot <- ggplot(df, aes(x = enrichmentScore, y = Description, color = p.adjust)) +
    geom_point(aes(size = setSize),alpha = 0.8) +
    scale_color_continuous(low = "red", high = "blue") +
    scale_size(range = c(4,10)) +
    facet_grid(~enrichment, scales = "free") +
    xlab("Enrichment Score") +
    theme(axis.text.y = element_text(size = 9)) +
    labs(x = "Enrichment Score",
         y = "",
         color = "P adjust",
         size = "Set Size") +
    theme_bw() +
    theme(text = element_text(family = "Times New Roman"),
          axis.text.x = element_text(size = 10, angle = 90),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 18),
          strip.text = element_text(size = 14, face = "bold"),
          title= element_text(size = 20),
          legend.text = element_text(size = 12)) 
  return(plot)
}

# Set up output directory for custom GSEA plots
customGSEAplots <- "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Custom_GSEA"
if (!dir.exists(customGSEAplots)) {
  dir.create(customGSEAplots)
}

# Plot dotplots
GOdotplot <- plotDot(curGO)
ggsave("Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Custom_GSEA/TCombo_vs_TControl_Curated_GO_dotplot.tiff", GOdotplot, width = 12, height = 10, dpi = 300)

KEGGdotplot <- plotDot(curKEGG)
ggsave("Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Custom_GSEA/TCombo_Vs_TControl_Curated_KEGG_dotplot.tiff", KEGGdotplot, width = 12, height = 10, dpi = 300)

################################################################################

# Prepare for TPM heatmap
# Read in the gene lengths
lengths <- read.csv("Data_Files/Rattus_norvegicus_GTF/Rattus_norvegicus_exon_lengths_in_BP.csv")
rownames(lengths) <- lengths$X

# Get ENSEMBL gene IDs
counts$ENSEMBL <- gsub("^(.*?)\\s-.*", "\\1", rownames(counts))
rownames(counts) <- counts$ENSEMBL
counts$ENSEMBL <- NULL

# Take the gene lengths for the all the genes in our counts data
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

# Save as a csv
# Set up output diretory architecture for curated data
curatedDataDir <- "Data_Files/TCombo_vs_TControl_Curated/Log10_TPM_Values"
if (!dir.exists(curatedDataDir)){
  dir.create(curatedDataDir)
}
write.csv(logTPM, file = "Data_Files/TCombo_vs_TControl_Curated/Log10_TPM_Values/TCombo_vs_TControl_Log10_TPM_Values.csv")

# Unhash the following line if you are regenerating figures
# logTPM <- read.csv("Data_Files/TCombo_vs_TControl_Curated/Log10_TPM_Values/TCombo_vs_TControl_Log10_TPM_Values.csv")

################################################################################
# Set base font
par(family = "Times New Roman")

# Set up stats directory
statsDir <- "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/Heatmap_Statistics/"
if (!dir.exists(statsDir)){
  dir.create(statsDir)
}

# Fibroblast Heatmap
fibroblasts <- read.csv("Data_Files/TCombo_vs_TControl_Curated/Custom_Heatmap_Gene/Fibroblasts_Combo.csv")
fibroblasts$Specific.Cell <- factor(fibroblasts$Specific.Cell, levels = c("Fibroblast",
                                                                          "Myofibroblast",
                                                                          "ECM fibroblast",
                                                                          "Inflammatory CAF"))

# Set rownames
fibroblasts <- fibroblasts %>%
  column_to_rownames(var = "Genes")

# Get a list of genes we want to keep
fibroGenes <- rownames(fibroblasts)


# Get associated stats
res <- read.csv("Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_dds_results.csv")
rownames(res) <- res$X
res$X <- NULL
res$Symbols <- trimws(res$Symbols)
res.filt <- res[res$Symbols %in% fibroGenes,]
rownames(res.filt) <- res.filt$Symbols
res.filt <- res.filt[fibroGenes,]
write.csv(res.filt, file = "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/Heatmap_Statistics/Fibroblast_Stats_Combo.csv")

# Isolate these genes from the logTPM dataframe
fibroCounts <- logTPM[logTPM$Symbols %in% fibroGenes,]
fibroCounts$X <- NULL
rownames(fibroCounts) <- fibroCounts$Symbols
fibroCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
fibroCounts <- fibroCounts[fibroGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate Zscore on log10 TPM data
fibroMat <- t(apply(fibroCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
n <- 7
m <- 7
ref <- "Untreated"
treatment <- "Combination"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#E78AC3")
names(sample_colors) <- c(ref, treatment)

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:4, c(4,4,6,6))

#Define the number of slices in the heatmap
slices <- n+m

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

#Create a heatmap annotation
fibroAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#E78AC3"), fontisze = 14, fontface = "bold", fontfamily = "Times"), 
                     labels = c(ref, treatment),
                     labels_gp = gpar(col = "black", fontsize = 16, fontface = 2, fontfamily = "Times")),
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

# Create directory for custom heatmaps
customHeatmapDir <- "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps"
if (!dir.exists(customHeatmapDir)) {
  dir.create(customHeatmapDir)
}


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
                       row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
tiff("Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/Fibroblast_Log10_TPM_Heatmap_Combo.tiff", width = 8, height = 10, units = "in", res = 300)
draw(fibroScaled)
dev.off()

################################################################################

immune <- read.csv("Data_Files/TCombo_vs_TControl_Curated/Custom_Heatmap_Gene/Immune_Activation_Combo.csv")
immune <- immune %>%
  column_to_rownames(var = "Genes") %>%
  mutate(Pathway = factor(Pathway, levels = c("ROS generating", "Immune cell recruitment", "Immune Suppressive", "Tumor promoting cytokines", "Innate immune activation")))

# Get a list of genes we want to keep
immuneGenes <- rownames(immune)

# Get associated stats
res <- read.csv("Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_dds_results.csv")
rownames(res) <- res$X
res$X <- NULL
res$Symbols <- trimws(res$Symbols)
res.filt <- res[res$Symbols %in% immuneGenes,]
rownames(res.filt) <- res.filt$Symbols
res.filt <- res.filt[immuneGenes,]
write.csv(res.filt, file = "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/Heatmap_Statistics/Immune_Stats_Combo.csv")


# Isolate these genes from the logTPM dataframe
immuneCounts <- logTPM[logTPM$Symbols %in% immuneGenes,]
immuneCounts$X <- NULL
rownames(immuneCounts) <- immuneCounts$Symbols
immuneCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
immuneCounts <- immuneCounts[immuneGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate Zscore on log10 TPM data
immuneMat <- t(apply(immuneCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
n <- 7
m <- 7
ref <- "Untreated"
treatment <- "Combination"

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:5, c(3,3,3,4,9))

#Define the number of slices in the heatmap
slices <- n+m

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))


#Create a heatmap annotation
immuneAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#E78AC3"), fontisze = 14, fontface = "bold", fontfamily = "Times"), 
                     labels = c(ref, treatment),
                     labels_gp = gpar(col = "black", fontsize = 16, fontface = 2, fontfamily = "Times")),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

# Choose colors
colors <- brewer.pal(8, "Paired")

# Create the row annotation
immuneGroupAnno <- rowAnnotation(
  Pathway = immune$Pathway,
  col = list(Pathway = c("ROS generating" = "#A6CEE3",
                         "Immune cell recruitment" = "#B2DF8A",
                         "Immune Suppressive" = "#FB9A99",
                         "Tumor promoting cytokines" = "#FDBF6F",
                         "Innate immune activation" = "#1F78B4")),
  show_annotation_name = FALSE
)

#Heatmap for scaled data
immuneScaled <- Heatmap(immuneMat,
                        column_labels = colnames(immuneMat), 
                        show_column_names = FALSE,
                        name = "Feature Scale",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        column_title = NULL,
                        row_title = NULL,
                        top_annotation = immuneAnno,
                        left_annotation = immuneGroupAnno,
                        column_split = hmSplit,
                        row_split = rowSplit,
                        row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
tiff("Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/Immune_Log10_TPM_Heatmap_Combo.tiff", width = 8, height = 10, units = "in", res = 300)
draw(immuneScaled)
dev.off()

################################################################################

# Cell fate differentiation heatmap
differ <- read.csv("Data_Files/TCombo_vs_TControl_Curated/Custom_Heatmap_Gene/Cell_Fate_Combo.csv")
differ <- differ %>%
  column_to_rownames(var = "Genes") %>%
  mutate(Cell.type = factor(Cell.type, levels = c("Enterocyte", "Enteroendocrine", "Crypt fibroblast", "Goblet cell", "Stem cell")))

# Specify the elements to keep
differGenes <- rownames(differ)

# Get associated stats
res <- read.csv("Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_dds_results.csv")
rownames(res) <- res$X
res$X <- NULL
res$Symbols <- trimws(res$Symbols)
res.filt <- res[res$Symbols %in% differGenes,]
rownames(res.filt) <- res.filt$Symbols
res.filt <- res.filt[differGenes,]
write.csv(res.filt, file = "Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/Heatmap_Statistics/Differ_Stats_Combo.csv")


# Isolate these genes from the logTPM dataframe
differCounts <- logTPM[logTPM$Symbols %in% differGenes,]
differCounts$X <- NULL
rownames(differCounts) <- differCounts$Symbols
differCounts$Symbols <- NULL

# Make sure the fibroCounts dataframe is in the same order as thr fibroblast dataframe
differCounts <- differCounts[differGenes,]

# Calculate Z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate Z-score on log10 TPM data
differMat <- t(apply(differCounts, 1, cal_z_score))

# Define slices (n = reference, m = treatment) and labels
n <- 7
m <- 7
ref <- "Untreated"
treatment <- "Combination"

# Set splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:5, c(4,4,3,4,4))

#Define the number of slices in the heatmap
slices <- n+m

# Set global options for heatmaps
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

#Create a heatmap annotation
differAnno <- HeatmapAnnotation(
  Group = anno_block(gp = gpar(fill = c("#66C2A5", "#E78AC3"), fontisze = 14, fontface = "bold", fontfamily = "Times"), 
                     labels = c(ref, treatment),
                     labels_gp = gpar(col = "black", fontsize = 16, fontface = 2, fontfamily = "Times")),
  col = list(Group = sample_colors),
  show_annotation_name = FALSE
)

# Use the Dark2 palette for the colors of groups
colors <- brewer.pal(8, "Dark2")

# Create the row annotation
differGroupAnno <- rowAnnotation(
  `Cell Type` = differ$Cell.type,
  col = list(`Cell Type` = c("Enterocyte" = "#1B9E77",
                             "Enteroendocrine" = "#D95F02",
                             "Crypt fibroblast" = "#E7298A",
                             "Goblet cell" = "#E6AB02",
                             "Stem cell" = "#A6761D")),
  show_annotation_name = FALSE
)


# Heatmap for scaled data
differScaled <- Heatmap(differMat,
                        column_labels = colnames(differMat), 
                        show_column_names = FALSE,
                        name = "Feature Scale",
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        column_title = NULL,
                        row_title = NULL,
                        top_annotation = differAnno,
                        left_annotation = differGroupAnno,
                        column_split = hmSplit,
                        row_split = rowSplit,
                        row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
tiff("Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/Curated_Heatmaps/CellFate_Log10_TPM_Heatmap_Combo.tiff", width = 8, height = 10, units = "in", res = 300)
draw(differScaled)
dev.off()

################################################################################

# Read in the deseq2 output file from this run
res <- read.csv("Outputs/007_TControl_vs_NControl_Outputs/DESeq2/TControl_vs_NControl_dds_results.csv", header = TRUE, sep = ",")

# Get the gene symbols
res$Symbols <- gsub("^[^-]+-(.*)$", "\\1", res$X)

# Set thresholds
significance_threshold <- 0.05
fc_threshold <- 1

# Create a grouping variable
res$Groups <- ifelse(res$padj < significance_threshold & res$log2FoldChange*-1 > fc_threshold, "Enriched in Control Normal",
                     ifelse(res$padj < significance_threshold & res$log2FoldChange > fc_threshold, "Enriched in Control Tumor",
                            ifelse(res$padj < significance_threshold, "padj < 0.05", "ns")))
res$Groups <- factor(res$Groups, levels = c("Enriched in Control Normal", "Enriched in Control Tumor",
                                            "padj < 0.05", "ns"))

# Be aware there is a leading white space character before each gene name!!!!! Remove it!
res$Symbols <- sub(" ", "", res$Symbols)

# Now make a vector to label the genes we want
genes_to_label <- c("Mmp7", "Lgr5", "Prss22", "Nox1", "Defa6", "S100a9",
                    "Vip", "Ghrl", "Myh11", "Ptgis", "Car3", "Scn7a")

################################################################################

# Volcano plot
# Read in the results file
res <- read.csv("Outputs/010_TCombo_vs_TControl_Outputs/DESeq2/TCombo_vs_TControl_dds_results.csv")

# Get the gene symbols
res$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", res$X))

# Set thresholds
significance_threshold <- 0.05
fc_threshold <- 1

# Create a grouping variable
res$Groups <- ifelse(res$padj < significance_threshold & res$log2FoldChange*-1 > fc_threshold, "Down regulated in Combination Tumor",
                     ifelse(res$padj < significance_threshold & res$log2FoldChange > fc_threshold, "Up regulated in Combination Tumor",
                            ifelse(res$padj < significance_threshold, "padj < 0.05", "ns")))
res$Groups <- factor(res$Groups, levels = c("Down regulated in Combination Tumor", "Up regulated in Combination Tumor",
                                            "padj < 0.05", "ns"))

# Label a subset of genes
genesToLabel <- c("Reg3b", "Csf3", "S100a9", "Il17a", "Nos2", "Mmp7", "Cxcl1", "COX2", # Combo Down
                  "Aqp8", "Pck1", "Grem2", "Ghrl", "Cnn3", "Myh11") # Combo Up

# Plot custom volcano
Volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Groups)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_text_repel(data = res[res$Symbols %in% genesToLabel,], aes(label = Symbols), nudge_y = 0.5,
                  color = "black",
                  bg.color = "white",
                  bg.r = 0.15,
                  box.padding = unit(0.60, "lines")) +
  scale_color_manual(name = "",
                     values = c("Down regulated in Combination Tumor" = "steelblue", "Up regulated in Combination Tumor" = "firebrick2", 
                                "padj < 0.05" = "darkgrey", "ns" = "black")) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = fc_threshold*-1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = c(-10, -5, -1, 0, 1, 5),  # Specify only the breaks you want
                     limits = c(-10, 6)) +
  labs(x = "Log2 Fold Change",
       y = "-log10(padj)") +
  theme_classic() +
  theme(text = element_text(family = "Times New Roman"),
        legend.position = "bottom",
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24, face = "bold")) 
ggsave("Outputs/010_TCombo_vs_TControl_Outputs/CustomFigures/TCombo_vs_TControl_VolcanoPlot.tiff", Volcano, dpi = 300, width = 12, height = 12)





