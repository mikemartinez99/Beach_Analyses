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
library(magick)
library(circlize)
library(extrafont)


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

# Set up output directory architecture for GSEA
gseaPath <- "Outputs/014_TIfetroban_vs_All_Outputs/GSEA"
if (!dir.exists(gseaPath)) {
  dir.create(gseaPath)
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

# Unhash this line if you are regenerating figures
# dds <- readRDS("Outputs/014_TIfetroban_vs_All_Outputs/Rds_Files/TIfetroban_vs_TControl_dds.Rds")

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

# Prepare data for GSEA
# Extract Ensembl IDs
rownames(results) <- results$Row.names
results$Row.names <- NULL
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
write.csv(GSEA_input, file = "Outputs/014_TIfetroban_vs_All_Outputs/GSEA/TIfetroban_vs_TControl_GSEA_Input_List.csv")


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
GO.simp <- clusterProfiler::simplify(GO, cutoff = 0.6, by = "p.adjust", select_fun = min)

# Save Rds
saveRDS(GO.simp, file = "Outputs/014_TIfetroban_vs_All_Outputs/Rds_Files/TIfetroban_vs_TControl_Simplified_gseGO.rds")

# Save gseGO results as a dataframe
GO.df <- as.data.frame(setReadable(GO.simp, org.Rn.eg.db, "ENTREZID"))
write.csv(GO.df, file = "Outputs/014_TIfetroban_vs_All_Outputs/GSEA/TIfetroban_vs_TControl_Simplified_gseGO.csv")

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

# Save Rds
saveRDS(KEGG, file = "Outputs/014_TIfetroban_vs_All_Outputs/Rds_Files/TIfetroban_vs_TControl_gseKEGG.rds")

# Save gseGO results as a dataframe
KEGG.df <- as.data.frame(setReadable(KEGG, org.Rn.eg.db, "ENTREZID"))
write.csv(KEGG.df, file = "Outputs/014_TIfetroban_vs_All_Outputs/GSEA/TIfetroban_vs_TControl_gseKEGG.csv")

################################################################################
# Curated GSEA figures
# Read in the curated GO and KEGG results
curGO <- read.csv("Data_Files/Ifetroban_Study/TIfetroban_vs_TControl_Curated/Curated_GSEA/TIfetroban_vs_TControl_Curated_GO.csv")
curKEGG <- read.csv("Data_Files/Ifetroban_Study/TIfetroban_vs_TControl_Curated/Curated_GSEA/TIfetroban_vs_TControl_Curated_KEGG.csv")

reference = "Control"
treatment = "Ifetroban"

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
curGO <- formatGSEA(curGO, "Control", "Ifetroban")
curGO$Description <- factor(curGO$Description, levels = curGO$Description)

curKEGG <- formatGSEA(curKEGG, "Control", "Ifetroban")
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
customGSEAplots <- "Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Custom_GSEA"
if (!dir.exists(customGSEAplots)) {
  dir.create(customGSEAplots)
}

# Plot dotplots
GOdotplot <- plotDot(curGO)
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Custom_GSEA/TIfetroban_vs_TControl_Curated_GO_dotplot.tiff", GOdotplot, width = 14, height = 10, dpi = 300)

KEGGdotplot <- plotDot(curKEGG)
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Custom_GSEA/TIfetroban_Vs_TControl_Curated_KEGG_dotplot.tiff", KEGGdotplot, width = 14, height = 10, dpi = 300)

################################################################################

# All group PCA
# Variance-stabilization transformation
vsd <- vst(dds)

# Run PCA
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- ifelse(PCA$Group == "Control", "Untreated", "Ifetroban")
PCA$Group <- factor(PCA$Group, levels = c("Untreated", "Ifetroban"))

# Assign colors
custom_colors <- c("Untreated" = "#66C2A5", "Ifetroban" = "#E78AC3")

# Plot the explortatory PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text_repel(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  labs(color = "",
       fill = "") +
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
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_PCA.tiff", plot, width = 10, height = 10, dpi = 300)

################################################################################

# Create a directory to hold custom volcano
volcDir <- "Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Custom_Volcanos"
if (!dir.exists(volcDir)) {
  dir.create(volcDir)
}

# Read in the ifetroban vs Control results
res <- read.csv("Outputs/014_TIfetroban_vs_All_Outputs/DESeq2/TIfetroban_vs_TControl_dds.csv")
rownames(res) <- res$Row.names
res$Row.names <- NULL

# Get the gene symbols
res$Symbols <- trimws(gsub("^[^-]+-(.*)$", "\\1", rownames(res)))

# Set thresholds
significance_threshold <- 0.05
fc_threshold <- 1

# Create a grouping variable
res$Groups <- ifelse(res$padj < significance_threshold & res$log2FoldChange*-1 > fc_threshold, "Downregulated in Ifetroban",
                     ifelse(res$padj < significance_threshold & res$log2FoldChange > fc_threshold, "Upregulated in Ifetroban",
                            ifelse(res$padj < significance_threshold, "padj < 0.05", "ns")))
res$Groups <- factor(res$Groups, levels = c("Downregulated in Ifetroban", "Upregulated in Ifetroban",
                                            "padj < 0.05", "ns"))

# Omit NA values
res <- na.omit(res)

# Plot custom volcano
Volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Groups)) +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(name = "",
                     values = c("Downregulated in Ifetroban" = "steelblue", "Upregulated in Ifetroban" = "firebrick2", 
                                "padj < 0.05" = "darkgrey", "ns" = "black")) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = fc_threshold*-1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change",
       y = "-log10(padj)") +
  theme_classic() +
  scale_x_continuous(breaks = c(-5, -1, 0, 1, 5),  # Specify only the breaks you want
                     limits = c(-8, 5)) +
  labs(x = "Log2 Fold Change",
       y = "-log10(padj)") +
  theme(text = element_text(family = "Times New Roman"),
        legend.position = "bottom",
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24, face = "bold")) 
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Custom_Volcanos/TIfetroban_vs_TControl_VolcanoPlot.tiff", Volcano, dpi = 300, width = 12, height = 12)

# Labelled Volcano plot
# Order the results so we can see the top hits
orderedResults <- res[order(res$padj, decreasing = FALSE),][1:8,]$Symbols
orderedFCpos <- res[order(res$log2FoldChange, decreasing = TRUE),][1:8,]$Symbols
orderedFCneg <- res[order(res$log2FoldChange, decreasing = FALSE),][1:8,]$Symbols

# Make a vector of genes to label
genesToLabel <- c(orderedResults, orderedFCpos, orderedFCneg)

# Remove Tdh from the label since it is not significant
indexRemove <- c("Tdh")
genesToLabel <- genesToLabel[genesToLabel != indexRemove]
indexRemove <- c("Cyp2c24")
genesToLabel <- genesToLabel[genesToLabel != indexRemove]

padj_thresh <- -log10(0.05)
pos_fc_thresh <- 1
neg_fc_thresh <- -1

# Make custom colors for plotting
custom_colors <- c("Enriched in Ifetroban Tumor" = "firebrick2",
                   "Enriched in Control Tumor" = "steelblue",
                   "padj < 0.05" = "grey",
                   "ns" = "black")


# Plot volcano
Volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Groups)) +
  geom_point() +
  geom_text_repel(data = res[res$Symbols %in% genesToLabel,], aes(label = Symbols), nudge_y = 0.5,
                  box.padding = unit(0.35, "lines")) +
  geom_vline(xintercept = pos_fc_thresh, linetype = "dashed") +
  geom_vline(xintercept = neg_fc_thresh, linetype = "dashed") +
  geom_hline(yintercept = padj_thresh, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20))
ggsave("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Custom_Volcanos/TIfetroban_vs_TControl_Labelled_Volcano_Plot.tiff", Volcano, width = 12, height = 12, dpi = 300)

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
write.csv(logTPM, file = "Data_Files/Ifetroban_Study/Log10_TPM_Values/TIfetroban_vs_TControl_Log10_TPM_Values.csv")

# Unhash this line if you are editing figures
# logTPM <- read.csv("Data_Files/Ifetroban_Study/Log10_TPM_Values/TIfetroban_vs_TControl_Log10_TPM_Values.csv")

# Create a directory to hold stats information
statsDir <- "Outputs/014_TIfetroban_vs_All_Outputs/Stats"
if (!dir.exists(statsDir)) {
  dir.create(statsDir)
}


################################################################################
# Create a directory to hold curated heatmaps
hmDir <- "Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Curated_Heatmaps"
if (!dir.exists(hmDir)) {
  dir.create(hmDir)
}

# Unhash the following line if you are regenerating figures
logTPM <- read.csv("Data_Files/Ifetroban_Study/Log10_TPM_Values/TIfetroban_vs_TControl_Log10_TPM_Values.csv")

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
n <- 5
m <- 5
ref <- "Untreated"
treatment <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#E78AC3")
names(sample_colors) <- c(ref, treatment)

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:4, c(4,4,4,4))

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
                       row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
tiff("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Curated_Heatmaps/Ifetroban_Fibroblast_Curated_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
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
immuneCounts$X <- NULL
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
ref <- "Untreated"
treatment <- "Ifetroban"

# Assign sample group colors
sample_colors <- c("#66C2A5", "#E78AC3")
names(sample_colors) <- c(ref, treatment)

#Set heatmap splitting pattern
hmSplit <- rep(1:2, c(n, m))
rowSplit <- rep(1:4, c(5,4,3,5))

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
                        row_names_gp = gpar(fontsize = 14, fontfamily = "Times"))
tiff("Outputs/014_TIfetroban_vs_All_Outputs/CustomFigures/Curated_Heatmaps/Ifetroban_Immune_Curated_Heatmap.tiff", width = 8, height = 10, units = "in", res = 300)
draw(immuneScaled)
dev.off()

################################################################################









