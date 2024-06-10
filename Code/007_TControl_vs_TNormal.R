# The purpose of this script is to analyze control pirc tumor relative to control pirc normal
# Also included in this code is a custom volcano plot highlighting several key differentially regulated genes

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
library(magick)
library(ComplexHeatmap)
library(circlize)

# Generate output directory
opPath <- "Outputs/007_TControl_vs_NControl_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Set up output directory architecture for Rds files
rdsPath <- "Outputs/007_TControl_vs_NControl_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Set up output diretory architecture for DESeq2 files
ddsPath <- "Outputs/007_TControl_vs_NControl_Outputs/DESeq2"
if (!dir.exists(ddsPath)) {
  dir.create(ddsPath)
}

# Set up output diretory architecture for GSEA files
gseaPath <- "Outputs/007_TControl_vs_NControl_Outputs/GSEA"
if (!dir.exists(gseaPath)) {
  dir.create(gseaPath)
}

# CSet up output diretory architecture for custom figures
figDir <- "Outputs/007_TControl_vs_NControl_Outputs/CustomFigures"
if (!dir.exists(figDir)){
  dir.create(figDir)
}

################################################################################

# Path to counts and metadata
normalCounts <- "Data_Files/Filtered_Normal_Counts/Final_Filtered_Normal_Counts.csv"
normalMetadata <- "Data_Files/Filtered_Normal_Counts/Final_Filtered_Normal_Metadata.csv"
tumorCounts <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Counts.csv"
tumorMetadata <- "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Metadata.csv"

# Read in the normal counts and associated metadata from `002_Exploratory_PCA.R`
normal <- read.csv(normalCounts, header = TRUE, sep = ",", row.names = 1)
colnames(normal) <- sub("\\.", "-", colnames(normal))
meta <- read.csv(normalMetadata, header = TRUE, sep = ",", row.names = 1)
meta$SampleID <- rownames(meta)
tumor <- read.csv(tumorCounts, header = TRUE, sep = ",", row.names = 1)
colnames(tumor) <- sub("\\.", "-", colnames(tumor))
tumormeta <- read.csv(tumorMetadata, header = TRUE, sep = ",", row.names = 1)
tumormeta$SampleID <- rownames(tumormeta)

# !!!!!
# Because some groups in the normal tissue dataset have high within group variability, running DESeq2 with all groups included might inflate the 
# per-gene dispersion estimate for other groups. Therefore, we will subset the data beforehand, then run DESeq2 on the subset of groups.
normalGroup <- meta[meta$Group == "Control",]$SampleID
tumorGroup <- tumormeta[tumormeta$Group == "Control",]$SampleID

# Subset meta and counts based on the `groups` vector
normcounts <- normal[,colnames(normal) %in% normalGroup]
meta <- meta[rownames(meta) %in% colnames(normcounts),]
meta$Group <- "Control Normal"
tumcounts <- tumor[,colnames(tumor) %in% tumorGroup]
tumormeta <- tumormeta[rownames(tumormeta) %in% colnames(tumcounts),]
tumormeta$Group <- "Control Tumor"

# Merge the counts
counts <- merge(normcounts, tumcounts, by = 0)
rownames(counts) <- counts$Row.names
counts$Row.names <- NULL

# Merge the metadata
met <- rbind(meta, tumormeta)

# Prepare DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = met,
                              design = ~ Group)

# Specify comarison name
comparison <- "Control Tumor vs Control Normal"
reference <- "Control Normal"
treatment <- "Control Tumor"


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

# Save the dds object
saveRDS(dds, "Outputs/007_TControl_vs_NControl_Outputs/Rds_Files/TControl_vs_NControl_dds.Rds")

# Variance stabilize transform the data
vsd <- vst(dds)

# Run PCA and return data for custom plotting
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- factor(PCA$Group, levels = c(treatment, reference))

# Plot PCA
plot <- ggplot(PCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text(size = 4, aes(label = name), hjust = 1, vjust = 1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
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
ggsave("Outputs/007_TControl_vs_NControl_Outputs/DESeq2/TControl_vs_NControl_PCA.tiff", plot, width = 12, height = 8, dpi = 100)

################################################################################

# Let's get the results and the normalized counts
results <- as.data.frame(results(dds))
counts <- as.data.frame(counts(dds, normalized = TRUE))
results <- merge(results, counts, by = 0)
rownames(results) <- results$Row.names
results$Row.names <- NULL

write.csv(results, file = "Outputs/007_TControl_vs_NControl_Outputs/DESeq2/TControl_vs_NControl_dds_results.csv")

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
write.csv(GSEA_input, file = "Outputs/007_TControl_vs_NControl_Outputs/GSEA/TControl_vs_NControl_GESA_Input_List.csv")

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
saveRDS(GO, file = "Outputs/007_TControl_vs_NControl_Outputs/Rds_Files/TControl_vs_NControl_Simplified_gseGO.Rds")

# Save results as a dataframe and write to csv
GO.df <- as.data.frame(setReadable(GO, "org.Rn.eg.db", "ENTREZID"))
write.csv(GO.df, file = "Outputs/007_TControl_vs_NControl_Outputs/GSEA/TControl_vs_NControl_Simplified_gseGO_results.csv")


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
saveRDS(KEGG, file = "Outputs/007_TControl_vs_NControl_Outputs/Rds_Files//TControl_vs_NControl_gseKEGG.Rds")

# Save results as a dataframe and write to csv
KEGG.df <- as.data.frame(setReadable(KEGG, "org.Rn.eg.db", "ENTREZID"))
write.csv(KEGG.df, file = "Outputs/007_TControl_vs_NControl_Outputs/GSEA/TControl_vs_NControl_gseKEGG_results.csv")

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


# Plot custom volcano
Volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Groups)) +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(name = "",
                     values = c("Enriched in Control Normal" = "steelblue", "Enriched in Control Tumor" = "firebrick2", 
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
        axis.title.y = element_text(size = 24)) +
  geom_text_repel(data = subset(res, Symbols %in% genes_to_label), 
                  aes(label = Symbols), 
                  nudge_y = 0.5, 
                  nudge_x = 0.5,
                  box.padding = 1,
                  segment.size = 0.5,
                  segment.color = "black",
                  size = 6,
                  color = "black",
                  fontface = 4) 
ggsave("Outputs/007_TControl_vs_NControl_Outputs/CustomFigures/TControl_vs_NControl_Curated_VolcanoPlot.tiff", Volcano, dpi = 300, width = 10, height = 10)

################################################################################

# Generate curated GO dotplots

# Read in the curated GO and KEGG results
curGO <- read.csv("Data_Files/TControl_vs_NControl_Curated/Curated_GSEA/TControl_vs_NControl_Curated_GO.csv")
curKEGG <- read.csv("Data_Files/TControl_vs_NControl_Curated/Curated_GSEA/TControl_vs_NControl_Curated_KEGG.csv")

# Function to format the dataframe
formatGSEA <- function(x, reference, treatment) {
  formatted <- x %>%
    mutate(enrichment = ifelse(enrichmentScore > 0, paste("Enriched in", treatment, sep = " "), paste("Enriched in", reference, sep = " "))) %>%
    mutate(GeneRatio = length(strsplit(as.character(core_enrichment), "/")) / setSize) %>%
    mutate(enrichment = factor(enrichment, levels = c(paste("Enriched in", reference, sep = " "), paste("Enriched in", treatment, sep = " ")))) %>%
    arrange(enrichmentScore) 
  
  return(formatted)
}

# Format the dataframes and set factors
curGO <- formatGSEA(curGO, "Tumor", "Normal")
curGO$Description <- factor(curGO$Description, levels = curGO$Description)

curKEGG <- formatGSEA(curKEGG, "Tumor", "Normal")
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
         y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, angle = 90),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 18),
          strip.text = element_text(size = 14, face = "bold"),
          title= element_text(size = 20),
          legend.text = element_text(size = 12)) 
  return(plot)
}

# Set up output directory for custom GSEA plots
customGSEAplots <- "Outputs/007_TControl_vs_NControl_Outputs//CustomFigures/Custom_GSEA"
if (!dir.exists(customGSEAplots)) {
  dir.create(customGSEAplots)
}

# Plot dotplots
GOdotplot <- plotDot(curGO)
ggsave("Outputs/007_TControl_vs_NControl_Outputs/CustomFigures/Custom_GSEA/TControl_vs_NControl_Curated_GO_dotplot.tiff", GOdotplot, width = 12, height = 10, dpi = 300)

KEGGdotplot <- plotDot(curKEGG)
ggsave("Outputs/007_TControl_vs_NControl_Outputs/CustomFigures/Custom_GSEA/TControl_Vs_NControl_Curated_KEGG_dotplot.tiff", KEGGdotplot, width = 12, height = 10, dpi = 300)

# Function to plot barplot
plotbar <- function(df) {
  plot <- ggplot(df, aes(x = enrichmentScore, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(~enrichment, scales = "free") +
    scale_fill_continuous(low = "red", high = "blue") +
    labs(x = "Enrichment Score",
         y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, angle = 90),
          axis.title.x = element_text(size = 20),
          axis.text.y = element_text(size = 18),
          strip.text = element_text(size = 16, face = "bold"),
          title= element_text(size = 20),
          legend.text = element_text(size = 12))
  return(plot)
}


GOdotplot <- plotbar(curGO)
ggsave("Outputs/007_TControl_vs_NControl_Outputs/CustomFigures/Custom_GSEA/TControl_vs_NControl_Curated_GO_barplot.tiff", GOdotplot, width = 12, height = 10, dpi = 300)

KEGGdotplot <- plotbar(curKEGG)
ggsave("Outputs/007_TControl_vs_NControl_Outputs/CustomFigures/Custom_GSEA/TControl_Vs_NControl_Curated_KEGG_barplot.tiff", KEGGdotplot, width = 12, height = 10, dpi = 300)








