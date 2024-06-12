# The purpose of this script is to separate the tumor and normal counts and to generate an exploratory PCA

# Clear environment
rm(list = ls())

# Load libraries
# Load libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(extrafont)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(dendextend)


# Create an output directory
opPath <- "Outputs/002_Exploratory_PCA_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Read in the master counts and metadata file
counts <- read.csv("Data_Files/MASTER_COUNTS/PREVENT_MASTER_COUNTS.csv")
rownames(counts) <- counts$X
counts$X <- NULL
meta <- read.csv("Data_Files/PREVENT_METADATA.csv")
rownames(meta) <- meta$Sample

# We already checked that everything was in the same order in `001_RNASeq_LibraryQC.R` but we will check again 
colnames(counts) <- sub("\\.", "-", colnames(counts))
sampleOrder <- rownames(meta)
counts <- counts[,sampleOrder]
colnames(counts)[!(colnames(counts) %in% rownames(meta))]
all(colnames(counts) == rownames(meta))

# Let's look at just the tumor samples
tumorSamples <- meta[meta$Tissue == "Tumor",]$Sample
tumorCounts <- counts[,tumorSamples]
tumorMeta <- meta[tumorSamples,]
colnames(tumorCounts)[!(colnames(tumorCounts) %in% rownames(tumorMeta))]
all(colnames(tumorCounts) == rownames(tumorMeta))

# Filter out other samples
remove <- c("T5-3190", "T4-2547", "T21-2852", 
            "T25-4406", "T26-2844", "T35", "T30-4249")

# Isolate the tumor counts
tumorCounts <- tumorCounts[,!(colnames(tumorCounts) %in% remove)]
tumorMeta <- tumorMeta[!rownames(tumorMeta) %in% remove,]

# Generate a directory to hold the filtered tumor counts
tumorDir <- "Data_Files/Filtered_Tumor_Counts"
if (!dir.exists(tumorDir)) {
  dir.create(tumorDir)
}

# Save the tumor counts
write.csv(tumorCounts, file = "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Counts.csv")
write.csv(tumorMeta, file = "Data_Files/Filtered_Tumor_Counts/Final_Filtered_Tumor_Metadata.csv")

# Set up directory architecture
rdsDir <- "Outputs/002_Exploratory_PCA_Outputs/Rds_files"
if (!dir.exists(rdsDir)) {
  dir.create(rdsDir)
}

# For PREVENT analyses
keep <- c("Control", "Naproxen", "EPA", "Combo")
tumorMeta <- tumorMeta[tumorMeta$Group %in% keep, ]

# Filter the tumor counts based on the meta
tumorCounts <- tumorCounts[,colnames(tumorCounts) %in% rownames(tumorMeta)]

# Prepare the DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = tumorCounts,
                              colData = tumorMeta,
                              design = ~ Group)

# Set a reference level (in this case, it doesn't really matter)
dds$Group <- relevel(dds$Group, ref = "Control")

# Pre-filter for low counts
smallestGroup <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroup
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Save the .Rds obejct
saveRDS(dds, file = "Outputs/002_Exploratory_PCA_Outputs/Rds_files/Filtered_TumorSamples_ExploratoryPCA_DESeq2_Object.Rds")

# Run this line if you are re-generating a pca
# dds <- readRDS(file = "Outputs/002_Exploratory_PCA_Outputs/Rds_files/Filtered_TumorSamples_ExploratoryPCA_DESeq2_Object.Rds")

#Create a summarized experiment
vsd <- vst(dds)

# Run PCA
PCA <- plotPCA(vsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(PCA, "percentVar"))
PCA$Group <- factor(PCA$Group, levels = c("Control", "Naproxen", "EPA", "Combo"))

# Set custom colors
custom_colors <- c("Control" = "#66C2A5", "Naproxen" = "#FC8D62", "EPA" = "#8DA0CB", "Combo" = "#E78AC3")

# Plot
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
ggsave("Outputs/002_Exploratory_PCA_Outputs/Final_Filtered_TumorSamples_PCA.tiff", plot, width = 10, height = 10, dpi = 300)


################################################################################

# Clear and re-run for normal
rm(list = ls())

# Read in the master counts and metadata file
counts <- read.csv("Data_Files/MASTER_COUNTS/PREVENT_MASTER_COUNTS.csv")
rownames(counts) <- counts$X
counts$X <- NULL
meta <- read.csv("Data_Files/PREVENT_METADATA.csv")
rownames(meta) <- meta$Sample

# We already checked that everything was in the same order in `001_RNASeq_LibraryQC.R` but we will check again 
colnames(counts) <- sub("\\.", "-", colnames(counts))
sampleOrder <- rownames(meta)
counts <- counts[,sampleOrder]
colnames(counts)[!(colnames(counts) %in% rownames(meta))]
all(colnames(counts) == rownames(meta))

# Let's look at just the normal samples
normalSamples <- meta[meta$Tissue == "Normal",]$Sample
normalCounts <- counts[,normalSamples]
normalMeta <- meta[normalSamples,]
colnames(normalCounts)[!(colnames(normalCounts) %in% rownames(normalMeta))]
all(colnames(normalCounts) == rownames(normalMeta))

# Filter samples
normalremove <- c("W1-4248")
normalCounts <- normalCounts[,!(colnames(normalCounts) %in% normalremove)]
normalMeta <- normalMeta[!rownames(normalMeta) %in% normalremove,]

# Generate a directory to hold the filtered tumor counts
normalDir <- "Data_Files/Filtered_Normal_Counts"
if (!dir.exists(normalDir)) {
  dir.create(normalDir)
}

# Save the normal counts
write.csv(normalCounts, file = "Data_Files/Filtered_Normal_Counts/Final_Filtered_Normal_Counts.csv")
write.csv(normalMeta, file = "Data_Files/Filtered_Normal_Counts/Final_Filtered_Normal_Metadata.csv")

# Has these if you want to include TP252 and TP252 Naproxen combo
normalMeta <- normalMeta[normalMeta$Group != "TP252",]
normalMeta <- normalMeta[normalMeta$Group != "TP252_Naproxen",]
normalCounts <- normalCounts[,colnames(normalCounts) %in% rownames(normalMeta)]
colnames(normalCounts) == rownames(normalMeta)

# Prepare the DESeq2 dataset object
normaldds <- DESeqDataSetFromMatrix(countData = normalCounts,
                                    colData = normalMeta,
                                    design = ~ Group)

# Set a reference level (in this case, it doesn't really matter)
normaldds$Group <- relevel(normaldds$Group, ref = "WT")

# Pre-filter for low counts
smallestGroup <- 6
keep <- rowSums(counts(normaldds) >= 10) >= smallestGroup
normaldds <- normaldds[keep,]

# Run DESeq2
normaldds <- DESeq(normaldds)

# Save the .Rds obejct
saveRDS(normaldds, file = "Outputs/002_Exploratory_PCA_Outputs/Rds_files/Filtered_NormalSamples_ExploratoryPCA_DESeq2_Object.Rds")

# Unhash this line out if you are regenerating figures
normaldds <- readRDS("Outputs/002_Exploratory_PCA_Outputs/Rds_files/Filtered_NormalSamples_ExploratoryPCA_DESeq2_Object.Rds")

# Variance-stabilization transformation
normalvsd <- vst(normaldds)

# Run PCA and return data for custom plotting
normalPCA <- plotPCA(normalvsd, intgroup = "Group", returnData = TRUE) 
percentVar <- round(100* attr(normalPCA, "percentVar"))
normalPCA$Group <- factor(normalPCA$Group, levels = c("WT", "Control", "Naproxen", "EPA", "Combo"))

# Set custom colors
custom_colors <- c("WT" = "#E41A1C", "Control" = "#66C2A5", "Naproxen" = "#FC8D62", "EPA" = "#8DA0CB", "Combo" = "#E78AC3")


normalplot <- ggplot(normalPCA, aes(PC1, PC2, fill = Group, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.10, aes(fill = group)) +
  geom_text(size = 3, aes(label = name), hjust = 1, vjust = 1.5) +
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
ggsave("Outputs/002_Exploratory_PCA_Outputs/Final_Filtered_NormalSamples_PCA.tiff", normalplot, width = 10, height = 10, dpi = 200)





