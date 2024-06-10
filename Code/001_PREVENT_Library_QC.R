# The purpose of this script is to generate a cohesive counts matrix based for the PREVENT analysis

# Load libraries
library(dplyr)
library(tidyverse)
library(DESeq2)
library(genefilter)
library(ggplot2)
library(ggh4x)
library(grDevices)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(RColorBrewer)

# Set the path to the raw counts files from htseq
countsDir <- "Data_Files/Raw_Counts_Files/"

# Make a new directory to hold the master counts files
dirPath <- "Data_Files/MASTER_COUNTS/"
if (!dir.exists(dirPath)) {
  dir.create(dir_path)
} 

# Create an output directory to hold the outputs of this script
opPath <- "Outputs/001_PREVENT_Library_QC_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Combine the files into a cohesive counts matrix
# Ensure you are running this from the Rproj file so relative file paths work
files <- list.files("Data_Files//Raw_Counts_Files/", pattern = "*.counts$")
files
i<-0
while (i<length(files)){
  if (i == 0){
    print(paste0("Read in file : ",files[i+1]))
    fileName=paste0(countsDir,"/",files[i+1])
    df1<-read.table(fileName,sep="\t",stringsAsFactors = FALSE,header=F)
    colnames(df1)<-c("geneID",strsplit(files[i+1],".tsv"))
    i<-i+1
  }else{
    print(paste0("Read in file : ",files[i+1]))
    fileName=paste0(countsDir,"/",files[i+1])
    df2<-read.table(fileName,sep="\t",stringsAsFactors = FALSE,header=F)
    colnames(df2)<-c("geneID",strsplit(files[i+1],".tsv"))
    df1<-merge(df1,df2,by.x="geneID",by.y="geneID",sort=FALSE)
    i<-i+1
  }  
}
head(df1)
dim(df1)
df1<-df1[(!stringr::str_starts(df1[["geneID"]],"__")),]
dim(df1)

# Remove .counts from the colnames
rownames(df1) <- df1$geneID
df1$geneID <- NULL
new_colnames <- sub("\\.counts$", "", colnames(df1))
colnames(df1) <- new_colnames

# Read in the metadata and get samples coordinated
meta <- read.csv("Data_Files/PREVENT_METADATA.csv", header = TRUE, sep = ",")
meta <- meta[-15,]
rownames(meta) <- meta$Sample
colnames(df1)[!(colnames(df1) %in% rownames(meta))]
write.csv(meta, file = "Data_Files/PREVENT_METADATA.csv")

# Let's get the rownames and column names in the same order
sampleOrder <- rownames(meta)
df1 <- df1[,sampleOrder]

# Check that the samples are in the same order between the counts and the metadata
all(colnames(df1) == rownames(meta))

# Prepare the gene symbols
df1$Ensembl <- rownames(df1)
df1$Symbol <- mapIds(org.Rn.eg.db, key = df1$Ensembl, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first")
rownames(df1) <- paste(df1$Ensembl, df1$Symbol, sep = " - ")
df1$Ensembl <- NULL
df1$Symbol <- NULL

# Write as a csv
write.csv(df1, file = "Data_Files/MASTER_COUNTS/Prevent_Master_Counts.csv")

################################################################################

# Assess the library sizes graphically
df.m <- reshape2::melt(df1, id.vars =NULL)
colnames(df.m) <- c("Sample", "Value")
df.m <- merge(df.m, meta, by = "Sample")

df.m$Group <- factor(df.m$Group, levels = c("WT", "Control", "Naproxen", "EPA", "Combo", "TP252", "TP252_Naproxen", "Ifetroban"))
groups <- unique(df.m$Group)
colors <- setNames(brewer.pal(8, "Paired"), groups)


# Plot
QC <- ggplot(df.m,aes(factor(Sample),log10(Value), fill=Group)) +
  geom_violin() +
  facet_nested(~Tissue + Group, scale = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        strip.text= element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = colors) +
  labs(x = "Sample",
       y = "log10 Counts",
       title = "Library QC")
QC
ggsave("Outputs/001_PREVENT_Library_QC_Outputs/Prevent_All_Samples_Library_Size_Violin.tiff", QC, width = 20, height = 12, dpi = 100)









