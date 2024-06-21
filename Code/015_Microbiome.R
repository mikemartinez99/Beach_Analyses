# The purpose of this script is to analyze the 16S microbiome samples sent to UChicago for control and combination tumors

# Clear the environment
rm(list = ls())

# Load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)
library(vegan)
library(ggpubr)
library(EnhancedVolcano)
library(dada2)
library(metagMisc)
library(microViz)

# Generate output directory
opPath <- "Outputs/015_Microbiome_Outputs"
if (!dir.exists(opPath)) {
  dir.create(opPath)
}

# Generate subdirectory for Rds files
rdsPath <- "Outputs/015_Microbiome_Outputs/Rds_Files"
if (!dir.exists(rdsPath)) {
  dir.create(rdsPath)
}

# Generate subdirectory for figures
figPath <- "Outputs/015_Microbiome_Outputs/Figures"
if (!dir.exists(figPath)) {
  dir.create(figPath)
}

################################################################################

# Read in the phyloseq object provided to us by UChicago
ps <- readRDS("Data_Files/Microbiome/UChicago/MMF.16S.402_RyanBeach.finalPhy_merged_rdp.rds")

# Assign a more meaningful group designation than 1 or 2
ps@sam_data$GroupName <- ifelse(ps@sam_data$group == "02", "Combination", "Control")
ps@sam_data$GroupName <- factor(ps@sam_data$GroupName, levels = c("Control", "Combination"))
meta <- as.data.frame(ps@sam_data)
meta$Site <- rownames(meta)

# What is the minimum and maximum sequencing depth?
min(sample_sums(ps))
max(sample_sums(ps))

################################################################################

# Plot rarefection curves to see if we should rarefy for analysis or not
otu <- otu_table(ps)
otu.df <- as.data.frame(otu)
sample_names <- rownames(otu)

# Calculate rarefaction curves. Tidy = TRUE will return a dataframe for custom plotting
otu.rarecurve <- rarecurve(otu.df, step = 100, label = TRUE, tidy = TRUE)

# Append group information to otu.rarecurve
matched_indices <- match(otu.rarecurve$Site, meta$Site)
otu.rarecurve$Group <- meta$GroupName[matched_indices]
otu.rarecurve$Sample_ID <- otu.rarecurve$Site
otu.rarecurve$Group <- factor(otu.rarecurve$Group, levels = c("Control", "Combination"))

# Plot rarefaction curves
curves <- ggplot(otu.rarecurve, aes(x = Sample, y = Species, color = Group)) +
  theme_bw() + 
  geom_point(size = 0.5) +
  geom_text(data = otu.rarecurve %>%
              group_by(Site) %>%
              filter(Species == max(Species)),
            aes(label = Sample_ID), size = 3.4, hjust = -0.1) +
  scale_x_continuous(breaks = seq(0, max(otu.rarecurve$Sample), by = 5000)) +
  labs(x = "Number of Sequence Reads",
       y = "Number of ASVs",
       title = "Rarefaction Curves") +
  geom_hline(yintercept = max(otu.rarecurve$Species), linetype = "dashed", color = "black") +
  geom_hline(yintercept = mean(otu.rarecurve$Species), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 30582, linetype = "solid", color = "black") +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22, face = "bold"),
        title = element_text(size = 18),
        legend.text = element_text (size = 18))
ggsave("Outputs/015_Microbiome_Outputs/Figures/Grouped_Rarefaction_Curves.png", curves, width = 14, height = 10)

################################################################################

# Set an output directory for alpha diversity
alphaOP <- "Outputs/015_Microbiome_Outputs/Figures/Alpha_Diversity"
if (!dir.exists(alphaOP)) {
  dir.create(alphaOP)
}

# Alpha diversity estimation
richness <- estimate_richness(ps, measures = c("Observed", "Shannon"))

# Set a new column to store Sample names and remove leading "X"
richness$Sample <- rownames(richness)
richness$Sample <- sub("^X", "", richness$Sample)
rownames(richness) <- richness$Sample

# Now, let's match up the GroupName data by the rownames
richness <- merge(richness, meta[,23:25], by = 0)

# Factor GroupName
richness$GroupName <- factor(richness$GroupName, levels = c("Control", "Combination"))

# Plot
nonrare_alphaDiv <- ggplot(richness, aes(x = GroupName, y = Shannon, fill = GroupName)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.format") +
  theme_bw() +
  scale_fill_manual(values = c("Combination" = "firebrick3", "Control" = "skyblue")) +
  labs(y = "Estimate",
       x = "",
       title = "Shannon Index") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        title = element_text(size = 20))
ggsave("Outputs/015_Microbiome_Outputs/Figures/Alpha_Diversity/Non_rarefied_Shannon_Diversity.png", nonrare_alphaDiv)

################################################################################

# Create an output directory for beta diversity
betaOP <- "Outputs/015_Microbiome_Outputs/Figures/Beta_Diversity"
if (!dir.exists(betaOP)) {
  dir.create(betaOP)
}

# Rarefy to an even depth for beta diversity calculations
set.seed(3699)
ps.30k <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = TRUE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

# Save the rarefied phyloseq object
saveRDS(ps.30k, file = "Data_Files/Microbiome/Rarefied_30582_Phyloseq_Object.Rds")

# Calculate PCoA (Bray-Curtis)
pcoa_bc <- ordinate(ps.30k, "PCoA", "bray")

# Plot ordination
betaDiv <- plot_ordination(ps.30k, pcoa_bc, type = "samples", color = "GroupName") +
  geom_point(size = 1) +
  geom_text(aes(label = sampleid), size = 4, hjust = 1.3, vjust = 1) +  
  theme_bw() +
  stat_ellipse(level = 0.99, aes(group = GroupName), linetype = 2) +
  stat_ellipse(level = 0.95, aes(group = GroupName), linetype = 1)
ggsave("Outputs/015_Microbiome_Outputs/Figures/Beta_Diversity/Rarefied_30k_BetaDiversity.tiff", betaDiv, width = 8, height = 8, dpi = 300)

# Get OTU table from rarefied data
rarOTUs = data.frame(otu_table(ps.30k))

# Generate distance matrix
ps.30k.bc.dist <- phyloseq::distance(ps.30k, method = "bray")

# Get metadata from sample_data slot of the phyloseq object
ps.bc.dist.df <- data.frame(sample_data(ps.30k))

# Conduct PERMANOVA test
adonis2(ps.30k.bc.dist~GroupName, data = ps.bc.dist.df)

################################################################################

# Generate directory for community compositions
commDir <- "Outputs/015_Microbiome_Outputs/Community_Compositions"
if (!dir.exists(commDir)) {
  dir.create(commDir)
}

# Generate directory for tabular glommed data
glomData <- "Outputs/015_Microbiome_Outputs/Tabular_Community_Composition_Glomms"
if (!dir.exists(glomData)) {
  dir.create(glomData)
}

# Transform the sample counts to relative abundances
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))

# Save RDS file for non-rarefied relative abundances
saveRDS(ps.rel, file = "Data_Files/Microbiome/Nonrarefied_ps_RelAbund_Transformed_Phyloseq_Object.Rds")

# This function takes a relative abundance transformed phyloseq object and a taxlevel and gloms at that taxlevel, melts to a dataframe and saves
taxLevel.glom <- function(x, rank) {
  glommed <- tax_glom(ps.rel, taxrank = rank, NArm = TRUE)
  glommedData <- psmelt(glommed)
  fileName <-paste(rank, "relative_Abundances.csv", sep = "_")
  write.csv(glommedData, file = paste("Outputs/015_Microbiome_Outputs/Tabular_Community_Composition_Glomms", fileName, sep = "/"))
  return(glommedData)
}

# Get a vector of taxranks
taxRanks <- rank_names(ps.rel)

# Empty list to hold the glommed data
abundances <- list()

# For each taxonomic level, run the taxLevel.glom function
for(taxLevel in taxRanks) {
  level_relAbund <- taxLevel.glom(ps.rel, taxLevel)
  abundances[[taxLevel]] <- level_relAbund
}

################################################################################

# Specify a phylum ouput directory
phylumDir <- paste(commDir, "Phylum", sep = "/")
if (!dir.exists(phylumDir)) {
    dir.create(phylumDir)
}

# Create a custom palette for phylum
customPal <- tax_palette(
  data = ps.rel, rank = "Phylum", pal = "kelly", n = 20, add = c(Other = "black")
)
tax_palette_plot(customPal)

# Classify everything with a relative abundance less than or equal to 0.15 as "other"
phylumData <- abundances[["Phylum"]]
phylumData$Phylum <- as.character(phylumData$Phylum)
PhylumMax <- ddply(phylumData, ~Phylum, function(x) c(Max = max(x$Abundance)))
Other <- PhylumMax[PhylumMax$Max <= 0.15,]$Phylum
phylumData[phylumData$Phylum %in% Other,]$Phylum <- "Other"

# Factor GroupName
phylumData$GroupName <- factor(phylumData$GroupName, levels = c("Control", "Combination"))

# Factor the phylum names so "other" appears last
phylumData$Phylum <- factor(phylumData$Phylum, levels = c("Bacteroidetes", "Firmicutes", "Verrucomicrobia", "Other"))

# Now we should be able to plot our relative abundance barplot
phylumplot <- ggplot(phylumData, aes(x = sampleid, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_grid(~GroupName, scales = "free_x") +
  theme_bw() +
  labs(x = "Sample ID",
       y= "Relative Abundance (%)",
       title = "Phylum-Level Relative Abundance") +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "bottom") +
  scale_fill_manual(values = customPal) +
  scale_y_continuous(labels = scales::percent)
ggsave("Outputs/015_Microbiome_Outputs/Community_Compositions/Phylum/Phylum_Relative_Abundance.tiff", phylumplot, dpi = 200)

################################################################################











