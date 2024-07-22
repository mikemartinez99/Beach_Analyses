# Introduction 
Colorectal cancer (CRC) is a leading cause of cancer-associated mortality. Our recent study in Pirc rats revealed that the combination of eicosapentaenoic acid and sodium naproxen exhibited additive tumor protection. While the results from this preclinical study were impressive, showing tumor suppression in excess of 95%, potential mechanisms underlying this tumor protection remain unclear. The purpose of this study was to determine the transcriptomic changes resulting in the tumor suppression afforded by the combination of EPA and naproxen supplementation in the diet.

# Methods
Six week-old Pirc rats were fed either a control AIN-93G diet, or AIN-93G diet supplemented with 2% EPA combined with 200 ppm naproxen, for a total of 20 weeks. Global transcriptomic analysis was performed using RNA sequencing (RNA-Seq) of colon tumors and tumor-adjacent normal tissue from 26-week-old Pirc rats. 100-bp paired-end sequencing was performed on an Illumina NovaSeq to a depth of 25 million reads per sample. Reads were trimmed with Trimmomatic-v0.39, and reads were mapped to the mRatBN7.2.105 rat reference genome using HISAT2. Differentially expressed genes were calculated using DESeq2, and gene set enrichment analysis (GSEA) were performed using RStudio.
Results: DESeq2 identified 775 differentially expressed genes (DEGs, 2-fold change > 1, FDR < 0.05) between large tumors from AIN-93G fed rats versus suppressed tumors from rats fed the combination of EPA and naproxen (500 up-regulated, 275 down-regulated). GSEA analysis revealed activation of inflammatory pathways including myeloid cell differentiation, cytokine signaling, TNF-alpha pathway, and IL-17 signaling were enriched in untreated, control tumors. In contrast, suppressed tumors from animals supplemented with EPA and naproxen were enriched in pathways including TGF-b signaling, cell fate specification, PPAR signaling, and protein digestion and absorption. Furthermore, suppressed tumors shared significant GSEA pathway overlap with normal tissue, leading us to believe the combination of EPA and naproxen may be preventing tumorigenesis in part by maintaining cells in a differentiated state.

# Conclusions 
Inflammation plays a key role in tumorigenesis and promoting cancer stem cell niches in tumors. According to our analysis, the combination of EPA and naproxen attenuated cancer-associated inflammation and reduced proliferation in part by augmenting iCAF formation. Future studies will utilize single-nuclei RNA-Seq (snRNA-Seq) to determine at a single cell resolution the transcriptomic changes of epithelial and stromal populations within the TME. This work was supported by National Cancer Institute task order 75N91019F00132.

# Code

# 001_PREVENT_Library_QC.R
Contains the code to collate htseq-generated counts tables into a cohesive counts file and perform basic library QC metrics. Data include samples derived from wild-type (WT) rats, normal and tumor tissue from Pirc rats, normal and tumor tissue from treated rats (treatment groups include the following: Naproxen, EPA, combination of Naproxen + EPA, TP252, Combination of TP252 + Naproxen, and Ifetroban.)

# 002_Exploratory_PCA.R
Code used to generate an exploratory principal component analysis of all groups

# 003_NControl_vs_WT.R
Analyses to compare the gene expression signatures of normal control Pirc tissue relative to WT rats

# 004 - 006
The following set of comparisons are all from normal pirc rat colon tissue.
Script 004_NNaproxen_vs_NControl compares naproxen treated Pirc rat relative-to non-treated Pirc rat
Script 005_NEPA_vs_NControl compares EPA treated Pirc rat relative-to non-treated pirc rat
Script 006_NCombo_vs_NControl compares combination (Naproxen + EPA) treated Pirc rat relative-to non-treated Pirc rat

# 007-012
The following set of comparisons are all from tumor pirc rat colon tissues (except for 007 which is tumor vs normal)
Script 007_TControl_vs_NControl compares non-treated Pirc rat tumor tissue relative to non-treated Pirc rat normal tissue.
Script 008_TNaproxen_vs_TControl compares naproxen-treated Pirc rat tumor tissue relative-to non-treated Pirc rat tumor tissue.
Script 009_TEPA_vs_TControl compares EPA-treated Pirc rat tumor tissue relative-to non-treated Pirc rat tumor tissue.
Script 010_TCombo_vs_TControl compares combination-treated (Naproxen + EPA) Pirc rat tumor tissue relative-to non-treated Pirc rat tumor tissue.
Script 011_TTP252_vs_TControl compares TP252-treated Pirc rat tumor tissue relative-to non-treated Pirc rat tumor tissue.
Script 012_TTP252_Naproxen_vs_TControl compares TP252 + Naproxen-treated Pirc rat tumor tissue relative-to non-treated Pirc rat tumor tissue.

# 013_GSEA_Cluster_Comparison.R
This script makes use of clusterProfilers cluster comparison methods to examine changes in KEGG terms across treatment groups.

# 014_Ifetroban_vs_All.R
This comparison is for the Ifetroban project and compares Ifetroban-treated Pirc Rat tumor tissue relative to non-treated Pirc rat tumor tissue.





