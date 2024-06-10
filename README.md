# Beach_Analyses
Analyses for NIH PREVENT and Ifetroban projects

# 001_PREVENT_Library_QC.R
Contains the code to collate htseq-generated counts tables into a cohesive counts file and perform basic library QC metrics. Data include samples derived from wild-type (WT) rats, normal and tumor tissue from Pirc rats, normal and tumor tissue from treated rats (treatment groups include the following: Naproxen, EPA, combination of Naproxen + EPA, TP252, Combination of TP252 + Naproxen, and Ifetroban.)

# 002_Exploratory_PCA.R
Code used to generate a PCA plot of all tumor samples and a PCA plot of all normal samples

# 003_NControl_vs_WT.R
Analyses to compare normal control Pirc tissue relative to WT rats

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





