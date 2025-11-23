
# This will install the main Bioconductor installer
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# This uses BiocManager to install the 'biomaRt' package
BiocManager::install("biomaRt")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Bt.eg.db")

BiocManager::install("enrichplot")

