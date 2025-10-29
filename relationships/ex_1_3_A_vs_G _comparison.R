# ---
# Exercise: Comparing Pedigree (A) and Genomic (G) Relationship Matrices
# ---
#
# Goal:
# 1. Load a real pedigree file and calculate the A matrix.
# 2. Load a real genotype file and calculate the G matrix.
# 3. Align both matrices to compare the same set of animals.
# 4. Compare the matrices by:
#    a) Plotting their off-diagonal elements (relationships).
#    b) Plotting their diagonal elements (inbreeding).
#    c) Calculating the correlation between them.
#
# ---
# Step 0: Install and Load Required Libraries
# ---

# Install packages if you don't have them (uncomment the lines below)
# install.packages("pedigreemm")
# install.packages("rrBLUP")

library(pedigreemm)
library(rrBLUP)

cat("--- Libraries loaded ---\n")

# ---
# Step 1: Load Pedigree Data and Calculate A Matrix
# ---
# We will read the pedigree file directly from the GitHub URL.
# The format is: Animal, Sire, Dam

ped_url <- "https://raw.githubusercontent.com/thembee/GS/refs/heads/main/relationships/pedigree.txt"
ped_df <- read.table(ped_url, col.names = c("Animal", "Sire", "Dam"))

cat(paste("Loaded pedigree for", nrow(ped_df), "animals.\n"))

# Create the pedigree object
# Note: The 'pedigree' function from pedigreemm handles sorting
# and renumbering (if needed) internally.
ped <- pedigree(ped_df$Sire, ped_df$Dam, ped_df$Animal)

# Calculate the full A matrix (Numerator Relationship Matrix)
# This matrix will have dimensions for ALL animals in the pedigree.
A_matrix_full <- as.matrix(getA(ped))

cat(paste("Calculated full A matrix. Dimensions:",
          nrow(A_matrix_full), "x", ncol(A_matrix_full), "\n"))
# print(head(rownames(A_matrix_full))) # See animal IDs

# ---
# Step 2: Load Genotype Data and Calculate G Matrix
# ---
# We will load the genotype file.
# The format is 0, 1, 2 (homozygous, heterozygous, homozygous).
# Animal IDs are the row names.

geno_url <- "https://raw.githubusercontent.com/thembee/GS/refs/heads/main/relationships/genotypes.txt"

# 'header=TRUE' (SNPs are headers), 'row.names=1' (Animal IDs are row names)
geno_df_012 <- read.table(geno_url, header = TRUE, row.names = 1, check.names = FALSE)

cat(paste("Loaded genotypes for", nrow(geno_df_012), "animals and",
          ncol(geno_df_012), "markers.\n"))

# ---
# IMPORTANT: Recode Genotypes for rrBLUP::A.mat
# ---
# The A.mat() function (VanRaden, 2008, Method 1) expects markers
# to be coded as -1, 0, 1.
# Our file is 0, 1, 2. We must subtract 1.
# 0 -> -1
# 1 ->  0
# 2 ->  1

M_matrix_012 <- as.matrix(geno_df_012)
M_matrix_m101 <- M_matrix_012 - 1

cat("Recoded genotypes from (0, 1, 2) to (-1, 0, 1).\n")

# Calculate the G matrix (Genomic Relationship Matrix)
G_matrix <- A.mat(M_matrix_m101)

cat(paste("Calculated G matrix. Dimensions:",
          nrow(G_matrix), "x", ncol(G_matrix), "\n"))
# print(head(rownames(G_matrix))) # See animal IDs

# ---
# Step 3: Align Matrices (CRITICAL STEP!)
# ---
# The A matrix is for *all* 460 animals in the pedigree.
# The G matrix is only for the *genotyped* animals (60 animals).
# We can only compare animals that are in BOTH matrices.

# Get the list of genotyped animals (from G matrix row names)
genotyped_animal_ids <- rownames(G_matrix)

# Subset the full A matrix to *only* include these genotyped animals,
# in the *same order* as they appear in the G matrix.
A_matrix_aligned <- A_matrix_full[genotyped_animal_ids, genotyped_animal_ids]

cat(paste("Aligned A matrix to match G matrix. Final dimensions:",
          nrow(A_matrix_aligned), "x", ncol(A_matrix_aligned), "\n"))

# ---
# Step 4: Compare A and G Matrices
# ---

# We will compare two things:
# 1. Off-diagonals: These represent the relationship/relatedness
#    coefficients between pairs of animals.
# 2. Diagonals: These represent (1 + F), where F is the
#    inbreeding coefficient.

# Extract the upper triangle (off-diagonal elements)
A_offdiag <- A_matrix_aligned[upper.tri(A_matrix_aligned)]
G_offdiag <- G_matrix[upper.tri(G_matrix)]

# Extract the diagonals and calculate inbreeding (F = diagonal - 1)
F_pedigree <- diag(A_matrix_aligned) - 1
F_genomic <- diag(G_matrix) - 1

# ---
# 4a. Calculate Correlations
# ---

cor_relationships <- cor(A_offdiag, G_offdiag)
cor_inbreeding <- cor(F_pedigree, F_genomic)

cat("\n--- Comparison Results ---\n")
cat(paste("Correlation of Off-Diagonal (Relationships):", round(cor_relationships, 4), "\n"))
cat(paste("Correlation of Diagonal (Inbreeding):     ", round(cor_inbreeding, 4), "\n\n"))

cat("Summary of Pedigree Inbreeding (F_ped):\n")
print(summary(F_pedigree))

cat("\nSummary of Genomic Inbreeding (F_gen):\n")
print(summary(F_genomic))

# ---
# 4b. Plot Comparisons
# ---
# Set up a 1x2 plotting window
par(mfrow = c(1, 2), pty = "s") # pty="s" makes plots square

# Plot 1: Off-Diagonal Relationships
plot(A_offdiag, G_offdiag,
     main = "Off-Diagonal Comparison (Relationships)",
     xlab = "Pedigree Relationships (A)",
     ylab = "Genomic Relationships (G)",
     pch = 20,
     col = "#00008B44") # Dark blue with transparency
# Add a 1:1 line
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# Add correlation to the plot
legend("topleft",
       legend = paste("cor =", round(cor_relationships, 3)),
       bty = "n",
       text.col = "red")

# Plot 2: Diagonal Inbreeding
plot(F_pedigree, F_genomic,
     main = "Diagonal Comparison (Inbreeding, F)",
     xlab = "Pedigree Inbreeding (F_ped)",
     ylab = "Genomic Inbreeding (F_gen)",
     pch = 20,
     col = "#8B000044") # Dark red with transparency
# Add a 1:1 line
abline(a = 0, b = 1, col = "blue", lty = 2, lwd = 2)
# Add correlation to the plot
legend("topleft",
       legend = paste("cor =", round(cor_inbreeding, 3)),
       bty = "n",
       text.col = "blue")

# Reset plotting window
par(mfrow = c(1, 1))

cat("\n--- Exercise Complete. Check the 'Plots' window. ---\n")
