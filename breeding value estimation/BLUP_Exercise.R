# --- Scenario Data Overview ---
#
# This exercise is based on the following 8 phenotype records
# from 7 animals across 3 herds.
#
# This table explains the y, X, and Z matrices below:
#
# ========================================================
# Obs # | Animal ID | Herd ID | Phenotype (y)
# --------------------------------------------------------
#   1   |     4     |    1    |      4.5
#   2   |     5     |    2    |      2.9
#   3   |     6     |    2    |      3.9
#   4   |     7     |    1    |      3.5
#   5   |     8     |    3    |      5.0
#   6   |     9     |    2    |      4.1
#   7   |     10    |    3    |      3.8
#   8   |     4     |    1    |      4.3
# ========================================================
#
# Note: Animal 4 has two records. Animals 1, 2, and 3
# are in the pedigree (see A matrix) but have no phenotype records.
#
# --- BLUP Exercise - Expanded Scenario ---
#
# Objective: Solve the Mixed Model Equations (MME) to get
# Best Linear Unbiased Estimates (BLUE) for fixed effects
# and Best Linear Unbiased Predictions (BLUP) for random effects.
#
# This is an exercise. You will need to fill in the missing code sections
# marked with # <-- FILL IN
#
# Scenario:
# - 8 animal observations (phenotypes)
# - 3 fixed herd effects
# - 10 animals in the pedigree
#
# --- Part 0: Define Data and Model Components ---

# Variance components
sigma2a <- 20 # Additive genetic variance
sigma2e <- 40 # Residual variance

# y - vector of 8 observations (phenotypes)
y <- matrix(c(
  4.5, # Obs 1
  2.9, # Obs 2
  3.9, # Obs 3
  3.5, # Obs 4
  5.0, # Obs 5
  4.1, # Obs 6
  3.8, # Obs 7
  4.3  # Obs 8
), ncol = 1, byrow = TRUE)

# X - incidence matrix for 3 fixed herd effects
# Q: Which observations belong to which herd?
X <- matrix(c(
  1, 0, 0, # Obs 1 -> Herd 1
  0, 1, 0, # Obs 2 -> Herd 2
  0, 1, 0, # Obs 3 -> Herd 2
  1, 0, 0, # Obs 4 -> Herd 1
  0, 0, 1, # Obs 5 -> Herd 3
  0, 1, 0, # Obs 6 -> Herd 2
  0, 0, 1, # Obs 7 -> Herd 3
  1, 0, 0  # Obs 8 -> Herd 1
), ncol = 3, byrow = TRUE)

# Z - incidence matrix for 10 random animal effects
# Q: Which animal has which observation?
# Q: Does any animal have more than one observation?
Z <- matrix(c(
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, # Obs 1 -> Animal 4
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, # Obs 2 -> Animal 5
  0, 0, 0, 0, 0, 1, 0, 0, 0, 0, # Obs 3 -> Animal 6
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, # Obs 4 -> Animal 7
  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, # Obs 5 -> Animal 8
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, # Obs 6 -> Animal 9
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, # Obs 7 -> Animal 10
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0  # Obs 8 -> Animal 4 (repeat)
), ncol = 10, byrow = TRUE)

# --- Part 1: BLUP using the Numerator Relationship Matrix (A) ---

# A - the 10x10 Numerator Relationship Matrix
A <- matrix(c(
  1.00, 0.00, 0.00, 0.50, 0.00, 0.50, 0.00, 0.25, 0.25, 0.37,
  0.00, 1.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.25, 0.25, 0.25,
  0.00, 0.00, 1.00, 0.00, 0.50, 0.00, 0.00, 0.50, 0.25, 0.25,
  0.50, 0.00, 0.00, 1.00, 0.00, 0.25, 0.00, 0.12, 0.25, 0.18,
  0.00, 0.50, 0.50, 0.00, 1.00, 0.25, 0.00, 0.37, 0.37, 0.37,
  0.50, 0.50, 0.00, 0.25, 0.25, 1.00, 0.00, 0.50, 0.37, 0.37,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00,
  0.25, 0.25, 0.50, 0.12, 0.37, 0.50, 0.00, 1.00, 0.31, 0.43,
  0.25, 0.25, 0.25, 0.25, 0.37, 0.37, 0.00, 0.31, 1.00, 0.31,
  0.37, 0.25, 0.25, 0.18, 0.37, 0.37, 0.00, 0.43, 0.31, 1.00
), ncol = 10, byrow = TRUE)

# --- Task 1.1: Calculate A-inverse ---
# Use the solve() function to find the inverse of A
# Ainv <- ... # <-- FILL IN

# Print Ainv to check
print("Inverse of A (Ainv):")
# print(Ainv)

# --- Task 1.2: Calculate alpha ---
# Alpha is the ratio of residual variance to additive variance
# alpha <- ... # <-- FILL IN

print("Alpha value:")
# print(alpha)

# --- Task 1.3: Set up the MME ---
# The MME are:
# | X'X     X'Z     | | b |   | X'y |
# | Z'X     Z'Z+Ainv*alpha | | a | = | Z'y |
#
# Let's build the Left-Hand Side (LHS) and Right-Hand Side (RHS)

# Top-left block: X'X
# Use t() for transpose and %*% for matrix multiplication
# TL <- ... # <-- FILL IN

# Top-right block: X'Z
# TR <- ... # <-- FILL IN

# Bottom-left block: Z'X
# BL <- ... # <-- FILL IN

# Bottom-right block: Z'Z + Ainv*alpha
# BR_A <- ... # <-- FILL IN

# Combine blocks into the full LHS matrix
# Use rbind() to combine rows and cbind() to combine columns
# LHS_A <- ... # <-- FILL IN

# Build the RHS vector
# RHS_A <- ... # <-- FILL IN

print("LHS (A-BLUP):")
# print(LHS_A)
print("RHS (A-BLUP):")
# print(RHS_A)


# --- Task 1.4: Solve the MME for A-BLUP ---
# Solve the system LHS * SOL = RHS for SOL
# SOL_A <- ... # <-- FILL IN

# Print the solutions
# The first part contains the BLUE estimates for fixed effects (b)
# The second part contains the BLUP estimates for random effects (a)
print("Solutions (A-BLUP):")
# print(SOL_A)


# --- Part 2: G-BLUP using the Genomic Relationship Matrix (G) ---

# G - the 10x10 Genomic Relationship Matrix
G <- matrix(c(
  1.00, 0.00, 0.00, 0.60, 0.00, 0.50, 0.02, 0.20, 0.15, 0.40,
  0.00, 1.00, 0.00, 0.00, 0.50, 0.50, 0.01, 0.25, 0.30, 0.22,
  0.00, 0.00, 1.00, 0.00, 0.55, 0.00, 0.02, 0.45, 0.30, 0.30,
  0.60, 0.00, 0.00, 1.00, 0.00, 0.25, 0.02, 0.20, 0.21, 0.20,
  0.00, 0.50, 0.55, 0.00, 1.00, 0.25, 0.02, 0.37, 0.45, 0.40,
  0.50, 0.50, 0.00, 0.25, 0.25, 1.00, 0.01, 0.50, 0.40, 0.35,
  0.02, 0.01, 0.02, 0.02, 0.02, 0.01, 1.00, 0.20, 0.05, 0.03,
  0.20, 0.25, 0.45, 0.20, 0.37, 0.50, 0.20, 1.00, 0.30, 0.44,
  0.15, 0.30, 0.30, 0.21, 0.45, 0.40, 0.05, 0.30, 1.00, 0.28,
  0.40, 0.22, 0.30, 0.20, 0.40, 0.35, 0.03, 0.44, 0.28, 1.00
), ncol = 10, byrow = TRUE)


# --- Task 2.1: Calculate G-inverse ---
# Use the solve() function to find the inverse of G
# Ginv <- ... # <-- FILL IN

print("Inverse of G (Ginv):")
# print(Ginv)


# --- Task 2.2: Set up the MME for G-BLUP ---
# The equations are the same, but we substitute Ginv for Ainv
#
# Q: Why is alpha the same as in Part 1?
# Q: Can we re-use any of the matrices from Part 1?
#    (Hint: X, y, and Z have not changed!)

# Bottom-right block: Z'Z + Ginv*alpha
# BR_G <- ... # <-- FILL IN

# Combine blocks into the full LHS matrix
# (Hint: You can re-use TL, TR, and BL from Part 1)
# LHS_G <- ... # <-- FILL IN

# The RHS is identical to the one from Part 1
RHS_G <- RHS_A # This line is kept as it's an instruction

print("LHS (G-BLUP):")
# print(LHS_G)


# --- Task 2.3: Solve the MME for G-BLUP ---
# Solve the system
# SOL_G <- ... # <-- FILL IN

# Print the solutions
# The first part contains the BLUE estimates for fixed effects (b)
# The second part contains the BLUP estimates for genomic breeding values (g)
print("Solutions (G-BLUP):")
# print(SOL_G)


# --- Part 3: Compare Solutions ---

# --- Task 3.1: Create a comparison table ---
# Combine the solutions into a data frame for easy comparison
# Make sure to label all 13 effects (3 fixed + 10 random)

# comparison <- ... # <-- FILL IN
# (Hint: use data.frame() and c() to create labels)

print("Comparison of A-BLUP and G-BLUP Solutions:")
# print(comparison)

# --- Discussion Questions ---
#
# 1. Compare the fixed effect estimates (Herd_1, Herd_2, Herd_3) from
#    A-BLUP and G-BLUP. Are they different? Why or why not?
#
# 2. Compare the random effect predictions (Animal_1 to Animal_10).
#    Which animals show the largest differences between A-BLUP and G-BLUP?
#
# 3. Look at animals 1, 2, 3, and 7. In the A matrix, they are "unrelated"
#    founders (relationship = 0). What are their relationships in the
#    G matrix? What does this tell you?
#
# 4. Animal 4 has two records. How does the model handle this?
#    (Hint: Look at the Z matrix and the Z'Z matrix)
#
# 5. Why might the genomic-based (G) solutions be preferred over the
#    pedigree-based (A) solutions in a modern breeding program?
#

