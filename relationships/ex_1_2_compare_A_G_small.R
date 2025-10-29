# R Exercise: Comparing Relationship Matrices A and G 

# Objective
# The goal of this exercise is to load two text files (`A.txt` and `G.txt`) as relationship matrices into R.
# You will then validate these matrices and perform a comparison to understand the differences
# and similarities between them. This is a common task when comparing, for example,
# a genomic relationship matrix (G) with a pedigree-based relationship matrix (A).
#
# This exercise is designed to be beginner-friendly. Most of the code is provided;
# you will just need to fill in the missing pieces (marked with `...`).

# Prerequisites
# You will need the `tidyverse` package collection, which includes `dplyr` (for data manipulation)
# and `ggplot2` (for plotting).

# Installation
# Run this line in your R console if you don't have the `tidyverse` package installed:
# install.packages("tidyverse")


# Part 1: Load Libraries and Data
# ---------------------------------
# First, let's load the packages and read our two matrix files directly from their URLs.
# Since these files have no headers or row names, we will use `header = FALSE`.

# Load the necessary libraries
library(tidyverse)

# Define the URLs for the data (as strings)
# (Fixed the URLs to be plain strings)
url_A <- "https://raw.githubusercontent.com/thembee/GS/refs/heads/main/relationships/A.txt"
url_G <- "https://raw.githubusercontent.com/thembee/GS/refs/heads/main/relationships/G.txt"

# Read the text files into data frames
# header = FALSE tells R the first row is *not* a header
df_A <- read.table(url_A, header = FALSE)
df_G <- read.table(url_G, header = FALSE)

# Convert the data frames to matrices for numerical operations
matrix_A <- as.matrix(df_A)
matrix_G <- as.matrix(df_G)

# Let's inspect the data
cat("Dimensions of A:\n")
print(dim(matrix_A))
cat("\nTop-left corner of A:\n")
print(matrix_A[1:5, 1:5])

cat("\nDimensions of G:\n")
print(dim(matrix_G))
cat("\nTop-left corner of G:\n")
print(matrix_G[1:5, 1:5])


# Part 2: Matrix Validation (Student Task)
# ---------------------------------
# To compare matrices element-by-element, they MUST have the exact same dimensions.
# Since there are no row or column names, we must also *assume* that the
# individuals in row 1 of matrix A are the same as in row 1 of matrix G, and so on.

# --- Task 1: Check Dimensions ---

# TO-DO: Run this check in your console.
# Are the dimensions identical? (all() checks if every comparison is TRUE)
# This is now the most important check!

all(dim(matrix_A) == dim(matrix_G))

# If this check is FALSE, you cannot proceed with the rest of the exercise.
# For this exercise, we will assume it is TRUE.


# Part 3: Compare Matrices (Student Task)
# ---------------------------------
# Now that we have confirmed our matrices have the same dimensions (`matrix_A` and `matrix_G`),
# we can compare them. We *assume* they are in the same order.
#
# We are typically interested in the off-diagonal elements (the relationships *between*
# individuals), as the diagonal represents self-relationship. We will extract
# the lower triangle of each matrix.

# --- Task 1: Create a Data Frame of Relationships ---

# We will pull out the relationship values and put them in a simple data frame.

# TO-DO: Extract the lower triangle elements from each matrix
# lower.tri() gives us a TRUE/FALSE map of the lower triangle.
# Fill in the ... with the correct matrix name.
A_values <- matrix_A[lower.tri(...)]
G_values <- ...[lower.tri(matrix_G)]

# TO-DO: Create a data frame
# Fill in the ... with the variables we just created.
compare_df <- data.frame(
  A = ...,
  G = ...
)

# TO-DO: Add a difference column
# We use mutate() from dplyr to add a new column.
# Fill in the ... to calculate the difference (A minus G).
compare_df <- compare_df %>%
  mutate(Difference = ... - ...)

# Inspect the resulting data (no to-do here, just run this)
head(compare_df)
summary(compare_df)


# --- Task 2: Visualize the Differences ---

# A histogram of the differences (`A - G`) is a great way to see if one
# matrix is systematically higher or lower than the other.

# TO-DO: Run this code block to generate the histogram.
# Read the comments to understand what each line does!

ggplot(compare_df, aes(x = Difference)) +
  # Create the histogram bars
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8, color = "black") +
  # Add a red dashed line at zero (no difference)
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  # Add labels
  labs(
    title = "Histogram of Differences (A - G)",
    x = "Difference in Relationship (A - G)",
    y = "Frequency"
  ) +
  # Use a clean, simple theme
  theme_minimal()


# --- Task 3: Visualize the Correlation ---

# A scatter plot of A vs. G is the best way to see how well the
# relationship values correlate.

# TO-DO: Run this code block to generate the scatter plot.
# We'll plot a random sample if the data is too large, to speed things up
sample_df <- compare_df %>%
  sample_n(min(50000, nrow(compare_df))) # Plot max of 50k points

ggplot(sample_df, aes(x = A, y = G)) +
  # Add the points. alpha = 0.1 makes them transparent to help see density.
  geom_point(alpha = 0.1) +
  # Add the y=x line. If points are on this line, A and G are equal.
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  # Add labels
  labs(
    title = "Comparison of A vs. G Relationship Values",
    subtitle = "Red line is y=x",
    x = "A Matrix Values (Pedigree)",
    y = "G Matrix Values (Genomic)"
  ) +
  theme_minimal()


# Part 4: Interpretation (Student Task)
# ---------------------------------
# Based on your summary statistics and plots, write a one-paragraph conclusion
# (as a comment below).

# * From the Histogram: Look at your plot. Are the bars centered around the
#   red line (zero)? Or are they mostly to one side? What does this tell you
#   about the average difference?
#
# * From the Scatter Plot: Look at your plot. Do the points form a tight line
#   along the red `y=x` line? Or are they spread out (like a wide cloud)?
#   What does this tell you about how well the A matrix values match the GOOD Matrix values?
#
# * From `summary(compare_df)`: Look at the `Min.`, `Mean.`, and `Max.` for
#   the `Difference` column. Does this confirm what you saw in the histogram?
#
#   * Hint: You can also run `cor(compare_df$A, compare_df$G)` in your console
#     to get the exact correlation coefficient (a number between -1 and 1).
#     A number close to 1 means they are very strongly related.

# (Write your answer here as a comment)
