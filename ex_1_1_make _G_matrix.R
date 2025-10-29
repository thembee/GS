##################################################
### Create G matrix according to VanRaden 2008 ###
###                                            ###
###         Marcin Pszczola 2025               ###
##################################################

## This script can be used to create G matrix according to VanRaden 2008
## The code contains TWO mistakes that will prevent it from running
## or will produce an incorrect G matrix.
##
## Your task: Find and fix both mistakes.
## Hint 1: One bug is in 'Step 5' and will cause an error.
## Hint 2: The second bug is in 'Step 6' and is a formula error.


# --- 1. Initial Setup ---

n=2 # number of animals
m=3 # number of markers

# --- 2. Allele Matrix ---
# This matrix represents the raw alleles for each animal.
# Each animal has two rows (one for each homologous chromosome).
# Allele '1' is coded as 1, Allele '0' (or '2') is coded as 0.
# Animal 1's alleles: (1, 1, 0) and (0, 1, 1)
# Animal 2's alleles: (1, 0, 0) and (0, 1, 0)
alleles <- matrix(nrow=n*2,ncol=m,byrow=TRUE,c(
1, 1, 0,
0, 1, 1,
1, 0, 0, 
0, 1, 0 ))

# --- 3. Genotype Matrix (M) ---
# Goal: Convert the 'alleles' matrix into a genotype 'M' matrix.
# The M matrix will have one row per animal.
# Genotypes are coded as -1 (aa), 0 (Aa), 1 (AA) [where 'A' is allele 1]

M<-NULL # Empty matrix declaration to store results

# Create a sequence to iterate over the 'alleles' matrix,
# taking two rows (one animal) at a time.
# s will be c(1, 3)
s<-seq(1,(n*2),2) 

for(i in 1:length(s)) {
  # For each animal i:
  # 1. Get its two allele rows: alleles[ s[i]:(s[i]+1) , ]
  # 2. Sum the columns: colSums(...) gives the count of allele '1' (0, 1, or 2)
  # 3. Subtract 1: This converts the (0, 1, 2) count to the (-1, 0, 1) genotype coding
  # 4. Bind this new genotype row to the M matrix
  M <- rbind(
       M,colSums( alleles[ s[i]:(s[i]+1) , ])-1 )
  }
# After the loop, M is an n x m matrix of genotypes:
# Animal 1: (0, 1, 0)
# Animal 2: (0, 0, -1)

# --- 4. Exploratory Calculations (Not used for final G) ---
# These lines are likely for demonstration or other calculations,
# they are not part of the main G matrix calculation below.
n.hom.loci<-M%*%t(M) # n x n matrix: dot product of animal genotypes
n.hom.ind<-t(M)%*%M  # m x m matrix: dot product of marker genotypes

# --- 5. Z Matrix (Centered Genotypes) ---
# Goal: Create matrix Z, where genotypes in M are centered by allele frequencies.

# 1. Calculate allele frequencies (p) for allele '1'
# colSums(alleles) = total count of allele '1' at each marker (e.g., 2, 3, 1)
# (n*2) = total number of alleles in the sample (4)
# freq (p) = (0.5, 0.75, 0.25)
freq<-colSums(alleles)/(n*2)

# 2. Calculate the P vector (correction factor for each marker)
# This is 2 * (p - 0.5) as per VanRaden
# P = (0, 0.5, -0.5)
P<-2*(freq-0.5)

# 3. Create the Z matrix by subtracting P from each row of M
#
# *** BUG 1 IS IN THIS LOOP ***
# Hint: You can't assign a value to a row in a matrix that
#       doesn't exist yet. How do you create an empty matrix?
#
for(i in 1:n){
   Z[i,]<-M[i,]-P
   }

# --- 6. Final G Matrix Calculation ---

# 1. Calculate the Numerator: Z * Z'
# (This line will fail until BUG 1 is fixed)
G.numerator<-Z%*%t(Z)

# 2. Calculate the Denominator
#
#
G.denominator<-sum(freq*(1-freq))

# 3. Calculate G
G<- G.numerator / G.denominator

