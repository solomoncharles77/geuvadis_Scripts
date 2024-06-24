
# # This script takes one argument -  the prefix name of the plink processed output
# aha <- commandArgs(trailingOnly = TRUE)
# 
# aha2 <- paste0("genoFiles/", aha, ".Geno.txt")
# aha3 <- paste0("phenoFiles/", aha, ".txt")

# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))


# Import and combine phenodata data --------------------------------------------
expr <- data.frame(fread("phenoFiles/geneExpr_qtlTools_Ready.bed"))
lcc <- data.frame(fread("phenoFiles/geuYri_LCsplice_qtlTools_Ready.bed.gz"))

allPhe <- rbind(expr, lcc)
allPhe <- allPhe[order(allPhe$X.chr, allPhe$start), ]
colnames(allPhe)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")

# Export files for QTLtools ---------------------------------------------
phenoFileName <- paste0("phenoFiles/isorform_LCC_combinedPheno_qtlTools_Ready.bed")
fwrite(allPhe, file = phenoFileName, sep = "\t")
system(paste0("bgzip -f ", phenoFileName))
system(paste0("tabix ", phenoFileName, ".gz"))



# Import and combine Covariate data --------------------------------------------

exprCov <- read.table("covFiles/yriGeno_Sex_10pc.txt", header = T)
lccCov <- read.table("covFiles/geuYri_Sex_10GenoPC_10SplicePC.txt", header = T)

# system(paste0("module load samtools"))
# system(paste0("module load tabix"))
# system(paste0("bgzip ", exprFileName))
# system(paste0("tabix ", exprFileName, ".gz"))
