
# # This script takes one argument -  the prefix name of the plink processed output
# aha <- commandArgs(trailingOnly = TRUE)
# 
# aha2 <- paste0("genoFiles/", aha, ".Geno.txt")
# aha3 <- paste0("phenoFiles/", aha, ".txt")

# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))


# Import and combine phenodata data --------------------------------------------
isoformQuant <- data.frame(fread("phenoFiles/isoformQuant_EUR_qtlTools_Ready.bed.gz"))
leafcutter <- data.frame(fread("phenoFiles/geneLeafcutter_EUR_qtlTools_Ready.bed.gz"))
exon <- data.frame(fread("phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed.gz"))
intron <- data.frame(fread("phenoFiles/geneIntronCount_EUR_qtlTools_Ready.bed.gz"))
exonIntronRatio <- data.frame(fread("phenoFiles/geneExonIntronRatio_EUR_qtlTools_Ready.bed.gz"))

allPhe <- rbind(isoformQuant, leafcutter, exon, intron, exonIntronRatio)
allPhe <- allPhe[order(allPhe$X.chr, allPhe$start), ]
colnames(allPhe)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")

# Export files for QTLtools ---------------------------------------------
phenoFileName <- paste0("phenoFiles/combinedPheno_EUR_qtlTools_Ready.bed")
fwrite(allPhe, file = phenoFileName, sep = "\t")
system(paste0("bgzip -f ", phenoFileName))
system(paste0("tabix ", phenoFileName, ".gz"))

