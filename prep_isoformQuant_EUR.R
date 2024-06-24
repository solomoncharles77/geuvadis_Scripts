
# # This script takes one argument -  the prefix name of the plink processed output
# aha <- commandArgs(trailingOnly = TRUE)
# 
# aha2 <- paste0("genoFiles/", aha, ".Geno.txt")
# aha3 <- paste0("phenoFiles/", aha, ".txt")

# load libraries ----------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))


# Import and process geno data --------------------------------------------
genoSamps <- read.table("genoFiles/eurGeno_samples.txt")
exprSamps <- read.table("exprFiles/matched_EUR_Pheno.txt")
head(genoSamps)

cat("\n")

cat("This is a view of the expression data \n")
cat("\n")

# Import and process Salmon quantified expr data --------------------------------------------
tx2gene <- read.csv("/scratch/vasccell/cs806/exprPhenoData/tx2gene.ensembl.v100.csv")
files <- file.path("/lustre/alice3/scratch/vasccell/cs806/geuvadis/salmonQuant", exprSamps$V2, "quant.sf")
names(files) <- paste0(exprSamps$V2)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = T)

cat("\n")
cat("Normalizing expression data \n")
cat("\n")

# normalize expr
sampleTable <- cbind(exprSamps, condition = factor("LCL"))

dds <- DESeqDataSetFromTximport(txi.salmon, colData = sampleTable, ~1)
ddsF <- dds[ rowSums(counts(dds)) > ncol(dds)/2, ]
vst=varianceStabilizingTransformation(ddsF)
normExpr <- as.data.frame(assay(vst))

# Rename to match name in vcf file, order samples and Re-add geneID columns
normExpr <- normExpr[, order(names(normExpr))]
normExpr <- cbind(txID = substr(rownames(normExpr), 1, 15), normExpr)

cat("This is a view of the normalized expression data \n")
cat("\n")
normExpr[1:10, 1:10]


# Prep bed for normalized genes -----------------------------------
# Add coordinate information 
exprCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Transcript_Coords_withStrand.txt"))
colnames(exprCoord)[1] <- "txID"
normExprCoord <- exprCoord[which(exprCoord$txID %in% normExpr$txID), ]
normExprCoord <- merge(normExprCoord, normExpr, by = "txID")
normExprCoord <- normExprCoord[, c(3:5,1:2,6, 8:ncol(normExprCoord))]
colnames(normExprCoord)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
normExprCoord <- normExprCoord[order(normExprCoord[,1], normExprCoord[,2]), ]
normExprCoord[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
exprFileName <- paste0("phenoFiles/isoformQuant_EUR_qtlTools_Ready.bed")
fwrite(normExprCoord, file = exprFileName, sep = "\t")

# system(paste0("module load samtools"))
# system(paste0("module load tabix"))
system(paste0("bgzip -f ", exprFileName))
system(paste0("tabix ", exprFileName, ".gz"))
