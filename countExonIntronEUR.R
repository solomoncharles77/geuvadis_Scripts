suppressMessages(library(Rsubread))
suppressMessages(library(DESeq2))
library(data.table)


# Exon count -----------------------------------------

samFiles <- readLines("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSam.txt")

exonCount <- featureCounts(files = paste0("/lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/", samFiles),
                           annot.ext = "/scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.100_exons.gtf",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "exon",
                           GTF.attrType = "gene_id",
                           strandSpecific = 0, isPairedEnd = TRUE,
                           nthreads = 28)

exonCountDF <- as.data.frame(exonCount$count)
colnames(exonCountDF) <- gsub("Aligned.out.bam", "", colnames(exonCountDF))

# export exonCount results
saveRDS(exonCount, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawExonCount_EUR.rds")

# Normalize with DESeq2 -------------------------------------------------
# create design matrix for DESeq2
exonCount <- readRDS("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawExonCount_EUR.rds")
exonCountDF <- as.data.frame(exonCount$counts)
colnames(exonCountDF) <- gsub("Aligned.out.bam", "", colnames(exonCountDF))

designMat <- data.frame(group = factor(rep("EUR_LCL",
                                           length(colnames(exonCountDF)))))
rownames(designMat) <- colnames(exonCountDF)
all(rownames(designMat) == colnames(exonCountDF))
# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = exonCountDF,
                              colData = designMat, design = ~ 1)
# filter genes whose sum of expression less than one across all samples
ddsF <- dds[ rowSums(counts(dds)) > length(colnames(exonCountDF)), ]

vst=varianceStabilizingTransformation(ddsF)
normExon <- as.data.frame(assay(vst))

# Rename to match name in vcf file, order samples and Re-add geneID columns
normExon <- normExon[, order(names(normExon))]
normExon <- cbind(gene_id = rownames(normExon), normExon)

cat("This is a view of the normalized expression data \n")
cat("\n")
normExon[1:10, 1:10]
write.csv(normExon, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normExonCount_EUR.csv", row.names = F)

# Prep bed for normalized genes -----------------------------------
# Add coordinate information 
exonCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
normExonCoord <- exonCoord[which(exonCoord$gene_id %in% normExon$gene_id), ]
normExonCoord <- merge(normExonCoord, normExon, by = "gene_id")
normExonCoord <- normExonCoord[, c(2:4,1,1,5, 7:ncol(normExonCoord))]
colnames(normExonCoord)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
normExonCoord$feature <- paste0(normExonCoord$feature, "_ExonCount")
normExonCoord <- normExonCoord[order(normExonCoord[,1], normExonCoord[,2]), ]
normExonCoord[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
exonFileName <- paste0("phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed")
fwrite(normExonCoord, file = exonFileName, sep = "\t")

system(paste0("bgzip -f ", exonFileName))
system(paste0("tabix ", exonFileName, ".gz"))
rm(list = ls())


##################################################################################################################
# Intron count -----------------------------------------

samFiles <- readLines("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSam.txt")

intronCount <- featureCounts(files = paste0("/lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/", samFiles),
                             annot.ext = "/scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.100_introns10.gtf",
                             isGTFAnnotationFile = TRUE,
                             GTF.featureType = "intron10",
                             GTF.attrType = "gene_id",
                             strandSpecific = 0, isPairedEnd = TRUE,
                             nthreads = 28)

intronCountDF <- as.data.frame(intronCount$count)
colnames(intronCountDF) <- gsub("Aligned.out.bam", "", colnames(intronCountDF))
# export intronCount results
saveRDS(intronCount, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawIntronCount_EUR.rds")


# Normalize with DESeq2 -------------------------------------------------
# create design matrix for DESeq2
intronCount <- readRDS("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawIntronCount_EUR.rds")
intronCountDF <- as.data.frame(intronCount$counts)
colnames(intronCountDF) <- gsub("Aligned.out.bam", "", colnames(intronCountDF))

designMat <- data.frame(group = factor(rep("EUR_LCL",
                                           length(colnames(intronCountDF)))))
rownames(designMat) <- colnames(intronCountDF)
all(rownames(designMat) == colnames(intronCountDF))
# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = intronCountDF,
                              colData = designMat, design = ~ 1)
# filter genes whose sum of expression less than one across all samples
ddsF <- dds[ rowSums(counts(dds)) > length(colnames(intronCountDF)), ]

vst=varianceStabilizingTransformation(ddsF)
normIntron <- as.data.frame(assay(vst))

# Rename to match name in vcf file, order samples and Re-add geneID columns
normIntron <- normIntron[, order(names(normIntron))]
normIntron <- cbind(gene_id = rownames(normIntron), normIntron)

cat("This is a view of the normalized expression data \n")
cat("\n")
normIntron[1:10, 1:10]
write.csv(normIntron, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normIntronCount_EUR.csv", row.names = F)

# Prep bed for normalized genes -----------------------------------
# Add coordinate information 
intronCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
normIntronCoord <- intronCoord[which(intronCoord$gene_id %in% normIntron$gene_id), ]
normIntronCoord <- merge(normIntronCoord, normIntron, by = "gene_id")
normIntronCoord <- normIntronCoord[, c(2:4,1,1,5, 7:ncol(normIntronCoord))]
colnames(normIntronCoord)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
normIntronCoord$feature <- paste0(normIntronCoord$feature, "_IntronCount")
normIntronCoord <- normIntronCoord[order(normIntronCoord[,1], normIntronCoord[,2]), ]
normIntronCoord[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
intronFileName <- paste0("phenoFiles/geneIntronCount_EUR_qtlTools_Ready.bed")
fwrite(normIntronCoord, file = intronFileName, sep = "\t")

system(paste0("bgzip ", intronFileName))
system(paste0("tabix ", intronFileName, ".gz"))
rm(list = ls())