suppressMessages(library(Rsubread))
suppressMessages(library(DESeq2))

# # Exon count -----------------------------------------
# samFiles <- readLines("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSamYRI.txt")
# exprSamps <- read.table("/lustre/alice3/scratch/vasccell/cs806/geuvadis/exprFiles/matched_YRI_Pheno.txt")
# exprSamps$V1 <- paste0(exprSamps$V1, "Aligned.out.sam")
# samFiles <- samFiles[samFiles %in% exprSamps$V1]
# samFiles <- sub(".sam", ".sam.sorted.bam", samFiles)
# 
# 
# exonCount <- featureCounts(files = paste0("/lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/geuvadisYRI_MappedReads/", samFiles),
#                        annot.ext = "/scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.100_exons.gtf",
#                        isGTFAnnotationFile = TRUE,
#                        GTF.featureType = "exon",
#                        GTF.attrType = "gene_id",
#                        strandSpecific = 0, isPairedEnd = TRUE,
#                        nthreads = 28)
# 
# exonCountDF <- as.data.frame(exonCount$counts)
# colnames(exonCountDF) <- gsub("Aligned.out.sam.sorted.bam", "", colnames(exonCountDF))
# # export exonCount results
# write.csv(exonCountDF, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawExonCount_YRI.csv", row.names = F)
# saveRDS(exonCount, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawExonCount_YRI.rds")
# 
# # Normalize with DESeq2 -------------------------------------------------
# # create design matrix for DESeq2
# exonCount <- readRDS("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawExonCount_YRI.rds")
# exonCountDF <- as.data.frame(exonCount$counts)
# colnames(exonCountDF) <- gsub("Aligned.out.sam.sorted.bam", "", colnames(exonCountDF))
# 
# designMat <- data.frame(group = factor(rep("YRI_LCL",
#                                            length(colnames(exonCountDF)))))
# rownames(designMat) <- colnames(exonCountDF)
# all(rownames(designMat) == colnames(exonCountDF))
# # create DESeq2 object
# dds <- DESeqDataSetFromMatrix(countData = exonCountDF,
#                               colData = designMat, design = ~ 1)
# # filter genes whose sum of expression less than one across all samples
# ddsF <- dds[ rowSums(counts(dds)) > length(colnames(exonCountDF)), ]
# 
# vst=varianceStabilizingTransformation(ddsF)
# normExon <- as.data.frame(assay(vst))
# 
# # Rename to match name in vcf file, order samples and Re-add geneID columns
# names(normExon) <- exprSamps$V2[match(names(normExon), exprSamps$V1)]
# normExon <- normExon[, order(names(normExon))]
# normExon <- cbind(gene_id = rownames(normExon), normExon)
# 
# cat("This is a view of the normalized expression data \n")
# cat("\n")
# normExon[1:10, 1:10]
# write.csv(norm, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normExonCount_YRI.csv")
# 
# # Prep bed for normalized genes -----------------------------------
# # Add coordinate information 
# exonCoord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
# normExonCoord <- exonCoord[which(exonCoord$gene_id %in% normExon$gene_id), ]
# normExonCoord <- merge(normExonCoord, normExon, by = "gene_id")
# normExonCoord <- normExonCoord[, c(2:4,1,1,5, 7:ncol(normExonCoord))]
# colnames(normExonCoord)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
# normExonCoord$feature <- paste0(normExonCoord$feature, "_ExonCount")
# normExonCoord <- normExonCoord[order(normExonCoord[,1], normExonCoord[,2]), ]
# normExonCoord[1:10, 1:10]
# cat("\n")
# cat("Exporting QTL ready files \n")
# cat("\n")
# 
# # Export files for QTLtools ---------------------------------------------
# exonFileName <- paste0("phenoFiles/geneExonCount_YRI_qtlTools_Ready.bed")
# fwrite(normExonCoord, file = exonFileName, sep = "\t")
# 
# system(paste0("bgzip ", exonFileName))
# system(paste0("tabix ", exonFileName, ".gz"))
# 
# rm(list = ls())


##################################################################################################################
# Intron count -----------------------------------------

samFiles <- readLines("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSamYRI.txt")
exprSamps <- read.table("/lustre/alice3/scratch/vasccell/cs806/geuvadis/exprFiles/matched_YRI_Pheno.txt")
exprSamps$V1 <- paste0(exprSamps$V1, "Aligned.out.sam")
samFiles <- samFiles[samFiles %in% exprSamps$V1]
samFiles <- sub(".sam", ".sam.sorted.bam", samFiles)

intronCount <- featureCounts(files = paste0("/lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/geuvadisYRI_MappedReads/", samFiles),
                           annot.ext = "/scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.100_introns10.gtf",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "intron10",
                           GTF.attrType = "gene_id",
                           strandSpecific = 0, isPairedEnd = TRUE,
                           nthreads = 28)

intronCountDF <- as.data.frame(intronCount$counts)
colnames(intronCountDF) <- gsub("Aligned.out.sam.sorted.bam", "", colnames(exonCountDF))
# export intronCount results
write.csv(intronCountDF, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawIntronCount_YRI.csv", row.names = F)
saveRDS(intronCount, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawIntronCount_YRI.rds")


# Normalize with DESeq2 -------------------------------------------------
# create design matrix for DESeq2
intronCount <- readRDS("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawIntronCount_YRI.rds")
intronCountDF <- as.data.frame(intronCount$counts)
colnames(intronCountDF) <- gsub("Aligned.out.sam.sorted.bam", "", colnames(intronCountDF))

designMat <- data.frame(group = factor(rep("YRI_LCL",
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
names(normIntron) <- exprSamps$V2[match(names(normIntron), exprSamps$V1)]
normIntron <- normIntron[, order(names(normIntron))]
normIntron <- cbind(gene_id = rownames(normIntron), normIntron)

cat("This is a view of the normalized expression data \n")
cat("\n")
normIntron[1:10, 1:10]
write.csv(normIntron, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normIntronCount_YRI.csv")

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
intronFileName <- paste0("phenoFiles/geneIntronCount_YRI_qtlTools_Ready.bed")
fwrite(normIntronCoord, file = intronFileName, sep = "\t")

system(paste0("bgzip ", intronFileName))
system(paste0("tabix ", intronFileName, ".gz"))
