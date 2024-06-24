

suppressMessages(library(data.table))



# Import and and harmonize exon and intron counts -------------------------------------------
normExon <- data.frame(fread("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normExonCount_EUR.csv"))
rownames(normExon) <- normExon$gene_id
normExon <- normExon[, -c(1)]
normExon[1:10, 1:10]

normIntron <- data.frame(fread("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normIntronCount_EUR.csv"))
rownames(normIntron) <- normIntron$gene_id
normIntron <- normIntron[, -c(1,2)]
normIntron[1:10, 1:10]

normExon <- normExon[rownames(normExon) %in% rownames(normIntron), ]
normIntron <- normIntron[rownames(normIntron) %in% rownames(normExon), ]

# Calculate log2 
log2NormExon <- log2(normExon)
log2NormIntron <- log2(normIntron)

# Calculate the difference between log2
exonIntronRatio <- log2NormExon - log2NormIntron
exonIntronRatio <- exonIntronRatio[, order(names(exonIntronRatio))]
exonIntronRatio <- cbind(gene_id = rownames(exonIntronRatio), exonIntronRatio)

# Export ratio
write.csv(exonIntronRatio, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/exonIntronRatio_EUR.csv", row.names = F)

# Prep bed for QTLtools -----------------------------------
# Add coordinate information 
coord <- data.frame(fread("/scratch/vasccell/cs806/exprPhenoData/VSMC_Gene_Coords_withStrand.txt"))
exonIntronRatioCoord <- coord[which(coord$gene_id %in% exonIntronRatio$gene_id), ]
exonIntronRatioCoord <- merge(exonIntronRatioCoord, exonIntronRatio, by = "gene_id")
exonIntronRatioCoord <- exonIntronRatioCoord[, c(2:4,1,1,5, 7:ncol(exonIntronRatioCoord))]
colnames(exonIntronRatioCoord)[1:6] <- c("#chr", "start", "end", "feature", "gene", "strand")
exonIntronRatioCoord$feature <- paste0(exonIntronRatioCoord$feature, "_ExonIntronRatio")
exonIntronRatioCoord <- exonIntronRatioCoord[order(exonIntronRatioCoord[,1], exonIntronRatioCoord[,2]), ]
exonIntronRatioCoord[1:10, 1:10]
cat("\n")
cat("Exporting QTL ready files \n")
cat("\n")

# Export files for QTLtools ---------------------------------------------
exonIntronRatioFileName <- paste0("phenoFiles/geneExonIntronRatio_EUR_qtlTools_Ready.bed")
fwrite(exonIntronRatioCoord, file = exonIntronRatioFileName, sep = "\t")

system(paste0("bgzip -f ", exonIntronRatioFileName))
system(paste0("tabix ", exonIntronRatioFileName, ".gz"))
rm(list = ls())



