suppressMessages(library(Rsubread))
suppressMessages(library(DESeq2))


samFiles <- readLines("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSam.txt")

clean <- featureCounts(files = paste0("/lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/", samFiles),
                       annot.ext = "/scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.100.gtf",
                       isGTFAnnotationFile = TRUE,
                       GTF.featureType = "exon",
                       GTF.attrType = "gene_id",
                       strandSpecific = 2, isPairedEnd = TRUE,
                       nthreads = 8)

count <- as.data.frame(clean$counts)
colnames(count) <- gsub("Aligned.out.sam", "", colnames(count))
# export featureCount results
write.csv(count, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/rawReadCount.csv", row.names = F)
saveRDS(clean, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/readFeatures.rds")

# export individual count files
for (i in seq_along(colnames(clean$counts))) {
  df <- cbind(clean$annotation, clean$counts[, i])
  colnames(df)[7] <- "count"
  write.table(df, file = paste0("/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/",
                                sub("Aligned.out.sam", "", colnames(clean$counts)[i]),
                                "_counts.txt"), sep = "\t", quote = FALSE)
}


# Create target file for SARTools
# Create an empty list to store data frames
df_list <- list()

# Loop over columns
for (i in seq_along(colnames(clean$counts))) {
  # Create a data frame for each iteration
  df <- data.frame(
    label = sub("Aligned.out.sam", "", colnames(clean$counts)[i]),
    files = paste0(sub("Aligned.out.sam", "", colnames(clean$counts)[i]), "_counts.txt"),
    group = "toAdd",
    sex = "toAdd",
    info = "toAdd"
  )
  
  # Add the data frame to the list
  df_list[[i]] <- df
}

# Combine all data frames into a single data frame
tg <- do.call(rbind, df_list)
write.table(tg, "/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/target.txt", quote = F, sep="\t")


# # Normalize with DESeq2 -------------------------------------------------
# # create design matrix for DESeq2
# designMat <- data.frame(group = factor(rep("SMC",
#                                            length(colnames(count)))))
# rownames(designMat) <- colnames(count)
# all(rownames(designMat) == colnames(count))
# # create DESeq2 object
# dds <- DESeqDataSetFromMatrix(countData = count,
#                                  colData = designMat, design = ~ 1)
# # filter genes whose sum of expression less than one across all samples
# ddsF <- dds[ rowSums(counts(dds)) > length(colnames(count)), ]
# 
# # get and export DESEq2 normalized count.
# est <- estimateSizeFactors(ddsF)
# norm <- counts(est,normalized=T)
# write.csv(norm, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/readCount/normReadVals.csv")
