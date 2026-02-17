library(data.table)

csEdCoord <- read.table("geuvadis_Scripts/codingSequenceEditingSitesCoordinates.txt")

# Import raw editing data
freqDF <- data.frame(fread("phenoFiles/rawGeuvadisRediKnown.txt.gz"))
colnames(freqDF)[1] <- "coordID"
rownames(freqDF) <- freqDF$coordID

cdsDF <- freqDF[freqDF$coordID %in% csEdCoord$hg38_ID, ]
cdsDF$coordID <- NULL

# Calculate mean editing level per sample (column), ignoring NA
cdsEditLevel <- colMeans(cdsDF, na.rm = TRUE)

# Convert to tidy data frame
kelDF <- data.frame(
  sample = names(cdsEditLevel),
  cds_editing = round(cdsEditLevel, 4),
  n_sites = colSums(!is.na(cdsDF)),  # number of non-NA sites per sample
  stringsAsFactors = FALSE)

head(kelDF)

############################################################
rnae <- data.frame(fread("phenoFiles/geuvadisMMSitesAEI_clean_AssocReady.txt"))
rnaeNov <- merge(rnae[, c(1,2,4)], kelDF, by.x = "V2", by.y = "sample")
head(rnaeNov)
cor(rnaeNov$V4, rnaeNov$cds_editing)

kelPheno <- rnaeNov[, c(2,1,4)]
fam <- read.table("../genotypeData/rawGeno/newImpute2025/huvecGeno.fam")
kelPheno <- kelPheno[match(fam$V2, kelPheno$V2), ]
colnames(kelPheno) <- c("FID",   "IID", "cdsEditLevel")
head(kelPheno)
write.table(kelPheno, "phenoFiles/geuvadisCodingSequenceEditingLevel_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")

png("resPlots/geuvadisAverageEditCodingSequences.png", width=1200, height=600)
hist(kelDF$cds_editing, col = "lightblue",
     main = "Average A-to-I Editing in Coding Sequences",
     xlab = "",
     breaks = 30)
dev.off()
