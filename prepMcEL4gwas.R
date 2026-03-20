library(data.table)

csEdCoord <- read.table("geuvadis_Scripts/codingSequenceEditingSitesCoordinates.txt")

# Import raw editing data
freqDF <- data.frame(fread("phenoFiles/filteredGeuvadisRediKnown.txt.gz"))
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

fam <- read.table("genoFiles/allGeuvadisSampGeno.fam")
kelPheno <- merge(fam[, c(1,2)], kelDF, by.x = "V2", by.y = "sample" )
kelPheno <- kelPheno[, c(2,1,3)]
colnames(kelPheno) <- c("FID",   "IID", "cdsEditLevel")
head(kelPheno)
write.table(kelPheno, "phenoFiles/geuvadisMcEL_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/geuvadisMcEL_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")


# normalize the data using rank-based inverse normal transformation-------
intTrans <- function(x){
  y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(y)
}

kelPheno$cdsEditLevel <-  intTrans(kelPheno$cdsEditLevel)
kelPhenoNT <- shapiro.test(kelPheno$cdsEditLevel)
ifelse(kelPhenoNT$p.value > 0.05, "normal", "not normal")

write.table(kelPheno, "phenoFiles/geuvadisMcEL_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/geuvadisMcEL_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")
