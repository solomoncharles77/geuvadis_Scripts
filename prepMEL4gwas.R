library(data.table)

knownEdit <- data.frame(fread("phenoFiles/filteredGeuvadisRediKnown.txt.gz"))
rownames(knownEdit) <- knownEdit$coordID
knownEdit$coordID <- NULL

# Calculate mean editing level per sample (column), ignoring NA
knownEditLevel <- colMeans(knownEdit, na.rm = TRUE)

# Convert to tidy data frame
kelDF <- data.frame(
  sample = names(knownEditLevel),
  overall_editing = round(knownEditLevel, 4),
  n_sites = colSums(!is.na(knownEdit)),  # number of non-NA sites per sample
  stringsAsFactors = FALSE)

head(kelDF)

############################################################
rnae <- data.frame(fread("phenoFiles/geuvadisRefSeqAEI_clean_AssocReady.txt"))
rnaeNov <- merge(rnae[, c(1,2,4)], kelDF, by.x = "V2", by.y = "sample")
head(rnaeNov)
cor(rnaeNov$V4, rnaeNov$overall_editing)

kelPheno <- rnaeNov[, c(2,1,4)]
fam <- read.table("../genotypeData/rawGeno/newImpute2025/huvecGeno.fam")
kelPheno <- kelPheno[match(fam$V2, kelPheno$V2), ]
colnames(kelPheno) <- c("FID",   "IID", "meanEditLevel")
head(kelPheno)
write.table(kelPheno, "phenoFiles/geuvadisMeanEditingLevel_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")

###############################################################
# prep covariate
genoPCA <- read.table(paste0("covFiles/geuvadisCovs.txt"), check.names=FALSE)
genoPCA$Covariates <- NULL
genoPCA <- t(genoPCA[1:3, ])
cov <- data.frame(fam[, c(1,2,5)], genoPCA)
colnames(cov)[1:3] <- c("FID",   "IID", "Sex")

write.table(cov, "covFiles/geuvadis_Sex_3gpc_covariates.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(cov, "covFiles/geuvadis_Sex_3gpc_covariates_wtHeader.txt", row.names = F, quote = F, sep = "\t")
########################################################################################################
########################################################################################################


png("resPlots/meanGlobalEditing.png", width=1200, height=600)
hist(kelDF$overall_editing, col = "lightyellow",
     main = "Mean A-to-I Editing in Known Sites",
     xlab = "",
     breaks = 30)
dev.off()
