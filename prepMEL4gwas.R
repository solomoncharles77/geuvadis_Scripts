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

fam <- read.table("genoFiles/allGeuvadisSampGeno.fam")
kelPheno <- merge(fam[, c(1,2)], kelDF, by.x = "V2", by.y = "sample" )
kelPheno <- kelPheno[match(fam$V2, kelPheno$V2), c(2,1,3)]
colnames(kelPheno) <- c("FID",   "IID", "meanEditLevel")
head(kelPheno)
write.table(kelPheno, "phenoFiles/geuvadisMEL_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/geuvadisMEL_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

# ###############################################################
# # prep covariate
# genoPCA <- read.table(paste0("covFiles/allGeuvadisSamp_genoPCs.txt"), check.names=FALSE)
# genoPCA$SampleID <- NULL
# genoPCA <- t(genoPCA[1:20, ])
# cov <- data.frame(fam[, c(1,2,5)], genoPCA)
# colnames(cov)[1:3] <- c("FID",   "IID", "Sex")
# 
# write.table(cov, "covFiles/geuvadis_Sex_3gpc_covariates.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(cov, "covFiles/geuvadis_Sex_3gpc_covariates_wtHeader.txt", row.names = F, quote = F, sep = "\t")
# ########################################################################################################
# ########################################################################################################

# normalize the data using rank-based inverse normal transformation-------
intTrans <- function(x){
  y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(y)
}

kelPheno$meanEditLevel <-  intTrans(kelPheno$meanEditLevel)
kelPhenoNT <- shapiro.test(kelPheno$meanEditLevel)
ifelse(kelPhenoNT$p.value > 0.05, "normal", "not normal")


write.table(kelPheno, "phenoFiles/geuvadisMEL_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/geuvadisMEL_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")
