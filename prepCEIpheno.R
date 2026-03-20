

# normalize the data using rank-based inverse normal transformation-------
intTrans <- function(x){
  y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(y)
}

rnae <- read.csv("RNAEditingIndexer/CEI/geuvadisRNASeqSUM/EditingIndex.csv")
rnae$Sample <- sub(".sorted", "", rnae$Sample)

# Harmonize names and samples in geno and expr ----------------------------
mapFile <- read.table("exprFiles/matched_allSampsPheno2.txt")
rnae$Sample <- mapFile$V2[match(rnae$Sample, mapFile$V1)]
rnae <- rnae[complete.cases(rnae), ]
rnae[1:10, c(1,3,5:10)]

rnaeRand <- rnae[rnae$StrandDecidingMethod == "Randomly", ] 
rnaeRef <- rnae[rnae$StrandDecidingMethod == "RefSeqThenMMSites", ] 
rnaeMMS <- rnae[rnae$StrandDecidingMethod == "MMSitesThenRefSeq", ]
write.table(data.frame(rnaeMMS$Sample),
            "genoFiles/geuvadisRNAseqSamps.txt",
            row.names = F, col.names = F, quote = F)

##########################################################################
fam <- read.table("genoFiles/allGeuvadisSampGeno.fam")

rnaeRand <- rnaeRand[, -c(1,2,4)]
# Reorder rnaeRand to match the order of pheno
rnaeRand <- rnaeRand[match(fam$V2, rnaeRand$Sample), ]
rnaeRand <- rnaeRand[complete.cases(rnaeRand), ]
rnaeRand[1:10, 1:10]
fam <- fam[match(rnaeRand$Sample, fam$V2), ]
fam <- fam[complete.cases(fam), ]
head(fam)
rnaeRand <- data.frame(fam[, c(1,2)], rnaeRand[, c(2:7)])
colnames(rnaeRand)[1:2] <- c("FID",   "IID")
head(rnaeRand)
write.table(rnaeRand, "phenoFiles/geuvadisRandomCEI_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(rnaeRand, "phenoFiles/geuvadisRandomCEI_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

rnaeRand[, c(3:8)] <- lapply(rnaeRand[, c(3:8)], function(x) intTrans(x))
normRand <- rnaeRand
normRandNT <- lapply(normRand[, 3:8], function(x) shapiro.test(x)) 
normRandNT <- do.call(rbind.data.frame, normRandNT)
normRandNT$status <- ifelse(normRandNT$p.value > 0.05, "normal", "not normal")
normRandNT

write.table(normRand, "phenoFiles/geuvadisRandomCEI_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(normRand, "phenoFiles/geuvadisRandomCEI_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

##########################################################################################################
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


rnaeRef <- rnaeRef[, -c(1,2,4)]
# Reorder rnaeRef to match the order of pheno
rnaeRef <- rnaeRef[match(fam$V2, rnaeRef$Sample), ]
rnaeRef <- rnaeRef[complete.cases(rnaeRef), ]
rnaeRef[1:10, 1:10]
fam <- fam[match(rnaeRef$Sample, fam$V2), ]
fam <- fam[complete.cases(fam), ]
head(fam)
rnaeRef <- data.frame(fam[, c(1,2)], rnaeRef[, c(2:7)])
colnames(rnaeRef)[1:2] <- c("FID",   "IID")
head(rnaeRef)
write.table(rnaeRef, "phenoFiles/geuvadisRefSeqCEI_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(rnaeRef, "phenoFiles/geuvadisRefSeqCEI_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

rnaeRef[, c(3:8)] <- lapply(rnaeRef[, c(3:8)], function(x) intTrans(x))
normRand <- rnaeRef
normRefNT <- lapply(normRand[, 3:8], function(x) shapiro.test(x)) 
normRefNT <- do.call(rbind.data.frame, normRefNT)
normRefNT$status <- ifelse(normRefNT$p.value > 0.05, "normal", "not normal")
normRefNT

write.table(normRand, "phenoFiles/geuvadisRefSeqCEI_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(normRand, "phenoFiles/geuvadisRefSeqCEI_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")


###################################################
rnaeMMS <- rnaeMMS[, -c(1,2,4)]
# Reorder rnaeMMS to match the order of pheno
rnaeMMS <- rnaeMMS[match(fam$V2, rnaeMMS$Sample), ]
rnaeMMS <- rnaeMMS[complete.cases(rnaeMMS), ]
rnaeMMS[1:10, 1:10]
fam <- fam[match(rnaeMMS$Sample, fam$V2), ]
fam <- fam[complete.cases(fam), ]
head(fam)
rnaeMMS <- data.frame(fam[, c(1,2)], rnaeMMS[, c(2:7)])
colnames(rnaeMMS)[1:2] <- c("FID",   "IID")
head(rnaeMMS)
write.table(rnaeMMS, "phenoFiles/geuvadisMMSitesCEI_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(rnaeMMS, "phenoFiles/geuvadisMMSitesCEI_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

rnaeMMS[, c(3:8)] <- lapply(rnaeMMS[, c(3:8)], function(x) intTrans(x))
normMMS <- rnaeMMS
normMMSNT <- lapply(normMMS[, 3:8], function(x) shapiro.test(x)) 
normMMSNT <- do.call(rbind.data.frame, normMMSNT)
normMMSNT$status <- ifelse(normMMSNT$p.value > 0.05, "normal", "not normal")
normMMSNT

write.table(normMMS, "phenoFiles/geuvadisMMSitesCEI_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(normMMS, "phenoFiles/geuvadisMMSitesCEI_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")
