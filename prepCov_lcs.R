
library(readxl)

genoPC <- read.table(paste0("covFiles/yriGeno_genoPCs.txt"), check.names=FALSE)
splicePC <- read.table("phenoFiles/geuYri_LCsplice_qtlTools_Ready.bed_IntronPCs.txt", header = T)

sex <- read_xlsx("geuvadis_Scripts/20130606_sample_info.xlsx")
sex <- sex[sex$Sample %in% colnames(genoPC), c(1,5)]
sex <- data.frame(t(sex))
colnames(sex) <- sex[1, ]
sex <- sex[-1, ]
sex <- cbind(Covariates = "Sex", sex)

cov <- rbind(sex, genoPC, splicePC)

write.table(cov, paste0("covFiles/geuYri_Sex_10GenoPC_10SplicePC.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
