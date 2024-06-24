
library(readxl)

pca <- read.table(paste0("covFiles/yriGeno_genoPCs.txt"), check.names=FALSE)

sex <- read_xlsx("geuvadis_Scripts/20130606_sample_info.xlsx")
sex <- sex[sex$Sample %in% colnames(pca), c(1,5)]
sex <- data.frame(t(sex))
colnames(sex) <- sex[1, ]
sex <- sex[-1, ]
sex <- cbind(SampleID = "Sex", sex)

cov <- rbind(sex, pca)

write.table(cov, paste0("covFiles/yriGeno_Sex_10pc.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
