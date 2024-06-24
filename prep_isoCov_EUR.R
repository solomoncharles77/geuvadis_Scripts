
library(readxl)

genoPCA <- read.table(paste0("covFiles/eurGeno_genoPCs.txt"), check.names=FALSE)
genoPCA <- genoPCA[1:3, ]
phenoPCA <- read.table(paste0("covFiles/isoformQuant_EUR_qtlTools_Ready.pca"), header = T)
phenoPCA <- phenoPCA[1:50, ]

# sex <- read_xlsx("geuvadis_Scripts/20130606_sample_info.xlsx")
# sex <- sex[sex$Sample %in% colnames(pca), c(1,5)]
# sex <- data.frame(t(sex))
# colnames(sex) <- sex[1, ]
# sex <- sex[-1, ]
# sex <- cbind(Covariates = "Sex", sex)

cov <- rbind(genoPCA, phenoPCA)

write.table(cov, paste0("covFiles/isoformQuant_geno3pc_pheno50pc.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
