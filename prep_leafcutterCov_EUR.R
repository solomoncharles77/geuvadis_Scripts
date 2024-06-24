


genoPCA <- read.table(paste0("covFiles/eurGeno_genoPCs.txt"), check.names=FALSE)
genoPCA <- genoPCA[1:3, ]
phenoPCA <- read.table(paste0("covFiles/geneLeafcutter_EUR_qtlTools_Ready.pca"), header = T)
phenoPCA <- phenoPCA[1:50, ]



cov <- rbind(genoPCA, phenoPCA)

write.table(cov, paste0("covFiles/geneLeafcutter_geno3pc_pheno50pc.txt"), sep = "\t", row.names = F, quote = F)

cat("Covariates file Ready  \n")
