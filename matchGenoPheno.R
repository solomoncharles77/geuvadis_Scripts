

mrnaSamps <- read.table("geuvadis_EUR_mRNA_detail_extract.txt", header = T)
genoSamps <- read.table("../colocalization/1000Genome/1kGMerge.fam")

mrnaSamps <- mrnaSamps[!duplicated(mrnaSamps$SAMPLE_NAME), c(1,2)]
mrnaSamps$SAMPLE_NAME <- sub("GEUV:", "", mrnaSamps$SAMPLE_NAME)

genoSamps2 <- genoSamps[genoSamps$V2 %in% mrnaSamps$SAMPLE_NAME, ]
mrnaSamps2 <- mrnaSamps[mrnaSamps$SAMPLE_NAME %in% genoSamps2$V2, ]

# write.table(genoSamps2, "genoFiles/matched_YRI_Geno_Pheno.txt", row.names = F, quote = F, col.names = F)
# write.table(mrnaSamps2, "exprFiles/matched_YRI_Pheno.txt", row.names = F, quote = F, col.names = F)

write.table(genoSamps2, "genoFiles/matched_EUR_Geno_Pheno.txt", row.names = F, quote = F, col.names = F)
write.table(mrnaSamps2, "exprFiles/matched_EUR_Pheno.txt", row.names = F, quote = F, col.names = F)


head(genoSamps2)
head(mrnaSamps2)
