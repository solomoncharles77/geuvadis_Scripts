

sampIndex <- read.delim("geuvadis_Scripts/geuvadis.sequence.index", skip = 23)
yri <- sampIndex[sampIndex$POPULATION == "YRI" & sampIndex$ANALYSIS_GROUP == "mRNA", c(3, 10, 19) ]
yri2 <- data.frame(wget = "wget -nc", yri$PAIRED_FASTQ)

write.table(yri, "geuvadis_YRI_mRNA_detail_extract.txt", row.names = F, quote = F)
write.table(yri2, "rawReads/geuvadis_YRI_mRNA_download.sh", row.names = F, quote = F, col.names = F)


sampIndex <- read.delim("geuvadis_Scripts/geuvadis.sequence.index", skip = 23)
eur <- sampIndex[sampIndex$POPULATION != "YRI" & sampIndex$ANALYSIS_GROUP == "mRNA", c(3, 10, 19)  ]
eur <- eur[-21, ]
eur2 <- data.frame(wget = "wget -nc", eur$PAIRED_FASTQ)

write.table(eur, "geuvadis_EUR_mRNA_detail_extract.txt", row.names = F, quote = F)
write.table(eur2, "rawReads/geuvadis_EUR_mRNA_download.sh", row.names = F, quote = F, col.names = F)



head(yri2)

unique(yri$ANALYSIS_GROUP)
