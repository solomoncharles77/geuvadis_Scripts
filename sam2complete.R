
cram <- read.delim("geuvadis_Scripts/eurCram.txt", header = F)
cram$V1 <- gsub("mappedReads/cramFiles/", "", cram$V1)
cram$V1 <- gsub(".sorted.bam.cram", "Aligned.out.bam", cram$V1)
sam <- read.delim("geuvadis_Scripts/allSam.txt", header = F)

sam2 <- data.frame(sam[!sam$V1 %in% cram$V1, ])

write.table(sam2, "geuvadis_Scripts/allSam_leftOver.txt", row.names = F, quote = F, col.names = F)

head(cram)
head(sam2)
