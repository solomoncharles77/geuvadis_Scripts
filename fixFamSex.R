
fam <- read.table("genoFiles/allGeuvadisSampGeno.fam")
si <- read.delim("geuvadis_Scripts/igsr_samples.tsv")

si <- si[si$Sample.name %in% fam$V2, c(1:2) ]
si$V5 <- ifelse(si$Sex == "M", 1, 2)
fam$V5 <- si$V5

write.table(fam, "genoFiles/allGeuvadisSampGeno.fam", row.names = F,
            quote = F, col.names = F)
