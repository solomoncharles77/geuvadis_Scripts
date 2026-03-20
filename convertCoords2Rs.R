library(data.table)

bim <- data.frame(fread("genoFiles/allGeuvadisSampGeno_RS.bim", header = F))
bim$id <- 1:nrow(bim)

rs <- data.frame(fread("genoFiles/geuvadis_rsLookup.txt", header = F))
rs <- rs[complete.cases(rs), ]
rs <- rs[!duplicated(rs$V2), ]

rsCoord <- merge(bim, rs, by.x = "V2", by.y = "V2", all.x = TRUE)
rsCoord <- rsCoord[order(rsCoord$id), ]
rsCoord$V1.y <- ifelse(is.na(rsCoord$V1.y), rsCoord$V2, rsCoord$V1.y)

rsCoord <- rsCoord[, c("V1.x", "V1.y", "V3", "V4", "V5", "V6")]
fwrite(rsCoord, "genoFiles/allGeuvadisSampGeno_RS.bim", sep = "\t", col.names = FALSE)



rsCoord[299455:299460, ]



