
tg <- read.table("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/target.txt")
tg
tg$group <- as.factor(c("S", "S", "C", "C", "S", "S", "S", "C", "C", "C")) # edit as required.


write.table(tg[, c(1:3)], "/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/target.txt", quote = F, sep="\t")

targetFile <- "/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/target.txt"
target <- read.table(targetFile, header=TRUE, sep="\t", na.strings="")
