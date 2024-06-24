library(tximport)


tx2gene <- read.csv("/scratch/vasccell/cs806/exprPhenoData/tx2gene.ensembl.v100.csv")

samples <- read.table("/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSeq3.txt", header = F)
files <- file.path("/lustre/alice3/scratch/vasccell/cs806/geuvadis/salmonQuant", samples$V1, "quant.sf")
names(files) <- paste0(samples$V1)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
saveRDS(txi.salmon, file = "/lustre/alice3/scratch/vasccell/cs806/geuvadis/salmonQuant/salmonTxi.rds")
