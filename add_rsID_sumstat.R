
library(tidyverse)
library(data.table)

qtlHeader <- read.delim("geuvadis_Scripts/QTLtools_nominal_header.txt", header = F)
colIDs <- qtlHeader$V2
colIDs <- sub(" \\| ve_by_pc1 \\| n_phe_in_grp", "", colIDs)
colIDs <- sub("phe_id \\| ", "", colIDs)

addRSID <- function(file_path, feature_name) {
  sumstat <- data.frame(fread(file_path, header = F, stringsAsFactors = F))
  colnames(sumstat) <- colIDs
  
  coordFileOut <- paste0("../colocalization/summStatsBEDs/", feature_name, "_coordID.txt")
  fwrite(sumstat[, c("var_id", "var_id")], coordFileOut, sep = "\t", quote = F, row.names = F, col.names = F)
  coord <- system(paste0("/home/c/cs806/tsv-utils-v2.1.2_linux-x86_64_ldc2/bin/tsv-join -f ", coordFileOut, " -k 1 -d 2 ../colocalization/dbSNP/dbSNP155.hg38.rsID_CoordID.txt"), intern = T)
  coord <- data.frame(fread(text = coord, header = F))
  colnames(coord) <- c("rsID", "var_id")
  
  sumstat2 <- merge(sumstat, coord, by = "var_id", all.x = TRUE)
  
  return(sumstat2)
}


isoformQuant <- addRSID("isoformQuant_EUR_cisEQTL/isoformQuant_EUR_cisEQTL_nominal_ALL.txt.gz", "isoformQuant")
fwrite(isoformQuant, "isoformQuant_EUR_cisEQTL/isoformQuant_EUR_cisEQTL_nominal_ALL_rsID.txt.gz", sep = "\t")
gc()

leafcutter <- addRSID("geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_nominal_ALL.txt.gz", "leafcutter")
fwrite(leafcutter, "geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_nominal_ALL_rsID.txt.gz", sep = "\t")
gc()

exonCount <- addRSID("geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_nominal_ALL.txt.gz", "exonCount")
fwrite(exonCount, "geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_nominal_ALL_rsID.txt.gz", sep = "\t")
gc()

intronCount <- addRSID("geneIntronCount_EUR_cisEQTL/geneIntronCount_EUR_cisEQTL_nominal_ALL.txt.gz", "intronCount")
fwrite(intronCount, "geneIntronCount_EUR_cisEQTL/geneIntronCount_EUR_cisEQTL_nominal_ALL_rsID.txt.gz", sep = "\t")
gc()

exonIntronRatio <- addRSID("geneExonIntronRatio_EUR_cisEQTL/geneExonIntronRatio_EUR_cisEQTL_nominal_ALL.txt.gz", "exonIntronRatio")
fwrite(exonIntronRatio, "geneExonIntronRatio_EUR_cisEQTL/geneExonIntronRatio_EUR_cisEQTL_nominal_ALL_rsID.txt.gz", sep = "\t")
gc()

comboFeature <- addRSID("combinedPheno_EUR_cisEQTL/combinedPheno_EUR_cisEQTL_nominal_ALL.txt.gz", "comboFeature")
fwrite(comboFeature, "combinedPheno_EUR_cisEQTL/combinedPheno_EUR_cisEQTL_nominal_ALL_rsID.txt.gz", sep = "\t")
gc()


