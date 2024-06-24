
options(scipen = 999, "warnPartialMatchDollar"=TRUE)
library(data.table)
library(rtracklayer)

# conda activate py27

annot <- data.frame(fread("../MostafaviFiles/gwas_eqtl/snp_annotations/filter_snps.txt"))
annotBed <- annot
annotBed$chr <- sapply( annotBed$SNP, function(x){unlist(strsplit(x, split = ":"))[1]})
annotBed$chr <- paste0("chr", annotBed$chr)
annotBed$position <- sapply( annotBed$SNP, function(x){unlist(strsplit(x, split = ":"))[2]})

# tt <- annot
# annot <- tt

annotBed <- annot[, c("chr", "position", "position", "rsID")]
colnames(annotBed) <- c("chr", "start", "end", "rsid")
# 
# Export BED file and liftover coordinates with Crossmap ------------------
outFileName <- "snpAnnot"
bedFile19 <- paste0("../colocalization/summStatsBEDs/", outFileName, "_hg19_BED.txt")
bedFile38 <- paste0("../colocalization/summStatsBEDs/", outFileName, "_hg38_BED.txt")

fwrite(annotBed, file = bedFile19, sep = "\t", col.names = F, row.names = F, quote = F)
#rm(annotBed_hg19)

cat("\n")
cat("Lifting over coordinates to hg38 ....... \n")
system(paste0("CrossMap.py bed ../colocalization/hg19ToHg38.over.chain ", bedFile19, " ", bedFile38))


annotBED_hg38 <- data.frame(fread(bedFile38))
colnames(annotBED_hg38) <- c("chr", "start", "hg38_bp", "rsID")

annot_hg38 <- merge(annot, annotBED_hg38[, -c(1,2)], by = "rsID")

head(annot_hg38)
