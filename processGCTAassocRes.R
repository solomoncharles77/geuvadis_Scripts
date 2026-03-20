
library(data.table)
library(qqman)
library(tidyverse)

tagDir <- "assocRes2/"
tagFiles <- list.files(tagDir, "*.mlma")
tagFiles <- tagFiles[!grepl("adjusted", tagFiles)]
tagFiles <- tagFiles[!grepl(".txt.gz", tagFiles)]
tagFiles <- tagFiles[!grepl(".clumped", tagFiles)]
tagFiles <- tagFiles[!grepl(".log", tagFiles)]
tagFiles <- tagFiles[!grepl(".txt", tagFiles)]
tagFiles <- tagFiles[!grepl(".csv", tagFiles)]

tagID <- c("Alu Editing Index",
           "Cytoplasmic Editing Index",
           "Mean Coding Sequence Editing Level",
           "Mean Editing Level")

# tagID <- c("Mean Coding Sequence Editing Level",
#            "Alu Editing Index",
#            "Cytoplasmic Editing Index",
#            "Mean Editing Level")
tagDF <- data.frame(tF = tagFiles, tD = tagID )
tagDF
#tagDF <- tagDF[2,]


lapply(1:nrow(tagDF), function(x){
  df = tagDF[x, 1]
  id = tagDF[x, 2]
  gwas <- data.frame(fread(paste0(tagDir, df)))
  # Rename columns
  gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                   c("CHR", "BP", "REF", "BETA", "SE", "P"), skip_absent = T)
  # Remove rows with missing p-values and clean up rsID
  gwas <- gwas[!is.na(gwas$P), ]
  gwas$SNP <- sub("_.*", "", gwas$SNP)
  
  # export manhattan plot
  png(paste0(tagDir, "manPlots/", df, "_manhattanPlot.png"), width=1200, height=600)
  manhattan(gwas, chr="CHR", bp="BP", snp="SNP", p="P",
            main=paste0(id))
  dev.off()
  
  png(paste0(tagDir, "qqPlots/", df, "_qqplot.png"), width=400, height=300) 
  qq <- qq(gwas$P)
  abline(h = -log10(5e-8), col = "red", lty = 2)
  title(paste0(id))
  dev.off()  
  
  return(NULL)
  
})

########################################################################################
# Extract hits for clumping

lapply(tagFiles, function(x){
  gwas <- data.frame(fread(paste0(tagDir, x)))
  
  # Rename columns
  gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                   c("CHR", "BP", "REF", "BETA", "SE", "P"), skip_absent = T)
  
  coordFileOut <- paste0("../colocalization/summStatsBEDs/gwasRes_coordID.txt")
  fwrite(gwas[, c("SNP", "SNP")], coordFileOut, sep = "\t", quote = F, row.names = F, col.names = F)
  coord <- system(paste0("/home/c/cs806/tsv-utils-v2.1.2_linux-x86_64_ldc2/bin/tsv-join -f ", coordFileOut, " -k 1 -d 2 ../colocalization/dbSNP/hg38.snp151_All.bed.rsID_CoordID.txt"), intern = T)
  coord <- data.frame(fread(text = coord, header = F))
  colnames(coord) <- c("rsID", "SNP")
  
  gwas <- merge(gwas, coord, by = "SNP", all.x = TRUE)
  gc()
  
  # Remove rows with missing p-values
  gwas <- gwas[!is.na(gwas$P), ]
  gwas <- gwas[order(gwas$P), ]
  hits05 <- gwas[gwas$P <= 5e-05, ]
  hits08 <- gwas[gwas$P <= 5e-08, ]
  
  write.table(hits05, paste0("assocRes/", x, "_e-05_Hits.txt"), row.names = F, quote = F, )
  write.table(hits08, paste0("assocRes/", x, "_e-08_Hits.txt"), row.names = F, quote = F, )
  fwrite(gwas, paste0(tagDir, x, ".txt.gz"), row.names = F, sep = "\t", quote = F)
})  

system2("/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/plinkClump.sh")

# Map to nearest gene
# Map to genes with biomartR
library("biomaRt")
mart <- useMart("ENSEMBL_MART_SNP")
dataset <- useDataset("hsapiens_snp", mart=mart)
library(otargen)

lapply(tagFiles, function(x){
  gwas <- data.frame(fread(paste0(tagDir, x, ".txt.gz")))
  clumpHits <- read.table(paste0("assocRes/", x, "_e-05_Hits.txt"), header = T)
  gwasTop <- clumpHits
  gwasTop$otID <- paste0(gwasTop$CHR, "_", gwasTop$BP, "_",
                         gwasTop$REF, "_", gwasTop$A1)
  
  # Define the list of rsIDs
  rsIDs = gwasTop$rsID
  otIDs <- gwasTop$otID
  
  # To get the ensembl gene id belonging to the SNPs
  annotHit <- getBM(attributes=c('refsnp_id', 'ensembl_gene_stable_id',
                                 "ensembl_transcript_stable_id"), 
                    filters = 'snp_filter', 
                    values = rsIDs, 
                    mart = dataset)
  
  annotGWAS <- merge(gwasTop, annotHit, by.x = "rsID", by.y = "refsnp_id")
  
  head(annotGWAS)
  write.csv(annotGWAS,paste0("assocRes/", x, "_e-05_Hits_nearestGeneEnsembl.csv"), row.names = F )
  
  
  otAnnot <- lapply(otIDs, function(x){
    cat(x, "\n")
    result <- variantEffectPredictorQuery(variantId = x)
    if (!is.null(result)) {
      result <- result[, c(6, 9:10, 14:17)]
      result <- result[result$target.biotype == "protein_coding", ]
      result <- result[order(abs(result$distanceFromTss)), ]
      result <- result[1, ]
      
    }else{
      result
    }
    
  })
  
  otAnnot <- do.call(rbind, otAnnot)
  
  annotGWAS2 <- if(nrow(otAnnot != 0)){ 
    merge(annotGWAS, otAnnot, by.x = "otID", by.y = "variantId", all.x = T)
  }else{ 
    annotGWAS}
  
  annotGWAS2$comGeneID <- ifelse(is.na(annotGWAS2$ensembl_gene_stable_id),
                                 annotGWAS2$target.id, annotGWAS2$ensembl_gene_stable_id)
  annotGWAS2$comTxID <- ifelse(is.na(annotGWAS2$ensembl_transcript_stable_id),
                               annotGWAS2$transcriptId, annotGWAS2$ensembl_transcript_stable_id)
  
  write.csv(annotGWAS2, paste0("resFiles/", x, "_e-05_Hits_mappedGenes.csv"), row.names = F)
  write.table(data.frame(unique(annotGWAS2$comGeneID)), paste0("resFiles/", x, "_e-05_Hits_nearestGenesList.txt"),
              row.names = F, col.names = F, quote = F)
  write.table(data.frame(unique(annotGWAS2$comTxID)), paste0("resFiles/", x, "_e-05_Hits_nearestTranscriptList.txt"),
              row.names = F, col.names = F, quote = F)
})


lapply(tagFiles, function(x){
  gwas <- data.frame(fread(paste0(tagDir, x)))
  clumpHits <- read.table(paste0("assocRes/", x, "_e-08_Hits.txt"), header = T)
  gwasTop <-clumpHits
  gwasTop$rsID <- gwasTop$SNP <- sub("_.*", "", gwasTop$SNP)
  gwasTop$otID <- paste0(gwasTop$CHR, "_", gwasTop$BP, "_",
                         gwasTop$REF, "_", gwasTop$A1)
  
  # Define the list of rsIDs
  rsIDs = gwasTop$rsID
  otIDs <- gwasTop$otID
  
  # To get the ensembl gene id belonging to the SNPs
  annotHit <- getBM(attributes=c('refsnp_id', 'ensembl_gene_stable_id',
                                 "ensembl_transcript_stable_id"), 
                    filters = 'snp_filter', 
                    values = rsIDs, 
                    mart = dataset)
  
  annotGWAS <- merge(gwasTop, annotHit, by.x = "rsID", by.y = "refsnp_id")
  
  head(annotGWAS)
  write.csv(gwasTop,paste0("assocRes/", x, "_e-08_Hits_nearestGeneEnsembl.csv"), row.names = F)
  
  
  otAnnot <- lapply(otIDs, function(x){
    cat(x, "\n")
    result <- variantEffectPredictorQuery(variantId = x)
    if (!is.null(result)) {
      result <- result[, c(6, 9:10, 14:17)]
      result <- result[result$target.biotype == "protein_coding", ]
      result <- result[order(abs(result$distanceFromTss)), ]
      result <- result[1, ]
      
    }else{
      result
    }
    
  })
  
  otAnnot <- do.call(rbind, otAnnot)
  
  annotGWAS2 <- if(nrow(otAnnot != 0)){ 
    merge(annotGWAS, otAnnot, by.x = "otID", by.y = "variantId", all.x = T)
  }else{ 
    annotGWAS}
  
  annotGWAS2$comGeneID <- ifelse(is.na(annotGWAS2$ensembl_gene_stable_id),
                                 annotGWAS2$target.id, annotGWAS2$ensembl_gene_stable_id)
  annotGWAS2$comTxID <- ifelse(is.na(annotGWAS2$ensembl_transcript_stable_id),
                               annotGWAS2$transcriptId, annotGWAS2$ensembl_transcript_stable_id)
  
  write.csv(annotGWAS2, paste0("resFiles/", x, "_e-08_Hits_mappedGenes.csv"), row.names = F)
  write.table(data.frame(unique(annotGWAS2$comGeneID)), paste0("resFiles/", x, "_e-08_Hits_nearestGenesList.txt"),
              row.names = F, col.names = F, quote = F)
  write.table(data.frame(unique(annotGWAS2$comTxID)), paste0("resFiles/", x, "_e-08_Hits_nearestTranscriptList.txt"),
              row.names = F, col.names = F, quote = F)
})




######################################################################
######################################################################
######################################################################
# # Apply Benjamini-Hochberg FDR correction
# library(qvalue)
# gwas$P_adj_BH <- p.adjust(gwas$P, method = "BH")
# 
# # calculate q-values using the 'qvalue' package (often preferred)
# q_obj <- qvalue(p = gwas$P)
# gwas$qvalue <- q_obj$qvalues
# 
# # Find significant SNPs at a 5% FDR level
# gwas_snps_bh <- gwas[gwas$P_adj_BH < 0.05, ]
# gwas_snps_qval <- gwas[gwas$qvalue < 0.05, ]


# library(otargen)
# result <- variantEffectPredictorQuery(variantId = "1_154624272_T_C")
# result <- result[, c(6, 9:10, 14:17)]
# result <- result[result$target.biotype == "protein_coding", ]
# result <- result[order(abs(result$distanceFromTss)), ]
# result <- result[1, ]
############################################################################


