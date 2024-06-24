
library(tidyverse)
library(data.table)

annot <- data.frame(fread("Mostafavi_snp_annotations.txt"))


# Function to bin distances
bin_distance <- function(distance) {
  bin <- floor(distance / 10000) * 10
  return(bin)
}


prep4Fractionplot <- function(file_path, feature_name) {
  sumstat <- data.frame(fread(file_path, header = F, stringsAsFactors = F))
  sumstat <- sumstat[sumstat$V14 < 5e-03, c(1, 9)]
  colnames(sumstat) <- c("phenotypeID", "dist")
  
  num_bins <- 10
  
  # Create bins for positive values
  pos <- sumstat %>%
    filter(dist > 0 & dist < 1000000) %>%
    mutate(bins = ntile(dist, num_bins)) %>%
    mutate(bins = paste0(bins, "0"))
  
  neg <- sumstat %>%
    filter(dist < 0 & dist > -1000000) %>%
    mutate(bins = ntile(dist, num_bins)) %>%
    mutate(bins = paste0("-", bins, "0"))
  
  sumstat2 <- rbind(pos, neg)
  sumstat2$bins2 <- sapply(sumstat2$dist, bin_distance)
  
  sumstat_agg <- aggregate(dist ~ bins2, data = sumstat2, FUN = length)
  sumstat_agg$fraction <- sumstat_agg$dist / nrow(sumstat2)
  
  sumstat_agg$bins2 <- as.numeric(sumstat_agg$bins2)
  sumstat_agg$feature <- feature_name
  
  return(sumstat_agg)
}

isoformQuant <- prep4Fractionplot("isoformQuant_EUR_cisEQTL/isoformQuant_EUR_cisEQTL_nominal_ALL.txt.gz", "isoformQuant")
leafcutter <- prep4Fractionplot("geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_nominal_ALL.txt.gz", "leafcutter")
exonCount <- prep4Fractionplot("geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_nominal_ALL.txt.gz", "exonCount")
intronCount <- prep4Fractionplot("geneIntronCount_EUR_cisEQTL/geneIntronCount_EUR_cisEQTL_nominal_ALL.txt.gz", "intronCount")
exonIntronRatio <- prep4Fractionplot("geneExonIntronRatio_EUR_cisEQTL/geneExonIntronRatio_EUR_cisEQTL_nominal_ALL.txt.gz", "exonIntronRatio")
comboFeature <- prep4Fractionplot("combinedPheno_EUR_cisEQTL/combinedPheno_EUR_cisEQTL_nominal_ALL.txt.gz", "comboFeature")


# # Plot single fractions ---------------------------------------------------
exonIntronRatio %>% 
  ggplot(aes(x = bins2, y = fraction)) +
  geom_point() +
  geom_line() +
  labs(x = "Distance to TSS (kb)", y = "Fraction of SNPs") +
  theme_bw() +
  scale_x_continuous(limits = c(-200, 200), breaks = seq(-200, 200, 50))


# Plot combined fractions -------------------------------------------------
rbind(isoformQuant, leafcutter, exonCount, intronCount, exonIntronRatio, comboFeature) %>% 
  group_by(feature) %>%
  ggplot(aes(x = bins2, y = fraction, colour = factor(feature))) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Distance to TSS (kb)", y = "Fraction of SNPs", color = "Feature") +
  theme_bw() +
  scale_x_continuous(limits = c(-200, 200), breaks = seq(-200, 200, 50))

ggsave("fractionPlot.pdf", path = "outputFilesnPlots",
       units = "in", width = 5.5, height = 3,       )

rbind(isoformQuant, leafcutter, exonCount, intronCount, exonIntronRatio, comboFeature) %>% 
  group_by(feature) -> fractionPlotDF

write.csv(fractionPlotDF, "outputFilesnPlots/fractionPlotData.csv", row.names = F)








d <- data.frame(fread("isoformQuant_EUR_cisEQTL/isoformQuant_EUR_cisEQTL_permute_ALL.txt.gz"))
annot <- data.frame(fread("../MostafaviFiles/gwas_eqtl/snp_annotations/filter_snps.txt"))

dAnnot <- merge(d[, c("V10", "V1")], annot, by.x = "V10", by.y = "SNP")
head(dAnnot)
