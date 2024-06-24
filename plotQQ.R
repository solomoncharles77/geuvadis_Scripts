
library(tidyverse)
library(data.table)

prep4QQlot <- function(file_path, feature_name) {
  sumstat <- data.frame(fread(file_path, header = F, stringsAsFactors = F))
  sumstat <- sumstat$V14
  return(sumstat)
}

isoformQuant <- prep4QQlot("isoformQuant_EUR_cisEQTL/isoformQuant_EUR_cisEQTL_nominal_ALL.txt.gz", "isoformQuant")
leafcutter <- prep4QQlot("geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_nominal_ALL.txt.gz", "leafcutter")
exonCount <- prep4QQlot("geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_nominal_ALL.txt.gz", "exonCount")
intronCount <- prep4QQlot("geneIntronCount_EUR_cisEQTL/geneIntronCount_EUR_cisEQTL_nominal_ALL.txt.gz", "intronCount")
exonIntronRatio <- prep4QQlot("geneExonIntronRatio_EUR_cisEQTL/geneExonIntronRatio_EUR_cisEQTL_nominal_ALL.txt.gz", "exonIntronRatio")
comboFeature <- prep4QQlot("combinedPheno_EUR_cisEQTL/combinedPheno_EUR_cisEQTL_nominal_ALL.txt.gz", "comboFeature")


df <- data.frame(
  p_value = c(isoformQuant, leafcutter, exonCount, intronCount, exonIntronRatio, comboFeature),
  group = rep(c("isoformQuant", "leafcutter", "exonCount", "intronCount", "exonIntronRatio", "comboFeature"), 
              times = c(length(isoformQuant), length(leafcutter), length(exonCount), length(intronCount), length(exonIntronRatio), length(comboFeature)))
)


ggplot(df, aes(sample = -log10(p_value), color = group)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "Q-Q Plot of P-Values",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles",
       color = "Group") 
  scale_color_manual(values = c("blue", "red", "green"))


# Create a Q-Q plot using qqplotr
ggplot(data.frame(sample = data), aes(sample = sample)) +
  stat_qq(size = 0.5) +
  stat_qq_line(col = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles")
