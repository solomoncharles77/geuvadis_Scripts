
library(data.table)
library(qqman)
library(tidyverse)

# Load data
gwas <- data.frame(fread("assocRes/meanEditLevelGWAS.mlma"))

# Rename columns
gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                 c("CHR", "BP", "REF", "BETA", "SE", "P"), skip_absent = T)
# Remove rows with missing p-values
gwas <- gwas[!is.na(gwas$P), ]
gwas$SNP <- sub("_.*", "", gwas$SNP)
gwas <- gwas[order(gwas$P), ]

# Manhattan plot
manhattan(gwas, chr="CHR", bp="BP", snp="SNP", p="P", 
          genomewideline = -log10(5e-8), 
          suggestiveline = -log10(5e-5),
          main="Manhattan Plot of GWAS Results")

png("assocRes/manPlots/meanEditLevel_manhattan_plot.png", width=1200, height=600)
manhattan(gwas, main="",
          suggestiveline = -log10(5e-05))
dev.off()

png("assocRes/qqPlots/meanEditLevel_qqplot.png", width=400, height=300) 
qq <- qq(gwas$P)
abline(h = -log10(5e-8), col = "red", lty = 2)
dev.off()  

fwrite(gwas, "resFiles/mean_geuvadis_A2I_editing_level_GWAS_Summary_Statistics.txt.gz")


################################################################################
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
# 
########################################################################################
# Load data
gwas <- data.frame(fread("assocRes/meanEditLevelGWAS.mlma"))

# Rename columns
gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                 c("CHR", "BP", "REF", "BETA", "SE", "P"), skip_absent = T)
# Remove rows with missing p-values
gwas <- gwas[!is.na(gwas$P), ]
gwas <- gwas[order(gwas$P), ]
suggestive_hits <- gwas[gwas$P <= 0.00005, ]
head(suggestive_hits)

write.table(suggestive_hits, "assocRes/meanEditLevelGWAS_suggestiveHits.txt", row.names = F, quote = F, )
clumpHits <- read.table("assocRes/meanEditLevelGWAS_clumped_suggestive_loci.clumped", header = T)
gwasTop <- merge(clumpHits[, c(3,6)], gwas, by = "SNP")
gwasTop$rsID <- gwasTop$SNP <- sub("_.*", "", gwasTop$SNP)
gwasTop$otID <- paste0(gwasTop$CHR, "_", gwasTop$BP, "_",
                       gwasTop$REF, "_", gwasTop$A1)

# Define the list of rsIDs
rsIDs = gwasTop$rsID
otIDs <- gwasTop$otID

# Map to genes with biomartR
library("biomaRt")
mart <- useMart("ENSEMBL_MART_SNP")
dataset <- useDataset("hsapiens_snp", mart=mart)

# # To show which marts are available
# listMarts()
# # You need the SNP mart
# mart <- useMart("ENSEMBL_MART_SNP")
# # Find homo sapiens
# listDatasets(mart)
# # This will be the dataset we want to use
# dataset <- useDataset("hsapiens_snp", mart=mart)
# # Show available filters
# listFilters(dataset)
# # Now list all available attributes
# listAttributes(dataset)

# To get the ensembl gene id belonging to the SNPs
annotHit <- getBM(attributes=c('refsnp_id', 'ensembl_gene_stable_id',
                               "ensembl_transcript_stable_id"), 
                  filters = 'snp_filter', 
                  values = rsIDs, 
                  mart = dataset)

annotGWAS <- merge(gwasTop, annotHit, by.x = "rsID", by.y = "refsnp_id")

head(annotGWAS)

write.csv(gwasTop, "assocRes/meanEditLevelGWAS_gwasTop.csv", row.names = F )


library(otargen)
result <- variantEffectPredictorQuery(variantId = "1_154624272_T_C")
result <- result[, c(6, 9:10, 14:17)]
result <- result[result$target.biotype == "protein_coding", ]
result <- result[order(abs(result$distanceFromTss)), ]
result <- result[1, ]


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
annotGWAS2 <- merge(annotGWAS, otAnnot, by.x = "otID", by.y = "variantId", all.x = T)

annotGWAS2$comGeneID <- ifelse(is.na(annotGWAS2$ensembl_gene_stable_id),
                               annotGWAS2$target.id, annotGWAS2$ensembl_gene_stable_id)
annotGWAS2$comTxID <- ifelse(is.na(annotGWAS2$ensembl_transcript_stable_id),
                             annotGWAS2$transcriptId, annotGWAS2$ensembl_transcript_stable_id)

write.csv(annotGWAS2, "resFiles/meanEditLevelGWAS_mappedGenes.csv", row.names = F)
write.table(data.frame(unique(annotGWAS2$comGeneID)), "resFiles/meanEditLevel_nearestGenesList.txt",
            row.names = F, col.names = F, quote = F)
write.table(data.frame(unique(annotGWAS2$comTxID)), "resFiles/meanEditLevel_nearestTranscriptList.txt",
            row.names = F, col.names = F, quote = F)

############################################################################
annotGWAS2 <- read.csv("resFiles/meanEditLevelGWAS_mappedGenes.csv")
gwas2 <- merge(annotGWAS2[, c(2,21)], gwas, by = "SNP", all.y = T)
gwas2$SNP <- ifelse(!is.na(gwas2$nearestGene.symbol), gwas2$nearestGene.symbol, gwas2$SNP)

genes_to_annotate <- c("ADAR","ADARB1")

toLabel <- gwas2[grepl("AD", gwas2$SNP), ]
png("assocRes/manPlots/meanEditLevel_manhattanPlot.png", width=1200, height=600)
manhattan(gwas2, main="",
          suggestiveline = -log10(5e-05),
          highlight = genes_to_annotate, 
          annotateTop = F)

text(
  x = toLabel$CHR,
  y = -log10(toLabel$P),
  labels = toLabel$SNP,
  pos = 3,      # Place text above the point
  cex = 2,    # !!! CUSTOM TEXT SIZE !!!
  font = 2,     # Bold font
  col = "blue4" # Custom color
)
dev.off()






#1. Calculate the necessary coordinates for the X-axis
gwas_with_coords <- gwas2 %>% 
  # Compute chromosome lengths
  group_by(CHR) %>% 
  summarise(chr_len = max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  # Join this info back to the original dataset
  left_join(gwas2, ., by=c("CHR"="CHR")) %>%
  # Calculate the final X-coordinate (BP within chromosome + cumulative length of previous chromosomes)
  arrange(CHR, BP) %>%
  mutate(X_COORD = BP + tot) %>%
  # Calculate the Y-coordinate
  mutate(Y_COORD = -log10(P))


# Identify the rows you want to label (based on your criteria: SNP contains "AD")
toLabel <- gwas_with_coords %>%
  filter(grepl("AD", SNP))

# Create the vector for highlighting (if you want the points colored differently)
genes_to_annotate <- toLabel$SNP


# Calculate Y-axis maximum height, adding extra space (e.g., 1.5) for the labels
y_max_value <- max(gwas_with_coords$Y_COORD, na.rm = TRUE) + 1.5

# Start PNG device
png("assocRes/manPlots/meanEditLevel_manhattanPlot.png", width=1200, height=600)

# Generate the Manhattan plot using the data with coordinates
manhattan(gwas_with_coords, 
          main="",
          suggestiveline = -log10(5e-05),
          highlight = genes_to_annotate, 
          annotateHighlight = FALSE, # Make sure this is FALSE or not present
          annotatePval = NULL,       # Make sure this is NULL or not present
          ylim = c(0, y_max_value)   # Ensure enough vertical space for labels
)

# Add the manual text labels using the CORRECT coordinates
text(
  x = toLabel$X_COORD,  # <<< Use the calculated X coordinate
  y = toLabel$Y_COORD,  # <<< Use the calculated Y coordinate
  labels = toLabel$SNP,
  pos = 3,      # Place text above the point
  cex = 2,      # Large text size
  font = 2,     # Bold font
  col = "blue4",# Custom color
  offset = 0.5  # Space between the point and the text
)

# Close the PNG device
dev.off()




############################################################################
# Load data
library(tidyverse)
library(ggfastman)
gwas <- data.frame(fread("assocRes/meanEditLevelGWAS.mlma"))
gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                 c("chr", "pos", "REF", "BETA", "SE", "pvalue"), skip_absent = T)
gwas$chr <- paste0("chr", gwas$chr)


fast_manhattan(gwas, build='hg38', speed = "f") +
  # add significance line
  geom_hline( aes(yintercept = -log10(5e-05)), color ="deeppink") + 
  coord_flip() +
  scale_x_discrete(position = "top") +
  scale_y_reverse() +
  theme_bw()



fast_manhattan(gwas, build='hg38', speed = "f") +
  geom_hline(aes(yintercept = -log10(5e-05)), color ="deeppink") + 
  coord_flip() +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_reverse(expand = c(0,0)) +
  labs(
    y = expression(-log[10](italic(p-value))),  # This will be on top after flip
    x = "Chromosome"  # This will be on left after flip
  ) +
  theme_bw()


# Store the original plot
p <- fast_manhattan(gwas, build='hg38', speed = "f") +
  geom_hline(aes(yintercept = -log10(5e-05)), color ="deeppink") +
  theme_bw()

# Get the built plot data to extract axis info
pb <- ggplot_build(p)

# Apply transformations with preserved axis info
p +
  coord_flip() +
  scale_x_discrete(
    position = "top",
    breaks = pb$layout$panel_params[[1]]$x$get_breaks(),
    labels = pb$layout$panel_params[[1]]$x$get_labels()
  ) +
  scale_y_continuous(  # Note: might need continuous instead of reverse
    trans = "reverse",
    breaks = pb$layout$panel_params[[1]]$y$get_breaks(),
    labels = pb$layout$panel_params[[1]]$y$get_labels()
  )





#############################################################################
expr <- data.frame(fread("../WGCNA2/trim_VSMC_Kallisto_NormExpr_All.csv"))
gwasExpr <- expr[expr$external_gene_name %in% mg$nearestGene.symbol, ]
rownames(gwasExpr) <- gwasExpr$external_gene_name
gwasExpr <- data.frame(t(gwasExpr[, -c(1:3)]))
gwasExpr <- log(gwasExpr + 1)
gwasExpr[1:10, 1:15]

gwasExprLong <- gwasExpr %>%
  pivot_longer(cols = everything(), names_to = "Gene", values_to = "Expression")

# Make the boxplot
gwasExprPlot <- ggplot(gwasExprLong, aes(x = Gene, y = Expression)) +
  geom_boxplot(outlier.size = 0.1, size = 0.1, fill = "lightblue") +
  labs(
    x = "Nearest gene to lead GWAS variant",
    y = bquote("log"[10]~"TPM"),
    title = ""
  ) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("gwasExprPlot.png", gwasExprPlot,
       path = "Figures",
       scale = 1,
       width = 80,
       height = 50,
       units = "mm",
       dpi = 300,
       limitsize = TRUE)

gwasExpr_filtered <- gwasExpr[rowSums(gwasExpr) != 0, ]
gwasExpr_filtered[, 1:10]
