
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(gplots)
library(openxlsx)

prep4Densityplot <- function(file_path, feature_name) {
  sumstat <- data.frame(fread(file_path, header = F, stringsAsFactors = F))
  sumstat <- sumstat[sumstat$V18 < 5e-08, c(1, 9)]
  colnames(sumstat) <- c("phenotypeID", "dist")
  
  num_bins <- 10
  
  # Create bins for positive values
  pos <- sumstat %>%
    filter(dist > 0 & dist < 100000) %>%
    mutate(bins = ntile(dist, num_bins)) %>%
    mutate(bins = paste0(bins, "0"))
  
  neg <- sumstat %>%
    filter(dist < 0 & dist > -100000) %>%
    mutate(bins = ntile(dist, num_bins)) %>%
    mutate(bins = paste0("-", bins, "0"))
  
  sumstat <- rbind(pos, neg)
  sumstat$feature <- feature_name
  
  return(sumstat)
}

isoformQuant <- prep4Densityplot("isoformQuant_EUR_cisEQTL/isoformQuant_EUR_cisEQTL_permute_ALL.txt.gz", "isoformQuant")
leafcutter <- prep4Densityplot("geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_permute_ALL.txt.gz", "leafcutter")
exonCount <- prep4Densityplot("geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_permute_ALL.txt.gz", "exonCount")
intronCount <- prep4Densityplot("geneIntronCount_EUR_cisEQTL/geneIntronCount_EUR_cisEQTL_permute_ALL.txt.gz", "intronCount")
exonIntronRatio <- prep4Densityplot("geneExonIntronRatio_EUR_cisEQTL/geneExonIntronRatio_EUR_cisEQTL_permute_ALL.txt.gz", "exonIntronRatio")
comboFeature <- prep4Densityplot("combinedPheno_EUR_cisEQTL/combinedPheno_EUR_cisEQTL_permute_ALL.txt.gz", "comboFeature")

# # Plot single densities ---------------------------------------------------
# leafcutter %>%
#   ggplot(aes(x = dist)) +
#   geom_density(alpha = 0.3) +  # Add transparency for better visualization
#   labs(title = "Density Plot", x = "Values", y = "Density") +
#   scale_fill_discrete(name = "txC")  # Add legend with proper title


# Plot combined densities -------------------------------------------------
rbind(isoformQuant, leafcutter, exonCount, intronCount, exonIntronRatio, comboFeature) %>% 
  group_by(feature) %>%
  ggplot(aes(x = dist, fill = factor(feature))) +
  geom_density(alpha = 0.3) +  # Add transparency for better visualization
  labs(title = "Density Plot", x = "Values", y = "Density") +
  scale_fill_discrete(name = "txC")  # Add legend with proper title


# Make upset ----------------------

dataList <- list(isoformQuant$phenotypeID, leafcutter$phenotypeID, exonCount$phenotypeID,
                   intronCount$phenotypeID, exonIntronRatio$phenotypeID, comboFeature$phenotypeID)
names(dataList) <- c("isoformQuant", "leafcutter", "exonCount", "intronCount", "exonIntronRatio", "comboFeature")

names(dataList) <- c("iso", "LCC", "exC", "inC", "exInRa", "cbF")

str(dataList)

m = make_comb_mat(dataList, mode = "intersect")
ss <- set_size(m)
cs <- comb_size(m)
UpSet(m,set_order = order(set_size(m), decreasing = FALSE),
      top_annotation = HeatmapAnnotation(
  "rGene Intersections" = anno_barplot(cs, 
                                       ylim = c(0, max(cs)*1.1),
                                       border = FALSE, 
                                       gp = gpar(fill = "black"), 
                                       height = unit(4, "cm")
  ), 
  annotation_name_side = "left", 
  annotation_name_rot = 90),
  
  left_annotation = rowAnnotation(
    "rGene per feature" = anno_barplot(-ss, 
                                      baseline = 0,
                                      axis_param = list(
                                        at = c(0, -2000, -4000 ),
                                        labels = c(0, 2000, 4000),
                                        labels_rot = 0),
                                      add_numbers = FALSE,
                                      border = FALSE, 
                                      gp = gpar(fill = "black"), 
                                      width = unit(3, "cm")
    ),
    set_name = anno_text(set_name(m), 
                         location = 0.5, 
                         just = "right",
                         width = max_text_width(set_name(m)) + unit(4, "mm"))
  ), 
  right_annotation = NULL,
  show_row_names = FALSE)

# ggsave("rQTL_UpsetPlot.pdf", path = "outputFilesnPlots",
#        units = "in", width = 5.5, height = 3,       )


# set_name(m)
# comb_name(m)
# set_size(m)
# comb_size(m)
# comb_degree(m)
# t(m)



senVenn <- venn(dataList, show.plot = FALSE)
saveOver <- attr(x = senVenn, "intersections")
names(saveOver) <- sub("*:", "_", names(saveOver))
names(saveOver) <- sub("*:", "_", names(saveOver))
names(saveOver) <- sub("*:", "_", names(saveOver))
names(saveOver) <- sub("*:", "_", names(saveOver))
names(saveOver) <- sub("*:", "_", names(saveOver))

#openxlsx::write.xlsx(saveOver, file = "outputFilesnPlots/OverlapGeneList.xlsx")


lapply(names(saveOver), function(x){
  write.table(saveOver[x], file=paste0("featOverlap/", x, ".txt"), row.names = F, col.names = F, quote = F)
})
  
