library(GGally)
library(tidyverse)
library(psych)

fam <- read.table("genoFiles/allGeuvadisSampGeno.fam")
fam <- fam[, c(2,5)]
fam$V5 <- ifelse(fam$V5 == 1, "Male", "Female")

aei <- read.table("phenoFiles/geuvadisMMSitesAEI_clean_AssocReady.txt")
aei <- aei[, c(2,4)]
head(aei)

cei <- read.table("phenoFiles/geuvadisMMSitesCEI_clean_AssocReady.txt")
cei <- cei[, c(2,4)]
head(cei)

mei <- read.table("phenoFiles/geuvadisMeanEditingLevel_clean_AssocReady.txt")
mei <- mei[, c(2,3)]
head(mei)

cds <- read.table("phenoFiles/geuvadisCodingSequenceEditingLevel_clean_AssocReady.txt")
cds <- cds[, c(2,3)]
head(cds)

indList <- list(aei, cei, mei, cds, fam)

edInd <- indList %>% reduce(full_join, by='V2')
colnames(edInd) <- c("Sample", "AEI", "CEI", "MEL", "McEL", "Sex")
edInd <- edInd[complete.cases(edInd), ]
head(edInd)
write.csv(edInd, "resFiles/allEditIndex.csv", row.names = F)

png("resPlots/geuvadisEdIndCorPlot.png", width=1200, height=600)
edInd %>% 
  select(AEI, CEI, MEL, McEL) %>% 
  ggpairs(columnLabels = c("AEI", "CEI", "MEL", "McEL"))
dev.off()

png("resPlots/geuvadisEdIndCorPlot2.png", width=1200, height=600)
ggpairs(edInd, 
        columns = 2:5,           # Numeric columns only
        aes(color = Sex, 
            fill = Sex),     # Group coloring if desired
        upper = list(continuous = wrap("cor", method = "pearson")),
        diag = list(continuous = wrap("barDiag", 
                                      bins = 20,           # Bin count
                                      color = "black",     # Border color
                                      alpha = 0.7)),       # Transparency
        lower = list(continuous = wrap("points", 
                                       alpha = 0.5))) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_bw()
dev.off()
