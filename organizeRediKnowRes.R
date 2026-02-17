library(data.table)

# Functions ---------------------------------------------------------------
# Define the merge function for Reduce
mergeFunc <- function(x, y) {
  merge(x, y, by = "coordID", all = TRUE) 
}

# Specify target directory and files --------------------------------------
tagDir <- "geuvadisRediKnown/"
tagFiles <- list.files(tagDir, "outTable", recursive = T)
tagID <- sub("/.*", "", tagFiles)

# Import and select only Frequency  ---------------------------------------
emptyFiles <- character(0)

freqList <- lapply(tagFiles, function(x) {
  df_path <- paste0(tagDir, x)
  df <- tryCatch(fread(df_path), error = function(e) return(NULL))
  df <- as.data.frame(df)
  
  # Check for problems: NULL, empty, or missing required columns
  if (is.null(df) || nrow(df) == 0 || 
      !all(c("Region", "Position", "Reference", "Frequency") %in% colnames(df))) {
    # Append filename to emptyFiles (using <<- to modify global variable)
    emptyFiles <<- c(emptyFiles, x)
    return(NULL)
  }else {
  
  df$coordID <- paste0(df$Region, ":", df$Position)
  df <- df[df$`Coverage-q30` >= 20  & df$Frequency > 0, ]
  df <- df[, c("coordID", "Frequency")]
  return(df)
  }
})
names(freqList) <- tagID


# aaa remove empty dataframes, rename to samples and merge to single dataframe -------
freqList <- Filter(function(df) !is.null(df) && nrow(df) > 0, freqList)

# export list of samples for complete samples 
cID <- names(freqList)
write.table(data.frame(cID), "geuvadis_Scripts/geuvadisRediKnownComplete.txt", row.names = F, col.names = F, quote = F)

freqList <- lapply(names(freqList), function(x){
  df <- freqList[[x]]
  df <- df[df$Frequency > 0, ]
  colnames(df)[2] <- x
  return(df)
})

gc()
# Use Reduce to iteratively merge the dataframes in the list
freqDF <- Reduce(mergeFunc, freqList)

fwrite(freqDF, "phenoFiles/rawGeuvadisRediKnown.txt.gz", row.names = F)

#########################################################################################################
# Samples with genotype data
fam <- read.table("genoFiles/all_hg38.fam")

# Import raw editing data
freqDF <- data.frame(fread("phenoFiles/rawGeuvadisRediKnown.txt.gz"))
rownames(freqDF) <- freqDF$coordID
freqDF$coordID <- NULL

# Filter out samples without geno data
matchedCols <- match(fam$V2, colnames(freqDF))
matchedCols <- matchedCols[!is.na(matchedCols)]
freqDF <- freqDF[, matchedCols]
freqDF[1:10, 1:10]

# freqDF_noNA <- na.omit(freqDF1)
# freqDF_noNA[, 1:10]


# Count missing values per row
rowNACount <- rowSums(is.na(freqDF))
# Count missing values per column
colNACount <- colSums(is.na(freqDF))


# Filter out rows with more than 50% missing values
freqDF1 <- freqDF[rowNACount <= 0.5 * ncol(freqDF), ]

# Filter out columns with excessive NAs:
freqDF1 <- freqDF1[, colSums(is.na(freqDF1)) <= 0.5 * nrow(freqDF1)]

# Preview the cleaned data (first 20 rows and columns)
freqDF1[1:20, 1:20]

freqDF1 <- cbind(coordID = rownames(freqDF1), freqDF1)

fwrite(freqDF1, "phenoFiles/filteredGeuvadisRediKnown.txt.gz", row.names = F)
