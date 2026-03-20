

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table("covFiles/allGeuvadisSampGeno.eigenvec", header = FALSE, skip=0)
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste0("genoPC", c(1:ncol(eigenvec)))
eigenvec <- data.frame(t(eigenvec))
eigenvec <- data.frame(SampleID = rownames(eigenvec), eigenvec)
write.table(eigenvec, "covFiles/allGeuvadisSamp_genoPCs.txt", sep = "\t", quote = F, col.names = colnames(eigenvec))

