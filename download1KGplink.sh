#!/bin/bash

# Important Note
# The 1KG version downloaded by this script is missing some 
# GEUVADIS mRNA samples that I am able to download in 
# colocalization folder

module load plink
module load plink2

cd genoFiles

# Download 1K genome
pgen=https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/vx09262b4k1kszy/all_hg38.pvar.zst?dl=1
sample=https://www.dropbox.com/s/2e87z6nc4qexjjm/all_hg38.psam?dl=1

wget $pgen
mv 'all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen

wget $pvar
mv 'all_hg38.pvar.zst?dl=1' all_hg38.pvar.zst

wget $sample
mv 'all_hg38.psam?dl=1' all_hg38.psam

# convert 1K genome to plink format
plink2 --pfile all_hg38 vzs\
      --max-alleles 2 \
      --make-bed \
      --allow-extra-chr 0 \
      --out all_hg38

