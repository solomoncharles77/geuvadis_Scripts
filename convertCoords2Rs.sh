#!/bin/bash

#SBATCH --job-name=SNP_Lookup
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR

module load plink
module load plink2


plink2 --bfile genoFiles/all_hg38 \
  --geno 0.05 \
  --hwe 1e-06 \
  --maf 0.01 \
  --keep genoFiles/matched_allSamps_Geno_Pheno2.txt \
  --make-bed \
  --out genoFiles/allGeuvadisSampGeno_RS \
  --allow-extra-chr \
  --autosome \
  --snps-only

awk 'BEGIN { OFS="\t" } { print $2 }' genoFiles/allGeuvadisSampGeno_RS.bim > genoFiles/geuvadis_coorID.txt
grep -wFf genoFiles/geuvadis_coorID.txt /scratch/vasccell/cs806/colocalization/dbSNP/hg38.snp151_All.bed.rsID_CoordID.txt > genoFiles/geuvadis_rsLookup.txt
Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/convertCoords2Rs.R

