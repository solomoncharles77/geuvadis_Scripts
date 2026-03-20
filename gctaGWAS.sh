#!/bin/bash

#SBATCH --job-name=gctaAssoc2
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR


module load gcta
module load plink
module load plink2
module load R

gcta --bfile genoFiles/allGeuvadisSampGeno --make-grm --out assocRes/allGeuvadisSampGeno_grm

# A2G GWAS with GCTA
gcta --mlma \
  --bfile genoFiles/allGeuvadisSampGeno \
  --grm assocRes/allGeuvadisSampGeno_grm \
  --pheno phenoFiles/geuvadisMMSitesAEI_Normalized_AssocReady.txt \
  --mpheno 2 \
  --out assocRes2/geuvadisAEI_A2I \
  --qcovar covFiles/geuvadis_Sex_3gpc_covariates.txt \
  --thread-num 28

###################################################
# CEI GWAS with GCTA
gcta --mlma \
  --bfile genoFiles/allGeuvadisSampGeno \
  --grm assocRes/allGeuvadisSampGeno_grm \
  --pheno phenoFiles/geuvadisMMSitesCEI_Normalized_AssocReady.txt \
  --mpheno 2 \
  --out assocRes2/geuvadisCEI_A2I \
  --qcovar covFiles/geuvadis_Sex_3gpc_covariates.txt \
  --thread-num 28


####################################################################################
# mean editing Level GWAS with GCTA
gcta --mlma \
  --bfile genoFiles/allGeuvadisSampGeno \
  --grm assocRes/allGeuvadisSampGeno_grm \
  --pheno phenoFiles/geuvadisMEL_Normalized_AssocReady_wtHeader.txt \
  --mpheno 1 \
  --out assocRes2/geuvadisMEL_A2I \
  --qcovar covFiles/geuvadis_Sex_3gpc_covariates.txt \
  --thread-num 28


####################################################################################
# coding Sequence editingLevel GWAS with GCTA
gcta --mlma \
  --bfile genoFiles/allGeuvadisSampGeno \
  --grm assocRes/allGeuvadisSampGeno_grm \
  --pheno phenoFiles/geuvadisMcEL_Normalized_AssocReady_wtHeader.txt \
  --mpheno 1 \
  --out assocRes2/geuvadisMcEL_A2I \
  --qcovar covFiles/geuvadis_Sex_3gpc_covariates.txt \
  --thread-num 28

#################################################################################

# module load R
Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/processGCTAassocRes.R
