
module load plink
module load plink2

# Extract and QC genotype data
plink2 --bfile ../colocalization/1000Genome/1kGMerge \
  --geno 0.05 \
  --hwe 1e-06 \
  --maf 0.01 \
  --keep genoFiles/matched_EUR_Geno_Pheno.txt \
  --make-bed \
  --out genoFiles/allGeuvadisSampGeno \
  --allow-extra-chr

Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/fixFamSex.R

# Prep covariates
plink --bfile genoFiles/allGeuvadisSampGeno \
  --pca \
  --out covFiles/allGeuvadisSampGeno \
  --allow-extra-chr

Rscript geuvadis_Scripts/autoPrepGenoPCA_allSamp.R
################################################################################


# Extract and QC genotype data
plink2 --bfile genoFiles/all_hg38 \
  --geno 0.05 \
  --hwe 1e-06 \
  --maf 0.01 \
  --keep genoFiles/matched_allSamps_Geno_Pheno2.txt \
  --make-bed \
  --out genoFiles/allGeuvadisSampGeno \
  --allow-extra-chr \
  --autosome \
  --snps-only
  
# Prep covariates
plink --bfile genoFiles/allGeuvadisSampGeno \
  --pca \
  --out covFiles/allGeuvadisSampGeno \
  --allow-extra-chr

Rscript geuvadis_Scripts/autoPrepGenoPCA_allSamp.R
