

module load plink2
module load bcftools
module load samtools
module load tabix

# Extract genotype data
plink2 --bfile ../colocalization/1000Genome/1kGMerge --keep genoFiles/matched_EUR_Geno_Pheno.txt --make-bed --out genoFiles/eurGeno

# QC genotype data
plink2 --bfile genoFiles/eurGeno --geno 0.05 --hwe 1e-06 --maf 0.01 --recode vcf --out genoFiles/eurGeno 

bcftools view genoFiles/eurGeno.vcf --threads 6 -Oz -o genoFiles/eurGeno.vcf.gz
bcftools annotate -x INFO genoFiles/eurGeno.vcf.gz | bcftools +fill-tags  > genoFiles/eurGeno_QCed_tem.vcf

bcftools query -l genoFiles/eurGeno.vcf.gz > genoFiles/eurGeno_samples.txt

bcftools view genoFiles/eurGeno_QCed_tem.vcf -S genoFiles/eurGeno_samples.txt | sed 's/chr//g' | bgzip > genoFiles/eurGeno_QCed.vcf.gz
bcftools index -t genoFiles/eurGeno_QCed.vcf.gz
rm genoFiles/eurGeno_QCed_tem.vcf

# Prep coveriates
plink2 --vcf genoFiles/eurGeno.vcf.gz --pca --out covFiles/eurGeno  
Rscript geuvadis_Scripts/autoPrepGenoPCA_Eur.R

# Prep gene expression bed file ---------------------------
Rscript geuvadis_Scripts/prepPheno_EUR.R

# bgzip phenoFiles/geneExpr_qtlTools_Ready.bed
# tabix -f phenoFiles/geneExpr_qtlTools_Ready.bed.gz

# Prep coveriates
plink2 --vcf genoFiles/eurGeno.vcf.gz --pca --out covFiles/eurGeno  
Rscript geuvadis_Scripts/autoPrepGenoPCA_Eur.R
Rscript geuvadis_Scripts/prepCov_EUR.R


for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/eurGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/isoformQuant_EUR_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/isoformQuant_EUR_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run QTLtools
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/eurGeno_Sex_10pc.txt \
        --window 50000 \
        --std-err \
        --grp-best \
        --out geuvadisEUR_cisEQTL/chr${chr}_geuvadisEUR_cisEQTL.txt.gz \
        --permute 100

done


QTLtools pca --bed phenoFiles/isoformQuant_EUR_qtlTools_Ready.bed.gz --scale --center --out covFiles/isoformQuant_EUR_qtlTools_Ready

#### TODO
## A script that calls other scripts for prepping geno and each phenotypes
