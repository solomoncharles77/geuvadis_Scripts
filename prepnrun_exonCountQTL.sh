

module load plink2
module load bcftools
module load samtools
module load tabix
module load bedtools2


# Prep gene expression bed file ---------------------------
Rscript geuvadis_Scripts/prep_geneExonCount_EUR.R

# bgzip phenoFiles/geneExpr_qtlTools_Ready.bed
# tabix -f phenoFiles/geneExpr_qtlTools_Ready.bed.gz

# Prep coveriates -------------------
# Get pheno PCA
QTLtools pca --bed phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed.gz --scale --center --out covFiles/geneExonCount_EUR_qtlTools_Ready
# Combine geno and pheno PCA
Rscript geuvadis_Scripts/prep_exonCov_EUR.R


for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/eurGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run QTLtools
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneExonCount_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneExonCount_EUR_cisEQTL/chr${chr}_geneExonCount_EUR_cisEQTL_nominal.txt.gz \
        --nominal 0.05

done



###### permute

for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/eurGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneExonCount_EUR_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run QTLtools
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneExonCount_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneExonCount_EUR_cisEQTL/chr${chr}_geneExonCount_EUR_cisEQTL_permute.txt.gz \
        --permute 100

done

zcat geneExonCount_EUR_cisEQTL/*_permute.txt.gz | gzip -c > geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_permute_ALL.txt.gz
zcat geneExonCount_EUR_cisEQTL/*_nominal.txt.gz | gzip -c > geneExonCount_EUR_cisEQTL/geneExonCount_EUR_cisEQTL_nominal_ALL.txt.gz

