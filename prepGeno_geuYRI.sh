

module load plink2
module load bcftools
module load samtools
module load tabix


# Extract genotype data
plink2 --bfile ../colocalization/1000Genome/yorubaSamps1kGMerge --keep genoFiles/matched_YRI_Geno_Pheno.txt --make-bed --out genoFiles/yriGeno

# QC genotype data
plink2 --bfile genoFiles/yriGeno --geno 0.05 --hwe 1e-06 --maf 0.01 --recode vcf --out genoFiles/yriGeno 

bcftools view genoFiles/yriGeno.vcf --threads 6 -Oz -o genoFiles/yriGeno.vcf.gz
bcftools annotate -x INFO genoFiles/yriGeno.vcf.gz | bcftools +fill-tags  > genoFiles/yriGeno_QCed_tem.vcf

bcftools query -l genoFiles/yriGeno.vcf.gz > genoFiles/yriGeno_samples.txt

bcftools view genoFiles/yriGeno_QCed_tem.vcf -S genoFiles/yriGeno_samples.txt | sed 's/chr//g' | bgzip > genoFiles/yriGeno_QCed.vcf.gz
bcftools index -t genoFiles/yriGeno_QCed.vcf.gz
rm genoFiles/yriGeno_QCed_tem.vcf

# Prep gene expression bed file ---------------------------
Rscript geuvadis_Scripts/prepPheno.R

bgzip phenoFiles/geneExpr_qtlTools_Ready.bed
tabix -f phenoFiles/geneExpr_qtlTools_Ready.bed.gz

# Prep coveriates
plink2 --vcf genoFiles/yriGeno.vcf.gz --pca --out covFiles/yriGeno  
Rscript geuvadis_Scripts/autoPrepGenoPCA.R
Rscript geuvadis_Scripts/prepCov.R


for chr in {21..22}
do
    echo "Processing chromosome ${chr}"
    # vcf files
    bcftools view genoFiles/yriGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz
    
    # exp files
    zgrep -w "^#chr" phenoFiles/geneExpr_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneExpr_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed
    
    QTLtools cis --vcf genoFiles/chr${chr}.vcf.gz --bed phenoFiles/chr${chr}.bed.gz --cov covFiles/yriGeno_Sex_10pc.txt --window 50000 --std-err --grp-best --out geuvadisYRI_cisEQTL/chr${chr}_geuvadisYRI_cisEQTL.txt.gz --permute 100

done
