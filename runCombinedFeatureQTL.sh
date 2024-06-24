
module load plink2
module load bcftools
module load samtools
module load tabix

Rscript geuvadis_Scripts/combinePhenoCov.R

for chr in {1..22}
do
    echo "Processing chromosome ${chr}"
    # vcf files
    bcftools view genoFiles/yriGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz
    
    # exp files
    zgrep -w "^#chr" phenoFiles/isorform_LCC_combinedPheno_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/isorform_LCC_combinedPheno_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed
    
    QTLtools cis --vcf genoFiles/chr${chr}.vcf.gz --bed phenoFiles/chr${chr}.bed.gz --cov covFiles/geuYri_Sex_10GenoPC_10SplicePC.txt --window 50000 --std-err --grp-best --out geuYri_IsoLcc_cisQTL/chr${chr}_geuYri_IsoLcc_cisQTL.txt.gz --permute 1000

done
