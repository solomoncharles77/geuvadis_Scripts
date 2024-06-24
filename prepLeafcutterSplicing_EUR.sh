


conda activate base
module load plink2
module load bedtools2
module load samtools
module load tabix


# Extract junction files
rm -i mappedReads/geuEurJunc/geuEURjuncFiles.txt

for bamfile in `ls mappedReads/*.sorted.bam`; do
    echo Converting $bamfile to $bamfile.junc
    sampID=$(basename "${bamfile}.junc")
    /home/c/cs806/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 -s XS $bamfile -o mappedReads/geuEurJunc/${sampID}
    echo mappedReads/geuEurJunc/${sampID} >> mappedReads/geuEurJunc/geuEURjuncFiles.txt
done


conda activate base
python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j mappedReads/geuEurJunc/geuEURjuncFiles.txt -m 50 -o mappedReads/geuEurJunc/juncClus -l 500000

python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py mappedReads/geuEurJunc/juncClus_perind.counts.gz -p 10

# Prep intron splicing bed file #########################
Rscript geuvadis_Scripts/prep_leafcutterQuant_EUR.R

# Prep coveriates -------------------
# Get pheno PCA
QTLtools pca --bed phenoFiles/geneLeafcutter_EUR_qtlTools_Ready.bed.gz --scale --center --out covFiles/geneLeafcutter_EUR_qtlTools_Ready
# Combine geno and pheno PCA
Rscript geuvadis_Scripts/prep_leafcutterCov_EUR.R


for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    # Process VCF files
    bcftools view genoFiles/eurGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz

    # Process expression files
    zgrep -w "^#chr" phenoFiles/geneLeafcutter_EUR_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geneLeafcutter_EUR_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed

    # Run QTLtools Nominal
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneLeafcutter_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneLeafcutter_EUR_cisEQTL/chr${chr}_geneLeafcutter_EUR_cisEQTL_nominal.txt.gz \
        --nominal 0.05
        
    # Run QTLtools Permute
    QTLtools cis \
        --vcf genoFiles/chr${chr}.vcf.gz \
        --bed phenoFiles/chr${chr}.bed.gz \
        --cov covFiles/geneLeafcutter_geno3pc_pheno50pc.txt \
        --window 1000000 \
        --std-err \
        --grp-best \
        --out geneLeafcutter_EUR_cisEQTL/chr${chr}_geneLeafcutter_EUR_cisEQTL_permute.txt.gz \
        --permute 100

done


zcat geneLeafcutter_EUR_cisEQTL/*_permute.txt.gz | gzip -c > geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_permute_ALL.txt.gz
zcat geneLeafcutter_EUR_cisEQTL/*_nominal.txt.gz | gzip -c > geneLeafcutter_EUR_cisEQTL/geneLeafcutter_EUR_cisEQTL_nominal_ALL.txt.gz

