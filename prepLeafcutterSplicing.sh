

conda activate base
module load plink2
module load bcftools
module load bedtools2
module load samtools
module load tabix


# Extract junction files
rm -i mappedReads/geuYriJunc/geuYRIjuncFiles.txt

for bamfile in `ls mappedReads/*Aligned.out.sam.sorted.bam`; do
    echo Converting $bamfile to $bamfile.junc
    sampID=$(basename "${bamfile}.junc")
    #samtools index $bamfile
    /home/c/cs806/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 -s XS $bamfile -o mappedReads/geuYriJunc/${sampID}
    echo $bamfile.junc >> mappedReads/geuYriJunc/geuYRIjuncFiles.txt
done


conda activate base
python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j mappedReads/geuYriJunc/geuYRIjuncFiles.txt -m 50 -o mappedReads/geuYriJunc/juncClus -l 500000

python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py mappedReads/geuYriJunc/juncClus_perind.counts.gz -p 10

# Prep intron splicing bed file #########################
Rscript geuvadis_Scripts/prepPheno_lcs.R

# Prep coveriates #####################
Rscript geuvadis_Scripts/prepCov_lcs.R


for chr in {1..22}
do
    echo "Processing chromosome ${chr}"
    # vcf files
    bcftools view genoFiles/yriGeno_QCed.vcf.gz --regions ${chr} -Oz -o genoFiles/chr${chr}.vcf.gz
    bcftools index -t genoFiles/chr${chr}.vcf.gz
    
    # exp files
    zgrep -w "^#chr" phenoFiles/geuYri_LCsplice_qtlTools_Ready.bed.gz > phenoFiles/chr${chr}.bed
    zgrep -w ^${chr} phenoFiles/geuYri_LCsplice_qtlTools_Ready.bed.gz | sort -nk2 >> phenoFiles/chr${chr}.bed
    cat phenoFiles/chr${chr}.bed | bgzip > phenoFiles/chr${chr}.bed.gz
    tabix -f phenoFiles/chr${chr}.bed.gz
    rm phenoFiles/chr${chr}.bed
    
    QTLtools cis --vcf genoFiles/chr${chr}.vcf.gz --bed phenoFiles/chr${chr}.bed.gz --cov covFiles/geuYri_Sex_10GenoPC_10SplicePC.txt --window 50000 --std-err --grp-best --out geuYri_LCC_cisQTL/chr${chr}_geuYri_LCC_cisQTL.txt.gz --permute 100

done

