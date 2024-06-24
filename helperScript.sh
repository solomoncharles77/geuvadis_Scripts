

##################################### 
# Split YRI VCF to one sample per file
for file in yriGeno.vcf; do
  for sample in `bcftools query -l $file`; do
    bcftools view -c1 -Ov -s $sample -o splitYriGeno/${sample}.vcf $file
  done
done


# rename fastq files to matcht genoIDs
while IFS=$' ' read -r orig new; do
    rename -vo "$orig" "$new" *.fastq.gz
done < ../exprFiles/matched_YRI_Pheno.txt

# List of 1000genome sample IDs
awk '{print $2}' matched_YRI_Pheno.txt > ../geuvadis_Scripts/matchedSampIDs.txt


##################################### 
# Split EUR VCF to one sample per file
for file in eurGeno.vcf; do
  for sample in `bcftools query -l $file`; do
    bcftools view -c1 -Ov -s $sample -o splitEurGeno/${sample}.vcf $file
  done
done


# rename fastq files to matcht genoIDs
while IFS=$' ' read -r orig new; do
    rename -vo "$orig" "$new" *.fastq.gz
done < ../exprFiles/matched_EUR_Pheno.txt

# List of 1000genome sample IDs
awk '{print $2}' matched_EUR_Pheno.txt > ../geuvadis_Scripts/matchedSampIDs_EUR.txt


# Extract junction files
for bamfile in `ls *.bam`; do
    echo Converting $bamfile to $bamfile.junc
    #samtools index $bamfile
    /home/c/cs806/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 -s XS $bamfile -o $bamfile.junc
    echo $bamfile.junc >> test_juncfiles.txt
done

python /home/c/cs806/leafcutter/clustering/leafcutter_cluster_regtools.py -j test_juncfiles.txt -m 50 -o testJunc -l 500000

pip install scikit-learn

python /home/c/cs806/leafcutter/scripts/prepare_phenotype_table.py testJunc_perind.counts.gz -p 10
