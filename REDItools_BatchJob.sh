#!/bin/bash

for i in {00..22}; do
  # Create a unique job name
  job_name="geuREDI_${i}"

  # Set the input file path
  input_file="bamList_${i}.txt"

  # Create a temporary script with the modified parameters
  tmp_script=$(mktemp)
  cat > "$tmp_script" <<EOF
#!/bin/bash

#SBATCH --job-name=${job_name}
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=${job_name}.out
#SBATCH --error=${job_name}.err
#SBATCH --export=NONE

cd \$SLURM_SUBMIT_DIR
module purge
source ~/miniforge3/bin/activate py27

EOF


# Run REDItoolsKnown
  cat >> "$tmp_script" <<EOF
# Run REDItoolsKnown
EOF
  while IFS= read -r sample; do
    cat >> "$tmp_script" <<EOF
        echo ${sample}
        mkdir huvecRediKnown/${sample}
        REDItoolKnown.py -i bamFiles/${sample}.sorted.bam \
        -f /scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        -l editDB/knownEditingSites_TABLE1_hg38_v3.tab.gz \
        -o geuvadisRediKnown/${sample}/ \
        -t 28
EOF
  done < /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/${input_file}

# Print the job script
  echo "Job Script for ${job_name}:"
  cat "$tmp_script" > job2rev.sh
  sbatch "$tmp_script"
  echo ""

  # Remove the temporary script
  rm "$tmp_script"
done

# # Organze edit level
# Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/organizeKnowEditLevel.R
# Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/prepMEL4gwas.R
# Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/prepMcEL4gwas.R