#!/bin/bash

#SBATCH --job-name=REDItools
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=80G
#SBATCH --time=140:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR

module purge
source ~/miniforge3/bin/activate py27
module load R
# Run REDItoolsKnown
while IFS= read -r sample; do
        echo ${sample}
        mkdir huvecRediKnown/${sample}
        REDItoolKnown.py -i mappedReads/${sample}.sorted.bam \
        -f /scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        -l /scratch/vasccell/cs806/rnaEditing/editDB/knownEditingSites_TABLE1_hg38_v3.tab.gz \
        -o geuvadisRediKnown/${sample}/ \
        -t 28

done < /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSeq2.txt

# Organze edit level
Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/organizeKnownEditLevel.R
Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/prepMEL4gwas.R
Rscript /scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/prepMcEL4gwas.R

# # Run REDItoolsNovel
# while IFS= read -r sample; do
#         echo ${sample}
#         mkdir rediNovel/${sample}
#         REDItoolKnown.py -i bamFiles/${sample}.sorted.bam \
#         -f /scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
#         -l editDB/vsmcAtoINovelEditingSites_Coord.tab.gz \
#         -o rediNovel/${sample}/ \
#         -t 28
#         
# done < /scratch/vasccell/cs806/rnaEditing/rnaEditing_scripts/samList.txt

