#!/bin/bash

#SBATCH --job-name=rnae
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR

# activate env
source /scratch/vasccell/cs806/rnaEditing/RNAE/conda/bin/activate rnae

# run geuvadisRNASeq
/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/RNAEditingIndex \
  -d /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads \
  -o /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/AEI/geuvadisRNASeqOUT \
  -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/AEI/geuvadisRNASeqLOG \
  -os /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/AEI/geuvadisRNASeqSUM \
  -f .bam \
  --genome hg38 \
  --verbose \
  --paired
  

# run geuvadisRNASeq CEI
/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/RNAEditingIndex \
  -d /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads \
  -o /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/CEI/geuvadisRNASeqOUT \
  -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/CEI/geuvadisRNASeqLOG \
  -os /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/CEI/geuvadisRNASeqSUM \
  -f .bam \
  -rb /lustre/alice3/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/oppRep3pHg38.bed.gz \
  --genome hg38 \
  --verbose \
  --paired

# run geuvadisRNASeq Tandem
/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/RNAEditingIndex \
  -d /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads \
  -o /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/tandem/geuvadisRNASeqOUT \
  -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/tandem/geuvadisRNASeqLOG \
  -os /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/tandem/geuvadisRNASeqSUM \
  -f .bam \
  -rb /lustre/alice3/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/tandRep3pHg38.bed.gz \
  --genome hg38 \
  --verbose \
  --paired

