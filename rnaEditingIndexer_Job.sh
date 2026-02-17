#!/bin/bash

#SBATCH --job-name=rnae
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=80G
#SBATCH --time=168:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR

# activate env
source /scratch/vasccell/cs806/rnaEditing/RNAE/conda/bin/activate rnae

# run huvecRNASeq
/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/RNAEditingIndex \
  -d /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads \
  -o /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/AEI/huvecRNASeqOUT \
  -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/AEI/huvecRNASeqLOG \
  -os /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/AEI/huvecRNASeqSUM \
  -f .bam \
  --genome hg38 \
  --verbose \
  --paired
  

# run huvecRNASeq CEI
/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/RNAEditingIndex \
  -d /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads \
  -o /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/CEI/huvecRNASeqOUT \
  -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/CEI/huvecRNASeqLOG \
  -os /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/CEI/huvecRNASeqSUM \
  -f .bam \
  -rb /lustre/alice3/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/oppRep3pHg38.bed.gz \
  --genome hg38 \
  --verbose \
  --paired

# run huvecRNASeq Tandem
/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/RNAEditingIndex \
  -d /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads \
  -o /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/tandem/huvecRNASeqOUT \
  -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/tandem/huvecRNASeqLOG \
  -os /lustre/alice3/scratch/vasccell/cs806/geuvadis/RNAEditingIndexer/tandem/huvecRNASeqSUM \
  -f .bam \
  -rb /lustre/alice3/scratch/vasccell/cs806/rnaEditing/RNAE/bin/AEI/RNAEditingIndexer/Resources/Regions/HomoSapiens/tandRep3pHg38.bed.gz \
  --genome hg38 \
  --verbose \
  --paired

