#!/bin/bash


# Open and edit slurm job details
echo "Opening job scripts"
gedit fastqcScript.job
gedit mapReads.job
gedit samtoolsScript.job
gedit countReads.job

# Submit job for read QC
echo "Submitting job for read QC."
sbatch fastqcScript.job
echo ""

# Submit job for read mapping
echo "Submitting job for read mapping."
sbatch mapReads.job
echo ""
echo "Waiting for read mapping to finish."


file_to_check="/lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSeq2.txt"

while [ ! -f "$file_to_check" ]; do
    echo "Waiting for the read mapping file to become available..."
    sleep 1h  # Sleep for 1 hour
done

echo "File is now available! Continuing with the script."

# Wait until the number of sam files matches the expected count
while [ "$(find /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/ -name "*Aligned.out.bam" | wc -l)" -lt "$(wc -l /lustre/alice3/scratch/vasccell/cs806/geuvadis/geuvadis_Scripts/allSeq2.txt | awk '{print $1}')" ]
do
sleep 1h
done

# Get a file list for samtools
find /lustre/alice3/scratch/vasccell/cs806/geuvadis/mappedReads/ -maxdepth 1 -name "*Aligned.out.bam" -exec basename {} \; > allSam.txt

# Sort the files and check if contents match
sort allSeq2.txt > allSeq2sorted.txt
sort allSam.txt > allSamsorted.txt
sed -i "s|Aligned.out.sam||g" allSamsorted.txt

if cmp --silent allSeq2sorted.txt allSamsorted.txt; then
echo "Read mapping completed."
echo ""
echo "Submitting Jobs for Samtools sorting and indexing"
sbatch samtoolsScript.job
echo ""
echo "Submitting Jobs for Rsubread feature count"
sbatch countReads.job
else
  echo "Please check that reads have been properly mapped."
fi

echo ""
echo "Submitting Stringtie assembly job"
sbatch runStrintie.job