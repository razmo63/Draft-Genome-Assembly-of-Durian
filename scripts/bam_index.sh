#!/bin/bash
#SBATCH -A uppmax2024-2-7          # Project ID
#SBATCH -M snowy                   # Specify the machine (optional)
#SBATCH -p core                    # Specify the partition or queue
#SBATCH -n 4                       # Number of cores
#SBATCH -t 4:00:00                 # Walltime
#SBATCH -J bam_index               # Job name
#SBATCH --mail-type=ALL            # Email notification
#SBATCH --mail-user ramo1628@student.uu.se    # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/bam_index_%j.out  # Output file for SLURM job logs

# Load the samtools module
module load bioinfo-tools
module load samtools/1.19

# Set the directory containing sample directories
ALIGNMENT_DIR="/home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads"

# Navigate to the aligned reads directory
cd "$ALIGNMENT_DIR" || exit

# Loop through each sample directory
for SAMPLE_DIR in */; do
    # Check if the sample directory is a directory
    if [ -d "$SAMPLE_DIR" ]; then
        # Navigate to the sample directory
        cd "$SAMPLE_DIR" || exit

        # Get the BAM file path
        BAM_FILE="Aligned.sortedByCoord.out.bam"

        # Check if the BAM file exists
        if [ -f "$BAM_FILE" ]; then
            # Create the index file (.bam.bai) using samtools
            samtools index "$BAM_FILE"
            echo "Index file created for $BAM_FILE"
        else
            echo "BAM file $BAM_FILE not found in $SAMPLE_DIR."
        fi

        # Navigate back to the aligned reads directory
        cd "$ALIGNMENT_DIR" || exit
    fi
done

echo "Indexing completed for all BAM files."
