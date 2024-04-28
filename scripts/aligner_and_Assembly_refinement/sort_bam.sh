#!/bin/bash

#SBATCH -A uppmax2024-2-7            # Project ID
#SBATCH -M snowy                     # Specify the machine (optional)
#SBATCH -p core                      # Specify the partition or queue
#SBATCH -n 2                         # Number of cores
#SBATCH -t 01:00:00                  # Walltime
#SBATCH -J sort_bam                  # Job name
#SBATCH --mail-type=ALL              # Email notification
#SBATCH --mail-user ramo1628@student.uu.se   # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/genome_alignment/sorted_bam.%j.out  # Output file

# Load samtools module
module load bioinfo-tools            # Load Bioinformatics tools module
module load samtools

# Set input and output paths
input_bam="/home/ramo1628/genome_project/outputs/genome_alignment/aligned_genome.bam"     # Path to the input BAM file
output_dir="/home/ramo1628/genome_project/outputs/genome_alignment"     # Output directory
sorted_bam="${output_dir}/aligned_genome.sorted.bam"               # Path to the sorted BAM file

# Sort the BAM file
samtools sort -o "$sorted_bam" "$input_bam"

# Index the sorted BAM file
samtools index "$sorted_bam"

# Move the index file to the same directory as the sorted BAM file
mv "${sorted_bam}.bai" "${output_dir}/"

