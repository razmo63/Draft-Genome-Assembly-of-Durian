#!/bin/bash

#SBATCH -A uppmax2024-2-7            # Project ID
#SBATCH -M snowy                     # Specify the machine (optional)
#SBATCH -p core                      # Specify the partition or queue
#SBATCH -n 4                         # Number of cores
#SBATCH -t 4:00:00                  # Walltime
#SBATCH -J assembly_refinement       # Job name
#SBATCH --mail-type=ALL              # Email notification
#SBATCH --mail-user ramo1628@student.uu.se                       # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/assembly_refinement/assembly_refinement.%j.out  # Output file

# Load necessary modules
module load bioinfo-tools            # Load Bioinformatics tools module
module load Pilon/1.24               # Load Pilon module for assembly refinement

# Set paths
genome_assembly="/home/ramo1628/genome_project/outputs/genome_assembly/assembly.fasta"  # Path to the genome assembly file
bam_file="/home/ramo1628/genome_project/outputs/genome_alignment/aligned_genome.sorted.bam"  # Path to the BAM file generated from the alignment

# Run Pilon to refine the assembly
java -jar $PILON_HOME/pilon.jar --genome "$genome_assembly" \
    --bam "$bam_file" \
    --output /home/ramo1628/genome_project/outputs/assembly_refinement/refined_assembly \
    --vcf \
    --tracks \
    --verbose \
    --threads 4
