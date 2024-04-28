#!/bin/bash

#SBATCH -A uppmax2024-2-7            # Project ID
#SBATCH -M snowy                     # Specify the machine (optional)
#SBATCH -p core                      # Specify the partition or queue
#SBATCH -n 2                         # Number of cores
#SBATCH -t 03:00:00                  # Walltime
#SBATCH -J genome_alignment          # Job name
#SBATCH --mail-type=ALL              # Email notification
#SBATCH --mail-user ramo1628@student.uu.se                  # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/genome_alignment/genome_alignment.%j.out             # Output file

# Load necessary modules
module load bioinfo-tools         # Load Bioinformatics tools module
module load bwa                   # Load BWA module for alignment
module load samtools              # Load samtools for file conversion

# Set paths
genome_assembly="/home/ramo1628/genome_project/outputs/genome_assembly/assembly.fasta"                          # Path to the genome assembly file
illumina_reads_1="/home/ramo1628/genome_project/data/raw_data/illumina/DNA/SRR6058604_scaffold_06.1P.fastq.gz"  # Path to the first Illumina read file
illumina_reads_2="/home/ramo1628/genome_project/data/raw_data/illumina/DNA/SRR6058604_scaffold_06.2P.fastq.gz"  # Path to the second Illumina read file
output_prefix="/home/ramo1628/genome_project/outputs/genome_alignment/aligned_genome"                           # Prefix for the output file

# Index the genome assembly
bwa index "$genome_assembly"            

# Map Illumina reads to the genome assembly and output to a SAM file
bwa mem -t 2 "$genome_assembly" "$illumina_reads_1" "$illumina_reads_2" -o "$output_prefix.sam"

# Convert SAM to BAM format
samtools view -bS "$output_prefix.sam" > "$output_prefix.bam"

# Remove intermediate SAM file
rm "$output_prefix.sam"


