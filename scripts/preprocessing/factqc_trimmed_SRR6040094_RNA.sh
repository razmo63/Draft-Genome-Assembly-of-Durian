#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:20:00
#SBATCH -J FastQC_illumina_RNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user raziyeh.mohseni.1628@student.uu.se
#SBATCH --output=/home/ramo1628/genome_project/outputs/data_preprocessing/trimmed_RNAseq/SRR6040094_trimmed/trimmed_RNA_fastQC.%j.out

# Load modules
module load bioinfo-tools
module load FastQC

# Your commands
fastqc /home/ramo1628/genome_project/data/raw_data/illumina/RNA/trimmed/SRR6040094_scaffold_06* -o /home/ramo1628/genome_project/outputs/data_preprocessing/trimmed_RNAseq/SRR6040094_trimmed
