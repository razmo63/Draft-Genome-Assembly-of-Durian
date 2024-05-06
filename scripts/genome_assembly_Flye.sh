#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 5:00:00
#SBATCH -J flye_genome_assembly
#SBATCH --mail-type=ALL 
#SBATCH --mail-user ramo1628@student.uu.se
#SBATCH --output=/home/ramo1628/genome_project/outputs/genome_assembly/flye_genome_assembly.%j.out

# Load modules
module load bioinfo-tools
module load Flye

# Your commands 
flye --pacbio-raw /home/ramo1628/genome_project/data/raw_data/pacbio_data/SRR6037732_scaffold_06.fq.gz --out-dir /home/ramo1628/genome_project/outputs/genome_assembly/ --threads 4
