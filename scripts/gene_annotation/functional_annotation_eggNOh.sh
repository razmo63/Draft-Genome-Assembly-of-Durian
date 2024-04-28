#!/bin/bash -l

# SLURM directives
#SBATCH -A uppmax2024-2-7          # Project ID
#SBATCH -M snowy                   # Specify the machine (optional)
#SBATCH -p core                    # Specify the partition or queue
#SBATCH -n 4                       # Number of nodes
#SBATCH -t 2:00:00                # Walltime
#SBATCH -J eggNOG_mapper           # Job name
#SBATCH --mail-type=ALL            # Email notification
#SBATCH --mail-user ramo1628@student.uu.se # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/functional_annotation_eggNOG/eggNOG_%j.out  # Output file for Slurm job logs

# Load necessary modules
module load bioinfo-tools
module load samtools
module load eggNOG-mapper/2.1.9

# Define input and output file paths
input_file="/home/ramo1628/genome_project/outputs/gene_annotation_braker/protein_sequences.fasta"
output_dir="/home/ramo1628/genome_project/outputs/functional_annotation_eggNOG"
gff_file="/home/ramo1628/genome_project/outputs/gene_annotation_braker/annotation_scaffold_06/augustus.hints.gff"

# Run eggNOG-mapper command with gene prediction using AUGUSTUS and functional annotation using HMMER
emapper.py -m diamond --itype proteins --cpu 4 -i $input_file -o $output_dir --decorate_gff $gff_file --decorate_gff_ID_field GeneID


