#!/bin/bash

#SBATCH -A uppmax2024-2-7            # Project ID
#SBATCH -M snowy                     # Specify the machine (optional)
#SBATCH -p core                      # Specify the partition or queue
#SBATCH -n 4                         # Number of cores
#SBATCH -t 5:00:00                   # Walltime
#SBATCH -J repeat_masking            # Job name
#SBATCH --mail-type=ALL              # Email notification
#SBATCH --mail-user ramo1628@student.uu.se      # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/repeat_masking/repeat_masking.%j.out

# Load RepeatMasker module
module load bioinfo-tools
module load RepeatMasker

# Set path to the improved draft genome assembly after Pilon
# Set path to the improved draft genome assembly after Pilon
input_assembly="/home/ramo1628/genome_project/outputs/assembly_refinement/refined_assembly.fasta"
output_dir="/home/ramo1628/genome_project/outputs/repeat_masking"

# Run RepeatMasker for repeat masking
# Here, we specify the output directory using -dir
# We use -gff to generate output in GFF format
# -pa specifies the number of processors to use
# -e specifies the search engine (NCBI in this case)
# -s performs softmasking
RepeatMasker -dir $output_dir -gff -pa 4 -e ncbi -s $input_assembly
