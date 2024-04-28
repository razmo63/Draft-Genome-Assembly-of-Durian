#!/bin/bash

#SBATCH -A uppmax2024-2-7            # Project ID
#SBATCH -M snowy                     # Specify the machine (optional)
#SBATCH -p core                      # Specify the partition or queue
#SBATCH -n 2                         # Number of cores
#SBATCH -t 01:00:00                  # Walltime
#SBATCH -J quast_quality_assessment  # Job name
#SBATCH --mail-type=ALL              # Email notification
#SBATCH --mail-user ramo1628@student.uu.se   # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/assembly_quality_assessment/after_pilon/quality_assessment_before_pilon.%j.out

# Load QUAST module
module load bioinfo-tools
module load quast/5.0.2

# Set path to the assembly before Pilon
input_assembly="/home/ramo1628/genome_project/outputs/assembly_refinement/refined_assembly.fasta"
output_directory="/home/ramo1628/genome_project/outputs/assembly_quality_assessment/after_pilon"

# Run QUAST for assembly quality assessment
quast.py $input_assembly -o "$output_directory"
