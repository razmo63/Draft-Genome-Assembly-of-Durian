#!/bin/bash -l

# SLURM directives
#SBATCH -A uppmax2024-2-7          # Project ID
#SBATCH -M snowy                   # Specify the machine (optional)
#SBATCH -p core                    # Specify the partition or queue
#SBATCH -n 4                       # Number of nodes
#SBATCH -t 2:00:00                # Walltime
#SBATCH -J braker_annotation       # Job name
#SBATCH --mail-type=ALL            # Email notification
#SBATCH --mail-user ramo1628@student.uu.se # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/gene_annotation_braker/braker_annotation_%j.out  # Output file for Slurm job logs

# Load necessary modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

# Copy GeneMark key file to home directory
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

# Set up AUGUSTUS_CONFIG_PATH and other paths
source $AUGUSTUS_CONFIG_COPY                   # Source the Augustus config copy

chmod a+w -R /home/ramo1628/genome_project/outputs/gene_annotation_braker/augustus_config/species/

# Path to the genome file
GENOME_FILE="/home/ramo1628/genome_project/outputs/repeat_masking/refined_assembly.fasta.masked"

# Path to the genome file
GENOME_FILE="/home/ramo1628/genome_project/outputs/repeat_masking/refined_assembly.fasta.masked"

OUTPUT_DIR="/home/ramo1628/genome_project/outputs/gene_annotation_braker/annotation_scaffold_06"


braker.pl --genome="$GENOME_FILE" \
              --bam /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6040092_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6040093_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6040094_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6040096_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6040097_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6156066_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6156067_scaffold_06/Aligned.sortedByCoord.out.bam, \
                    /home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/SRR6156069_scaffold_06/Aligned.sortedByCoord.out.bam \
              --cores=8 \
              --species=Durio_zibethinus \
              --workingdir="$OUTPUT_DIR" \
              --useexisting \
              --softmasking \
              --AUGUSTUS_CONFIG_PATH=/home/ramo1628/genome_project/outputs/gene_annotation/augustus_config \
              --AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin \
              --AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts \
              --GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy
