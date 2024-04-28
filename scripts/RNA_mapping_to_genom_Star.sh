#!/bin/bash

# SLURM job directives
#SBATCH -A uppmax2024-2-7            # Project ID
#SBATCH -M snowy                     # Specify the machine (optional)
#SBATCH -p core                      # Specify the partition or queue
#SBATCH -n 4                         # Number of cores
#SBATCH -t 24:00:00                  # Walltime
#SBATCH -J RNA_alignment             # Job name
#SBATCH --mail-type=ALL              # Email notification
#SBATCH --mail-user ramo1628@student.uu.se   # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/RNA_mapping_to_genome.%j.out

# Load necessary modules
module load bioinfo-tools
module load star/2.7.9a

# Set paths for the STAR genome index and raw RNA-seq data
genome_index="/home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/genome_index"  # Path to the STAR genome index
genome_file="/home/ramo1628/genome_project/outputs/repeat_masking/refined_assembly.fasta.masked"  # Path to the RepeatMasker masked genome file
RNAseq_path="/home/ramo1628/genome_project/data/raw_data/illumina/RNA_trimmed_06"  # Path to the raw RNA-seq data
RNAseq_samples="SRR6040092_scaffold_06 SRR6040093_scaffold_06 SRR6040094_scaffold_06 SRR6040096_scaffold_06 SRR6040097_scaffold_06 SRR6156066_scaffold_06 SRR6156067_scaffold_06 SRR6156069_scaffold_06"

# Create the STAR genome index directory if it doesn't exist
mkdir -p "$genome_index"

# Generate the STAR genome index using the genome FASTA file
# STAR options:
#   --runMode genomeGenerate: specifies the mode for generating the genome index
#   --runThreadN 8: number of threads to use
#   --genomeDir: path to the output directory for the genome index
#   --genomeFastaFiles: path to the genome FASTA file
#   --genomeSAindexNbases: determines the suffix array size for the genome index
STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir "$genome_index" \
     --genomeFastaFiles "$genome_file" \
     --genomeSAindexNbases 11  

# Loop through each RNA-seq sample
for sample in $RNAseq_samples; do
  
    # Create output directory for this RNA-seq sample
    mkdir -p "/home/ramo1628/genome_project/outputs/genome_mapping_RNA/aligned_reads/${sample}"  # create output directory for this RNA-Seq set
    
    # Print the paths to check if they are correctly constructed
    echo "Input files: $RNAseq_path/${sample}.1.fastq $RNAseq_path/${sample}.2.fastq"
    
    # Run STAR alignment for the current sample
    # STAR options:
    #   --readFilesIn: specifies the input FASTQ files for the RNA-seq sample
    #   --outFileNamePrefix: prefix for the output files
    #   --outSAMtype: type of output files (BAM and sorted by coordinate)
    #   --outSAMunmapped: specifies treatment of unmapped reads in the output SAM file
    STAR --runThreadN 8 \
         --genomeDir "$genome_index" \
         --readFilesIn "$RNAseq_path/${sample}.1.fastq" "$RNAseq_path/${sample}.2.fastq" \
         --outFileNamePrefix "/home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads/${sample}/" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within
done

