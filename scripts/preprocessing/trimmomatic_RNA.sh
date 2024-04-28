#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J trimmomatic_illumina_RNAseq
#SBATCH --mail-type=ALL 
#SBATCH --mail-user ramo1628@student.uu.se
#SBATCH --output=/home/ramo1628/genome_project/outputs/data_preprocessing/trimmed_RNAseq/trimmed_RNAseq.out

# Load modules
module load bioinfo-tools
module load trimmomatic

# Your commands 
java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar PE\
 -threads 2\
 -phred33\
 -basein /home/ramo1628/genome_project/data/raw_data/illumina/RNA/untrimmed/SRR6040095_scaffold_06.1.fastq.gz\
 -baseout /home/ramo1628/genome_project/outputs/data_preprocessing/trimmed_RNAseq/trimmed_RNAseq.fq.gz\
 ILLUMINACLIP:TruSeq2-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
