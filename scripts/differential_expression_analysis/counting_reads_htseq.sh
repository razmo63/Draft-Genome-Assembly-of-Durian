#!/bin/bash -l

#SBATCH -A uppmax2024-2-7          # Project ID
#SBATCH -M snowy                   # Specify the machine (optional)
#SBATCH -p core                    # Specify the partition or queue
#SBATCH -n 4                       # Number of cores
#SBATCH -t 2:00:00                # Walltime
#SBATCH -J htseq_count             # Job name
#SBATCH --mail-type=ALL            # Email notification
#SBATCH --mail-user ramo1628@student.uu.se  # Your email address for notifications
#SBATCH --output=/home/ramo1628/genome_project/outputs/counting_reads_Htseq/htseq_count_%j.out  # Output file for SLURM job logs

# Load necessary modules
module load bioinfo-tools
module load htseq/0.12.4
#module load htseq/2.0.2

# Set path to the directory containing sample directories
ALIGNMENT_DIR="/home/ramo1628/genome_project/outputs/RNA_mapping_to_genome/aligned_reads"
SAMPLES="SRR6040092_scaffold_06 SRR6040093_scaffold_06 SRR6040094_scaffold_06 SRR6040096_scaffold_06 SRR6040097_scaffold_06 SRR6156066_scaffold_06 SRR6156067_scaffold_06 SRR6156069_scaffold_06"

# Prepare the list of BAM files
BAM_FILES=""
for SAMPLE_NAME in $SAMPLES; do
    BAM_FILES+=" ${ALIGNMENT_DIR}/${SAMPLE_NAME}/Aligned.sortedByCoord.out.bam"
done

# Set path to the GTF file
GTF_FILE="/home/ramo1628/genome_project/outputs/gene_annotation_braker/annotation_scaffold_06/GeneMark-ET/genemark.gtf"

# Run HTSeq-count with the following options:
# -f bam: Specify the format of the input file as BAM format.
# -r pos: For paired-end reads, use the position of the reads to assign them to features.
# -s no: Disable strandedness information. This is useful when the data does not have strand-specific information.
# -i gene_id: Use the 'gene_id' attribute in the GTF file to define features (e.g., genes).
# -t exon: Use "exon" as the feature type
# -m union: Handle reads overlapping more than one feature by union.
# -a 4: Use 4 cores for parallel execution.
htseq-count -f bam -r pos -s no -i gene_id -t exon -m union -a 4 ${BAM_FILES} "${GTF_FILE}" > /home/ramo1628/genome_project/outputs/counting_reads_Htseq/counts_read.txt
