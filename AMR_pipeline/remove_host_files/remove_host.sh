#!/bin/bash
#SBATCH --job-name=host
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mem=80GB
#SBATCH --array=0-100 # Edit this if more than 100 samples

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################


function display_help() {

        echo
        echo "Usage: remove_host.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory with original reference .fa and .gtf files"
        echo "  2: .fa reference file name"
        echo "  3: .gtf file name"
	echo "  4: directory to fastq files with adapters removed (*R1_remove_adapter.fastq.gz and *R2_remove_adapter.fastq.gz)"
	echo "  5: directory to reference genome output from step2"
	echo "  6: overhang value selected"
        echo
	echo "Example: sbatch remove_host.sh /path/to/reference/ ref.fa ref.gtf /path/to/fastqs_no_adapters /path/to/overhang/ref n no_host"
	echo
        echo "arguments specified within the STAR command"
        echo "need to be edited manually in this file"
        echo
        echo

}


if [[ $1 == "--help" ]]; then
        display_help
	exit 0
fi

echo "Remove host"

# Load modules
module purge #clear loaded modules

module load tophat/2.0.13
module load bowtie2/2.5.1
module load samtools/1.17
module load bwa/0.7.17
module load star/2.7.10b
module load cufflinks/2.2.1

# Assign variables and echo for debugging
REFERENCE_DIR=$1
FASTA_REF=$2
GTF=$3
FILES_DIR=$4
OVERHANG_DIR=$5
OVERHANG_VALUE=$6

echo $REFERENCE_DIR
echo $FASTA_REF
echo $GTF
echo $FILES_DIR
echo $OVERHANG_DIR
echo $OVERHANG_VALUE

# Get files to be mapped
cd $FILES_DIR

FILES=($(ls -1 *R1_remove_adapter_fastq.gz))
fastq1=${FILES[$SLURM_ARRAY_TASK_ID]}
fastq2=$(echo $fastq1 | sed "s/R1_remove_adapter_fastq.gz/R2_remove_adapter_fastq.gz/")
align_out_name=$(echo $fastq1 | sed "s/R1_remove_adapter_fastq.gz/host_removed_/")

echo $fastq1
echo $fastq2
echo $align_out_name

# Perform mapping with STAR
STAR --genomeDir $OVERHANG_DIR \
	--readFilesCommand zcat \
	--readFilesIn $fastq1 $fastq2 \
	--runThreadN 24 \
	--outFilterScoreMinOverLread 0.3 \
	--outFilterMatchNminOverLread 0.3 \
	--outFilterMismatchNmax 2 \
	--outFileNamePrefix $align_out_name \
	--outReadsUnmapped Fastx \
	--alignIntronMax 1000000 \
	--sjdbGTFfile $REFERENCE_DIR/$GTF \
	--sjdbOverhang $OVERHANG_VALUE \
	--outSAMtype BAM Unsorted SortedByCoordinate \
	--outSAMattrRGline ID:$align_out_name PL:HiSeq SM:$align_out_name LB:RNA \
	--bamRemoveDuplicatesType UniqueIdentical \
	--outFilterIntronMotifs RemoveNoncanonical \
	--twopassMode Basic \
	--quantMode TranscriptomeSAM GeneCounts

