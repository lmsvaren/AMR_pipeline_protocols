#!/bin/bash
#SBATCH --job-name=adapters
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mem=30GB
#SBATCH --array=0-100 # Edit this if more than 100 samples

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################

function display_help() {

        echo
        echo "Usage: sbatch remove_adapters.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory that holds the adapter .fa file"
	echo "  2: path to fastq files (ends with *_1.fastq.gz)"
        echo "  3: path to bbduk.sh script"
        echo
	echo "Example: sbatch remove_adapters.sh /path/to/adapter/file /path/to/fastqs /path/to/script/bbduk.sh"
	echo
        echo "arguments specified within the bbduk.sh script"
        echo "need to be edited manually in this file"
        echo
        echo

}

if [[ $1 == "--help" ]]; then
        display_help
        exit 0
fi



echo "remove adapters"

# Load modules
module purge #clear loaded modules
module load java/17

# Assign variables and echo for debugging
REF_ADAPTERS=$1
FASTQ=$2
BBDUK_SCRIPT=$3

echo $REF_ADAPTERS
echo $FASTQ

# Select pass one files
cd $FASTQ
FILES=($(ls -1 *_1.fastq.gz))
fastq1=${FILES[$SLURM_ARRAY_TASK_ID]}
fastq2=$(echo $fastq1 | sed 's/_1.fastq.gz/_2.fastq.gz/')

# Set names for files with adapters removed
bbduk_R1=$(echo $fastq1 | sed 's/_1.fastq.gz/_R1_remove_adapter_fastq.gz/')
bbduk_R2=$(echo $fastq2 | sed 's/_2.fastq.gz/_R2_remove_adapter_fastq.gz/')

# Print for debugging
echo $fastq1
echo $fastq2
echo $bbduk_R1
echo $bbduk_R2

# Run bbduk.sh script to remove adapters script
sh $BBDUK_SCRIPT -Xmx20g \
	in1=$fastq1 in2=$fastq2 \
	out1=$bbduk_R1 out2=$bbduk_R2 \
	ref=$REF_ADAPTERS \
	ktrim=r k=23 mink=11 hdist=1 tpe tbo minlength=50


echo "Done"

