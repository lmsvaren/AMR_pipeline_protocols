#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mem=50GB
#SBATCH --array=0-100 # Edit this if more than 100 samples

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################

function display_help() {
	echo
	echo "Usage: salmon_quantification.sh [arguments]"
	echo
	echo "Order of arguments"
	echo "	1: directory path of the files to quantify"
	echo "	2: directory path where index file is"
	echo
	echo "Example: sbatch salmon_quantification.sh /path/to/fastqs /path/to/index"
	echo
	echo "arguments specified within the salmon command"
	echo "need to be edited manually in this file"
	echo
	echo
}

if [[ $1 == "--help" ]]; then
	display_help
	exit 0
fi

echo "AMR quantificaiton"

# Load modules
module purge #clear loaded modules
module load salmon/1.10.0

# Define variables
DIR=$1
TRANSCRIPT_INDEX=$2

echo $DIR
echo $TRANSCRIPT_INDEX

# Go to file directory
cd $DIR

# Select mate1.gz files
FILES=($(ls -1 *mate1.gz))

# Set fastq1 and fast2
fastq1=${FILES[$SLURM_ARRAY_TASK_ID]}
fastq2=$(echo $fastq1 | sed "s/mate1.gz/mate2.gz/")


# Set the name of the output folder
out_result=$(echo $fastq1 | sed "s/pass_host_removed_Unmapped.out.mate1.gz/Salmon_quant/")

# Run salmon
salmon quant -i $TRANSCRIPT_INDEX \
	-l A \
	-1 $fastq1 \
	-2 $fastq2 \
	--validateMappings -o $out_result


