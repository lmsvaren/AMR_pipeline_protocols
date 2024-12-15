#!/bin/bash
#SBATCH --job-name=refIndex
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --mem=70GB

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################

function display_help() {

        echo
        echo "Usage: make_host_reference.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory that holds the reference .fa and .gtf files"
        echo "  2: .fa file name"
        echo "  3: .gtf file name"
        echo "  4: path to save output"
	echo "  5...n: overhang value(s). There can be 1 to n values"
	echo
	echo "Example: sbatch make_host_reference.sh /path/to/reference/ ref.fa ref.gtf /path/to/output/ n..."
	echo 
	echo "You can create multiple indexes in one run by listing multiple overhang values"
        echo "arguments specified within the STAR command"
        echo "need to be edited manually in this file"
	echo
	echo

}

if [[ $1 == "--help" ]]; then
        display_help
        exit 0
fi



# Load modules
module purge #clears loaded modules
module load star/2.7.10b

# Assign variables
INPUT_DIR=$1
FASTA=$2
GTF=$3
OUTPUT_DIR=$4

# Output paths for debugging
echo $INPUT_DIR
echo $FASTA
echo $GTF
echo $OUTPUT_DIR

# Skip first 4 variables
shift 4

# Make a reference for every overhang value given
for value in "$@"; do
	
	echo "Overhang value: $value"
	genome_dir=$OUTPUT_DIR/${FASTA%.*}_star_${value}bp

	# Create output dir if not already present
	if [ ! -d "$genome_dir" ]; then
		mkdir $genome_dir
		echo "Directory made: $genome_dir"
	else
		echo "Directory exists: $genome_dir"
	fi
	
	echo "Running"	
	STAR --runMode genomeGenerate \
		--runThreadN 16 \
		--genomeDir $genome_dir \
		--genomeFastaFiles $INPUT_DIR/$FASTA \
		--sjdbGTFfile $INPUT_DIR/$GTF \
		--sjdbOverhang $value
done
