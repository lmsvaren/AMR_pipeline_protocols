#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mem=10GB

function display_help() {
	echo
	echo "Usage: move3.sh [arguments]"
	echo
	echo "Order of arguments"
	echo "	1: directory path of the salmon output directories"
	echo "	2: directory path to move the output to"
	echo
	echo "Example: sbatch move3.sh /path/to/results /path/to/move/to"
	echo
	echo
	echo
}

if [[ $1 == "--help" ]]; then
	display_help
	exit 0
fi

echo "Move quant files"

INPUT_DIR=$1
OUTPUT_DIR=$2

echo $INPUT_DIR
echo $OUTPUT_DIR

# Go to file directory
cd $INPUT_DIR

mv *Salmon_quant $OUTPUT_DIR


