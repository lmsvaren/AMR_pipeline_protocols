#!/bin/bash
#SBATCH --job-name=move1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mem=5GB

function display_help() {

        echo
        echo "Usage: move1.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory that holds files with adapter removed"
	echo "  2: path to directory to move files to"
        echo
	echo "Example: sbatch move1.sh /path/to/files /path/to/move/files/to"
	echo
        echo
        echo

}

if [[ $1 == "--help" ]]; then
        display_help
        exit 0
fi



echo "move files after remove adapters"

# Assign variables and echo for debugging
INPUT_DIR=$1
OUTPUT_DIR=$2

echo $INPUT_DIR
echo $OUTPUT_DIR

cd $INPUT_DIR
mv *remove_adapter_fastq.gz $OUTPUT_DIR
echo "Done"

