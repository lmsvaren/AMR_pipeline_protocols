#!/bin/bash
#SBATCH --job-name=gzip_move
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --mem=10GB

function display_help() {

        echo
        echo "Usage: move2.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory with files, host removed"
        echo "  2: path to move files to"
        echo
	echo "Example: sbatch move2.sh /path/to/files /path/to/move/to"
	echo
        echo
        echo

}


if [[ $1 == "--help" ]]; then
        display_help
	exit 0
fi

echo "gzip and move mate files"

INPUT_DIR=$1
OUTPUT_DIR=$2

cd $INPUT_DIR
ls *mate* | xargs -n 1 -P 4 gzip
mv *mate* $OUTPUT_DIR
