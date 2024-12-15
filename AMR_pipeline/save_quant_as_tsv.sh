#!/bin/bash
#SBATCH --job-name=sf_to_tsv         # Job name CHANGE
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=2                   # Run a single task
#SBATCH --partition=medium    # medium partition has 7 day limit
#SBATCH --mem=10gb                    # Total memory limit

## Load need modules
module purge

# Path to working directory
cd $(pwd)

# Directories with input and output
SALMON_QUANT_DIR=/path/to/salmon/output/directories/
OUTPUT_CSV_FILES=/path/to/save/tsv/files/

# Run the R script to save the quant.sf files
for dir in $SALMON_QUANT_DIR/*; do
	cd $dir
	directory=$(basename "$dir")
	save="${directory%_Salmon_quant}.tsv"
	mv quant.sf $save
	mv *tsv $OUTPUT_CSV_FILES
	echo $save
	cd ../
done



