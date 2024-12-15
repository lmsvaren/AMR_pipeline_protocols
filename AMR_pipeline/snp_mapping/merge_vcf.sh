#!/bin/bash
#SBATCH --job-name=mergeVCF # Job name CHANGE
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=2                   # Run a single task
#SBATCH --cpus-per-task=8            # Number of CPU cores per 
#SBATCH --partition=medium # medium partition has 7 day limit
#SBATCH --mem=100gb                    # Total memory limit

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################

function display_help() {

        echo
        echo "Usage: sh merge_vcf.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory with .vcf files"
	echo "  2: group1 files list (.txt file)"
	echo "  3: group2 files list (.txt file)"
	echo "  4: output file name for group1 (end with .vcf.gz)"
        echo "  5: output file name for group2 (end with .vcf.gz)"
        echo "  6: output directory name for bcftools isec results"
        echo
        echo "Example: sbatch remove_host.sh /path/to/vcf/files/ group1.txt group2.txt output1.vcf.gz output2.vcf.gz isec_output"
        echo
	echo

}

if [[ $1 == "--help" ]]; then
        display_help
        exit 0
fi

## Load need modules
module purge
module load samtools/1.17
module load bcftools/1.16
module load bwa/0.7.17
module load bowtie2/2.5.1


# Path to output directory to save the VCFs in
vcf_dir=$1
list_1=$2
list_2=$3
output_file_1=$4
output_file_2=$5
output_dir=$6

# Move to the input directory
cd $vcf_dir


# Merge all the indexed vcf files together
bcftools merge -l $list_1 -o $output_file_1 -Oz --threads 6
bcftools merge -l $list_2 -o $output_file_2 -Oz --threads 6

# Index the merged files
bcftools index $output_file_1
bcftools index $output_file_2

# Find shared variants
bcftools isec $output_file_1 $output_file_2 -p $output_dir

