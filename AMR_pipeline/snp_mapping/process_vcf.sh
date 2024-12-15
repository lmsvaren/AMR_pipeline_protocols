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
        echo "Usage: sh process_vcf.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: path to directory with .vcf files"
        echo
        echo "Example: sh process_vcf.sh /path/to/vcf/files/"
        echo "Run this BEFORE running merge_vcf.sh"
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

# Move to the input directory
cd $vcf_dir


# filter VCF
for file in *.vcf; do
        out_file=$(echo "$file" | sed 's/\.vcf$/_filter.vcf/')
        bcftools filter -i 'QUAL > 30 && DP > 10 && AF > 0.01' $file -o $out_file
done

# bgzip the filtered vcf
for file in *_filter.vcf; do bgzip "$file"; done

# Index the filtered vcf files
for file in *.vcf.gz; do bcftools index "$file"; done


