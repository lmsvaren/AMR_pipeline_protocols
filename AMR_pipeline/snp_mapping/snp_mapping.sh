#!/bin/bash
#SBATCH --job-name=SNP # Job name CHANGE
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=2                   # Run a single task
#SBATCH --cpus-per-task=8            # Number of CPU cores per 
#SBATCH --partition=medium # medium partition has 7 day limit
#SBATCH --mem=100gb                    # Total memory limit
#SBATCH --array=0-100 #CHANGE

## Load need modules
module purge
module load samtools/1.17
module load bcftools/1.16
module load bwa/0.7.17
module load bowtie2/2.5.1


# Path to reference index
# Change only the path up to AMR_pipeline
genome_fasta=/path/to/AMR_pipeline/snp_mapping/megares_database/megares_database_v3.00_index.fasta

# Path to output directory to save the VCFs in
output_dir=/path/to/vcf/output/
input_dir=/path/to/remove_host/

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################

# Move to the input directory
cd $input_dir

echo "genome_fasta", $genome_fasta

# Declare file names
FILES=($(ls -1 *mate1.gz))
mate1=${FILES[$SLURM_ARRAY_TASK_ID]}
mate2=$(echo $mate1 | sed 's/mate1.gz/mate2.gz/')
#mate2=$(echo $mate1 | sed 's/out1.fastq.gz/out2.fastq.gz')
echo "mate1", $mate1
echo "mate2", $mate2

# Delare main_base for 1 and 2
temp1=${mate1/_pass_host_removed_Unmapped.out.}

base=${temp1/mate1.gz} #remove _mate1.gz
#base=${temp1/out1.fastq.gz}
echo "temp1", $temp1
echo "base", $base

# Declare sam file name
sam_file="${base}_aligned.sam"
echo "sam_file", $sam_file

# Declare bam file names
bamfile="${base}_aligned.bam"
bam_unique_file="${bamfile%.bam}_unique10q.bam"
bam_unique_sorted_file="${bam_unique_file%.bam}_sorted.bam"
echo "bamfile", $bamfile
echo "bam_unique_file", $bam_unique_file
echo "bam_unique_sorted_file", $bam_unique_sorted_file

# Declare bcf file name
bcf_file="${bamfile%.bam}_unique10q_sorted.bcf"
echo "bcf_file", $bcf_file

# Declare vcf file name
variant_file="${base}_variant.vcf"
echo "variant file", $variant_file

###########################################################
# Generate sam file
bwa mem $genome_fasta $mate1 $mate2 -t 6 > $sam_file

#Generate bam file
samtools view -bT $genome_fasta $sam_file > $bamfile

# Filter bam file
samtools view -q 10 -b $bamfile > $bam_unique_file

# Sort filtered bam file
samtools sort $bam_unique_file -o $bam_unique_sorted_file

# Index sorted bam file
samtools index $bam_unique_sorted_file

# Generate bcf
bcftools mpileup -Ou -f $genome_fasta $bam_unique_sorted_file | bcftools call -mv -Ob -o $bcf_file

# Generate vcf
bcftools view $bcf_file > $variant_file

#bai file
bai="${bam_unique_sorted_file}.bai"
echo "bai", $bai

# Remove all files but vcf
rm $sam_file $bamfile $bam_unique_file $bam_unique_sorted_file $bcf_file $bai

# Move vcf file to designate directory
mv $variant_file $output_dir



