#!/bin/bash
#SBATCH --job-name=run_all
#SBATCH --nodes=1
#SBATCH --partition=medium # Change as needed according to your HPC
#SBATCH --mem=100GB        # Change as needed

#############################################################################################################
# Move to the location of this file
CWD=$(pwd)

# Set the output directory (step 7)
OUTPUT_DIR=/path/to/output/

# Remove adapter variables (step 8)
ORIGINAL_FASTQS=/path/to/fastqs/

# Remove host variables (step 9)
REFERENCE_DIR=/path/to/reference/directory/with/fa_and_gtf/files/
FASTA_REF_FILE=reference_genomic.fa
GTF_REF_FILE=reference_genomic.gtf

# Set overhang directory and value (step 10)
OVERHANG_DIR=/path/to/reference/with/overhang/values/
OVERHANG_VALUE=149 # Change this number to desired overhang value

##############################################################################################################
##############################################################################################################
###################################### DO NOT EDIT ANYTHING BELOW ############################################
##############################################################################################################
##############################################################################################################
remove_adapter_output=$OUTPUT_DIR/remove_adapters
remove_host_output=$OUTPUT_DIR/remove_host
salmon_quant_output=$OUTPUT_DIR/salmon_quant

mkdir $remove_adapter_output
mkdir $remove_host_output
mkdir $salmon_quant_output

echo "removed adapters output directory:" $remove_adapter_output
echo "removed host output directory: " $remove_host_output
echo "salmon output directory: " $salmon_quant_output
echo " "


# For debugging
echo "Current working directory: " $CWD
echo "Output directory: " $OUTPUT_DIR
echo "Original fastq files: " $ORIGINAL_FASTQS
echo "Reference directory: " $REFERENCE_DIR
echo "Fasta file name: " $FASTA_REF_FILE
echo "GTF file name: " $GTF_REF_FILE
echo "Overhang directory: " $OVERHANG_DIR
echo "Overhang value: " $OVERHANG_VALUE
echo " "

##### REMOVE ADAPTERS #####
echo "removing adapters"

# bbduk.sh and adapters.fa are in directory remove_adapter_files
remove_adapter_path=$CWD/remove_adapter_files
bbduk=$remove_adapter_path/bbmap/bbduk.sh
adapters=$remove_adapter_path/adapters.fa

job1=$(sbatch $remove_adapter_path/remove_adapters.sh $adapters $ORIGINAL_FASTQS $bbduk)
job1_id=$(echo $job1 | awk '{print $4}')
echo $job1_id

echo "adapters removed"
#######################
echo "move1"
job2=$(sbatch --dependency=afterany:$job1_id $remove_adapter_path/move1.sh $ORIGINAL_FASTQS $remove_adapter_output)
job2_id=$(echo $job2 | awk '{print $4}')
echo $job2_id

##### REMOVE HOST #####
echo "remove host"

remove_host_path=$CWD/remove_host_files
fasta_no_adapter=$remove_adapter_output

job3=$(sbatch --dependency=afterany:$job2_id $remove_host_path/remove_host.sh $REFERENCE_DIR $FASTA_REF_FILE $GTF_REF_FILE $fasta_no_adapter $OVERHANG_DIR $OVERHANG_VALUE)
job3_id=$(echo $job3 | awk '{print $4}')
echo $job3_id

echo "host removed"
#######################
echo "gzip and move"
job4=$(sbatch --dependency=afterany:$job3_id $remove_host_path/move2.sh $remove_adapter_output $remove_host_output)
job4_id=$(echo $job4 | awk '{print $4}')
echo $job4_id

##### QUANTIFICATION #####
echo "quantifying genes"
salmon_path=$CWD/salmon_quantification
index_dir=$salmon_path/AMR_salmon_index

job5=$(sbatch --dependency=afterany:$job4_id $salmon_path/salmon.sh $remove_host_output $index_dir)
job5_id=$(echo $job5 | awk '{print $4}')
echo $job5_id
echo "transcripts quantified"

########################
echo "move quant folders"
$(sbatch --dependency=afterany:$job5_id $salmon_path/move3.sh $remove_host_output $salmon_quant_output)

echo "All Done"


