#!/bin/bash
#SBATCH --job-name=getSRA
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --mem=20GB

################################################################################
################################################################################
############################ DO NOT EDIT THIS SCRIPT ###########################
################################################################################
################################################################################

function display_help() {
        echo
        echo "Usage: sratoolkit_download.sh [arguments]"
        echo
        echo "Order of arguments"
        echo "  1: Accession text file (located in current directory)"
	echo "  2: Output directory"
        echo
	echo "Example: sbatch sratoolkit_download.sh output_dir accession.txt"
	echo
        echo "arguments specified within the fastq-dump command"
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
module load sratoolkit/3.1.1


accessions=$1
OUTDIR=$2
sra_dir=$OUTDIR/"SRAs"/
fastq_dir=$OUTDIR/"fastqs"/
sras=${accessions%.txt}_sra.txt

echo $accessions
echo $sras
echo $OUTDIR
echo $sra_dir
echo $fastq_dir


mkdir $sra_dir
mkdir $fastq_dir

cat $accessions | xargs -I {} prefetch -O "$sra_dir" {}

sed "s|^|${sra_dir}|" $accessions > temp.txt
sed 's|\(.*\)/\([^/]*\)|\1/\2/\2.sra|' temp.txt > $sras
rm temp.txt

cat $sras | xargs -I {} -P 8 bash -c 'fastq-dump --split-files --gzip --skip-technical --readids --read-filter pass --dumpbase --outdir '"$fastq_dir"' {}'








