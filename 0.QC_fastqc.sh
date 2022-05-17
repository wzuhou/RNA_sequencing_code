#!/bin/sh

# Grid Engine options
#$ -N fastqc
#$ -cwd
#$ -l h_vmem=4G
#$ -pe sharedmem 8
#$ -M zwu33@ed.ac.uk
#$ -m baes



# If you plan to load any software modules, then you must first initialise the modules framework.
. /etc/profile.d/modules.sh

export OMP_NUM_THREADS=$NSLOTS

# Then, you must load the modules themselves
#module load java
module load roslin/fastqc/0.11.7

# Run FastQC
mkdir -p fastqc
#zcat ./example/*fastq.gz |fastqc stdin --outdir=./toydata/fastqc/
#fastqc ./example/toydata.fastq.gz
fastqc <path_to_file>/*fastq.gz -o fastqc/ -t 8

# Multiqc 
#module load anaconda
#source activate py37

module load roslin/conda/4.9.1
source /exports/applications/apps/community/roslin/conda/4.9.1/etc/profile.d/conda.sh
conda activate multiqc
export OMP_NUM_THREADS=1 ###Use the number of core you have requested
multiqc fastqc/ -o multiqc -n MultiQC_output

