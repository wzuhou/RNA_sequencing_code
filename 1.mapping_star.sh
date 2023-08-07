#!/bin/sh

#$ -N MappingSTAR
#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 10


. /etc/profile.d/modules.sh
module load roslin/star/2.5.3a
export OMP_NUM_THREADS=$NSLOTS

#################################################################
#Prepare the STAR index first: GRCg6a
#Usage sh 1.mapping_star.sh

#Input1
ref='<path_to_file>/GCF_ABC_genomic.fna'
anno='<path_to_file>/GCF_ABC_genomic.gtf'

#1. basic options to generate genome indices are as follows:
STAR \
--runThreadN $NSLOTS \
--runMode genomeGenerate \
--genomeDir ./GRCg6a \
--genomeSAindexNbases 10 \
--genomeFastaFiles ${ref} \
--sjdbGTFfile ${anno} \
--sjdbOverhang 149

##################################################################
#Mapping usage:
#For loop
#for prefix in `cat name.list25`;do sh 1.mapping_star.sh $prefix; done

#Input2
prefix=$1
listR1=<path to Fastq>/${prefix}/${prefix}_R1.fastq.gz 
listR2=<path to Fastq>/${prefix}/${prefix}_R2.fastq.gz
output=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/mapping/${prefix}

ls  ${listR1} ${listR2}
echo ${output}
mkdir -p /exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/mapping

#2. Mapping reads to reference 
STAR --genomeDir GRCg6a \
--runThreadN $NSLOTS \
--readFilesCommand zcat \
--readFilesIn ${listR1} ${listR2} \
--outFilterType BySJout \
--outSAMunmapped None \
--outFileNamePrefix ${output} \
--outSAMtype BAM SortedByCoordinate
