#!/bin/sh

#$ -N MappingSTAR
#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 10
#$ -P roslin_smith_grp

. /etc/profile.d/modules.sh
#To run this code
#for prefix in `cat name.list25`;do sh ZW.examp1.all25tot2.sh $prefix; done

module load roslin/star/2.5.3a
prefix=$1

listR1=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Fastq/${prefix}/${prefix}_R1.fastq.gz 
listR2=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Fastq/${prefix}/${prefix}_R2.fastq.gz
output=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/mapping/${prefix}

ls  ${listR1} ${listR2}
echo ${output}
#Prepare the STAR index first: GRCg6a

STAR --genomeDir GRCg6a \
--runThreadN 10 \
--readFilesCommand zcat \
--readFilesIn ${listR1} ${listR2} \
--outFilterType BySJout \
--outSAMunmapped None \
--outFileNamePrefix ${output} \
--outSAMtype BAM SortedByCoordinate
