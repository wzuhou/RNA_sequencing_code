#!/bin/bash
#make sample list
for prefix in `less LALO_RNA_list.txt`; do grep $prefix Fastq_list_1|cut -f 2 >${prefix}_list;done

#make mulFQ file for mapping
for prefix in `less LALO_RNA_list.txt`;do \
#SubDir=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/RNA_LALO/mhindle/LALO_RNASeq_newgenome/Masked_mapping/Rawdata/$prefix
#To make ${prefix}_mulFQ
#READS 1
for i in `less ${prefix}_list`; do \
printf  "${prefix}/paired_${i}1.fastq.gz,"
done >>${prefix}_mulFQ
#space
printf " ">>${prefix}_mulFQ
#READS 2
for j in `less ${prefix}_list`; do \
printf  "${prefix}/paired_${j}2.fastq.gz,"
done >>${prefix}_mulFQ
sed -i 's/\, / /' ${prefix}_mulFQ
sed -i 's/\,$/ /' ${prefix}_mulFQ
done
