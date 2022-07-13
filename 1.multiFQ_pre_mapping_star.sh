#!/bin/bash
#make Fastq_list_1
#ls */*.fq.gz >Fastq_list
ls */*1.fq.gz>Fastq_list_1
sed -i 's/\//\t/g' Fastq_list_1
sed -i 's/1\.fq\.gz/ /g' Fastq_list_1

#make sample list
for prefix in `less ../Sample_RNA2.txt`; do grep $prefix Fastq_list_1|cut -f 2 >${prefix}_list;done

#make mulFQ file for mapping
for prefix in `less ../Sample_RNA2.txt`;do \
#To make ${prefix}_mulFQ
#READS 1
for i in `less ${prefix}_list`; do \
printf  "${prefix}/${i}1.fq.gz,"
done >>${prefix}_mulFQ
#space
printf " ">> ${prefix}_mulFQ
#READS 2
for j in `less ${prefix}_list`; do \
printf  "${prefix}/${j}2.fq.gz,"
done >>${prefix}_mulFQ
sed -i 's/\, / /' ${prefix}_mulFQ
sed -i 's/\,$/ /' ${prefix}_mulFQ
done

#If error, remember to 
#rm *_mulFQ
