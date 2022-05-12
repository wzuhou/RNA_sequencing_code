#!/bin/sh

# Grid Engine options
#$ -N trimmo
#$ -cwd
#$ -l h_vmem=6G
#$ -pe sharedmem 6
#$ -l h_rt=10:30:00
#$ -M zwu33@ed.ac.uk
#$ -m baes

# If you plan to load any software modules, then you must first initialise the modules framework.
. /etc/profile.d/modules.sh

export OMP_NUM_THREADS=$NSLOTS

module load igmm/apps/trimmomatic/0.39
#module load java

#cp /exports/igmm/software/pkg/el7/apps/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa ./
dir=<file_dir>
input1=<read1.fastq.gz> #$1
input2=<read2.fastq.gz> #$2

echo "inputfile Forward R1: $input1 and Reverse R2 $input2"
#java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog
#<logFile>] >] [-basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> |
#<unpaired output 1> <paired output 2> <unpaired output 2> <step 1>

#Create folders for Trimmed reads
mkdir -p Trim
if [ ! -f TruSeq3-PE.fa ];then
cp /exports/igmm/software/pkg/el7/apps/trimmomatic/0.36/adapters/*.fa ./
fi
 
#Run Trimmomatic
trimmomatic PE -threads 6 -phred33 ${dir}/$input1 ${dir}/$input2 Trim/paired_$input1 Trim/unpaired_$input1 Trim/paired_$input2 Trim/unpaired_$input2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
