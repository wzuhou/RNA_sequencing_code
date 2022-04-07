
#!/bin/sh
###########################################
#                                         #
# Submit a job which uses some installed  #
# applications, using:                    #
# "module load <application>"             #
#                                         #
###########################################

# Grid Engine options
#$ -N map
#$ -cwd
#$ -l h_vmem=8G
#$ -M zwu33@ed.ac.uk
#$ -m baes
#$ -P roslin_smith_grp
#$ -pe sharedmem 4

# If you plan to load any software modules, then you must first initialise the modules framework.
. /etc/profile.d/modules.sh
module load R
# Choose the staging environment
export OMP_NUM_THREADS=$NSLOTS

# Then, you must load the modules themselves
module load igmm/apps/STAR/2.7.8a

#STAR Index should be done before this step

#STAR mapping
mkdir -p mapping
for prefix in `less LALO_RNA_list.txt`;do \

#prefix=$1 #Sample_list.txt e.g. K1_1 $i

#listR1=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Pipeline_RNA/data/K_data/raw_data/${prefix}/${prefix}_1.fq.gz
#listR2=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Pipeline_RNA/data/K_data/raw_data/${prefix}/${prefix}_2.fq.gz
output=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/RNA_LALO/mhindle/LALO_RNASeq_newgenome/Masked_mapping/mapping/${prefix}
ref=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/RNA_LALO/mhindle/LALO_RNASeq_newgenome/Masked_mapping/Mask_ref

#ls  ${listR1} ${listR2}
echo "STAR output ${output}"

# ${prefix}_mulFQ` should be created before this script
mytext=`cat ${prefix}_mulFQ`
echo "STAR --genomeDir ${ref} --runThreadN 4 --readFilesCommand zcat --readFilesIn \
$mytext \
--outFilterType BySJout --outSAMunmapped None --outReadsUnmapped Fastx --outFileNamePrefix ${output} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 7900000000"

done
