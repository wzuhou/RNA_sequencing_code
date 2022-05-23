#!/bin/bash
# Grid Engine options
#$ -N featurec
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=16:00:00
#$ -M zwu33@ed.ac.uk
#$ -m baes

. /etc/profile.d/modules.sh
date
featurecounts=/exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Install/subread-2.0.0-Linux-x86_64/bin/featureCounts
#OR
module load roslin/subread/1.5.2

# One sample
# sample=$1
# $featurecounts -p -C -t exon -g gene_id -T 12 -a star_index/Gallus_gallus.GRCg6a.95.gtf -o F_quantification/${sample} mapping/${sample}.Aligned.sortedByCoord.out.newheader.RG.bam >& F_quantification/$sample.featurecounts_sortedByCoord.log

#In batch
featureCounts -p -C -t exon -g gene_id -T 12 -a star_index/Gallus_gallus.GRCg6a.95.gtf -o /exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Code_shared/Output_samples /exports/cmvm/eddie/eb/groups/smith_grp/Zhou_Wu/mapping/*.Aligned.sortedByCoord.out.bam >& /exports/cmvm/eddie/eb/groups/smith_grp/Zhou_wu/Code_shared/Output_samples.featurecounts_sortedByCoord.log

