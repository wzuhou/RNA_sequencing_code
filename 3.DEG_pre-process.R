#!/bin/Rscript
################
#RNA-seq Work Flow
#Stage
################
{library(data.table)
library(stringr)
library(DESeq2)
library(dplyr)
library(tidyr)
}
setwd('D:/RNA_Seq_LALO/')
samples <- read.table("D:/RNA_Seq_LALO/Sample_strip_ID_Tissue.txt", header = TRUE,stringsAsFactors = T)
#head(samples)

a <- fread("D:/RNA_Seq_LALO/psiclass08022021_vote_Modified_Gene_name.count",header=T)
meta <- a[,1:6]
meta[1:6,1:6]
cts <- a[,7:ncol(a)]
rm(a)
head(cts)
{id<- names(cts) #change the ID of columns
id <- str_split(id, "/",simplify = T)[,1]
#id <- str_split(id,"[.]",simplify = T)[,1]
names(cts) <- id
}
cts <- as.matrix(cts)
row.names(cts)=meta$Geneid
#head(cts)
