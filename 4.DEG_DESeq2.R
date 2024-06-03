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
setwd('M:/')
source('./Codes/My_RNA_plot_functions.R')


#Making the dataframe used for interation between conditions
{dat <- crossing(x=1:3,y=1:3)
dat2 <- dat %>% 
  group_by(grp = paste0(pmin(x, y), pmax(x, y))) %>%
  filter(!x==y) %>%
  filter(row_number() == min(row_number())) %>%
  ungroup() %>%
  #filter(!grp%in%c("13","24"))%>%
  select(-grp)
}

#-----------------
#1.Read in Coldata
#-----------------
{ #Analyses in 

  samples <- read.table("Sample.txt", header = TRUE,stringsAsFactors = T)
  samples=samples[!samples$IID==1,] 
  
# Prepare the class  of columns
{sampleTable <- samples
    sampleTable$group <- as.factor(sampleTable$group )
    sampleTable$IID <- as.factor(sampleTable$IID )
    sampleTable$Zgroup <- as.factor(sampleTable$Zgroup )}
head(sampleTable)
row.names(sampleTable) <- sampleTable$ID #make the rownames= ID 
}

#Analyses in container
#-----------------
#2.Read in featurecounts output 
#-----------------

a <- fread("Output_featurecounts.txt",header=T)
meta <- a[,1:6]
cts <- a[,7:ncol(a)]
#rm(a)
head(cts)
{id<- names(cts) #change the ID od columns
  id2 <- str_split(id, "/",simplify = T)[,13]
  id <- str_split(id2,"A",simplify = T)[,1]
  names(cts) <- id}
cts <- as.matrix(cts)
row.names(cts)=meta$Geneid
head(cts)
#Make sure the order in colsata is the same as in count file
cts<- cts[,rownames(sampleTable)]#cts has to be matrix
all(rownames(sampleTable) == colnames(cts))
dim(cts)

#-----------------
#3. Run DESeq2 
#-----------------
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=sampleTable,
                              design =~group) #
#dds is now ready for DESeq() see DESeq2 vignette

#keep <- rowSums(counts(dds)) >= 10  #a minimal pre-filtering to keep only rows that have at least 10 reads total.
# table(keep)
# dds <- dds[keep,]
# dim(dds)

exprSet=assay(dds)
#Summary a table with the number of gene expressed (non zero) in each sample
SAMPLEgene <- as.data.frame(apply(exprSet,2, function(x) sum(x>0)))
 names(SAMPLEgene) <- paste0("expressed_gene_count_out_of_",dim(dds)[1])
 SAMPLEgene$Per <- round(SAMPLEgene[,1]/dim(dds)[1],digits=5)
# write.csv(SAMPLEgene,paste0("Ind_expressed_gene_count.csv"),row.names = T)

dim(dds)
dds <- DESeq(dds,minReplicatesForReplace=Inf)###for schijd

#  ###DEGs between PHENOTYPES

###########Gene of interest for highlight
GeneOfInterest<-c("IGF-1","HMGA2")
GeneOfInterest<-read.table('Gene_of_interest.txt',header=F)$V1


####Interate between conditions
for (i in 1:length(dat2$x)){
  #i=3
  Base=dat2$x[i]
  Case=dat2$y[i]
  print(paste('Base:', dat2$x[i],
    ';Case:', dat2$y[i])) # the bigger number is Base, and small number is Case
 #Case="3"
 #Base="1"
Pair=paste(Case,Base,sep="_")
wd=paste0("M:/ROSLIN/Klebsiella_variicola_RNA/",Pair)
setwd(wd)
dds$condition <- factor(dds$group, levels = c(Case,Base))
dds$condition <- relevel(dds$condition, ref = as.character(Base)) #as.char, otherwise give error in this command
res <- results(dds,contrast=c("group",Case,Base)) #Contrast: The level given last is the base level for the comparison.
resultsNames(dds)
res <- lfcShrink(dds,contrast=c("group",Case,Base), res=res, quiet=T,type="normal")####Shrink log2 fold changes

head(res)

resOrdered <- res[order(res$padj),]#Sort summary list by p-value
head(resOrdered)
DEG =as.data.frame(resOrdered)
#DEG = na.omit(DEG) 


# ####Output 1: all colums
{DEseq_DEG <- DEG
  nrDEG=DEseq_DEG #[,c(2,6)]
  nrDEG$Gene <- rownames(nrDEG)
  #colnames(nrDEG)=c('log2FoldChange','padj')
  #head(nrDEG)
  prefix=paste0("K_",Pair,"_DEG")
  print(prefix)
  write.table(nrDEG,paste0("1DESeq2_output_",prefix,".txt"),row.names = F,quote = F,sep="\t")
}

# #OUtput 2: Plot volcano and UP/Down regulated
source('./Codes/My_RNA_plot_functions.R')
{DEseq_DEG <- DEG
  DEseq_DEG=DEseq_DEG[,c(2,6)] #5: unadjusted P-value 6:Padj
  colnames(DEseq_DEG)=c('log2FoldChange','padj')
  prefix=paste("Conditions",Case,"VS",Base,sep="_")
  draw_v_U_D(exprSet,DEseq_DEG,paste0(prefix))}

# #rm(list=ls())
}#Iterative between combination





