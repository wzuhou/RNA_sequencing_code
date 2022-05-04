#Correlation plot

library(DESeq2)
library(edgeR, lib.loc="~/R/win-library/3.6")
library(limma)
library("RColorBrewer", lib.loc="~/R/win-library/3.6")
library("ggplot2", lib.loc="~/R/win-library/3.6")

exprSet=assay(dds)
#group_list=dds$pheno
group_list=colData(dds)
#rm(list = ls())

#Correlation plot
colnames(exprSet)
group_list
tmp=data.frame(group_list[,c(3,11,14,18)])

ID=data.frame(breed=group_list$ID)
rownames(tmp)=colnames(exprSet)


dim(exprSet)#Dimensions of an Object#expression martix
keep <- rowSums(counts(dds)) >= 480
table(keep)
dds <- dds[keep,]
exprSet=assay(dds)
dim(exprSet)
exprSet=log2(edgeR::cpm(exprSet)+1)#transform the value>log
dim(exprSet)
M=cor(log2(exprSet+1),method = "spearman")#correlation matrix

tmp=data.frame(group_list[,c(3,11,14,18)])

Num=brewer.pal(n = 12, name = "Paired")
names(Num)<-unique(samples$Num)
ann_colors = list(
  Tissue = c("Gonad"="#FF0000" , "Heart"="#00A08A","Hypothalamus"="#F98400", "Liver"="#5BBCD6"),
  LH_Sub_Stage = c(Arrival = "#0072B2",  Incubation = "#F0E442" ),
  storm=c(Extreme_Spring="#009E73",No_Storm="#CC79A7",Storm="#D3DDDC"),
  Conditions=c('1'="#81A88D",'2'="#FDD262",'3'="#D3DDDC",'4'="#C7B19C"))

pheatmap::pheatmap(M,
                   color = wesanderson::wes_palette("Zissou1", 50,type="continuous"),
                   annotation_colors = ann_colors,
                   annotation_row  = tmp,
                   #gaps_row = 15,
                   #gaps_col = 15,
                   cutree_rows = 4,
                   cutree_cols = 4,
                   #fontsize_row = fontsize,
                   fontsize_col= 9,
                   width=9,
                   filename = paste0("D:/RNA_Seq_LALO/correlation_allSamples_test.pdf"),
                   labels_col=tmp$Tissue,
                   treeheight_row=0,
                   cluster_rows = T)#plot with correlation matrix
dev.off()


choose_gene=head(rownames(need_DEG),50) ## 50 maybe better
# choose_matrix=exprSet[choose_gene,]
# choose_matrix=t(scale(t(choose_matrix)))

