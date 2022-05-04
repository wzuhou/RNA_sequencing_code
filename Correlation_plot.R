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
#pheatmap::pheatmap(cor(exprSet))
group_list
tmp=data.frame(group_list[,c(3,11,14,18)])
#tmp=data.frame(group_list)

ID=data.frame(breed=group_list$ID)
rownames(tmp)=colnames(exprSet)
#colnames(Breed)=group_list$Breed
#rownames(breed)=colnames(exprSet)
# WiDZ5DQy1>5DO`KFPTS&8CJGR*8_SZWi<d5D#!
pheatmap::pheatmap(cor(exprSet,method = "spearman"),
                   treeheight_row=0,
                   #filename = 'cor_allgenes.png',
                   labels_col=ID$breed,
                   #labels_row=breed$breed,
                   annotation_row = tmp)
dev.off()
exprSet=assay(dds)
dim(exprSet)#Dimensions of an Object#expression martix
#exprSet=exprSet[(apply(exprSet,1, function(x) sum(x>1) > 20)),]#apply(X, MARGIN, FUN, !-)#nothing changes, if all ind. are true
keep <- rowSums(counts(dds)) >= 480
table(keep)
dds <- dds[keep,]
exprSet=assay(dds)
dim(exprSet)
exprSet=log2(edgeR::cpm(exprSet)+1)#transform the value>log
dim(exprSet)
#exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:1000]),]#Median Absolute Deviation#top 500?#>x6TVPN;2nJGR;VVM3<F@k2n5D2bA?6xGR#,MADJGR;VVB30tM3<FA?#,1H1jW<2n8|D\JJS&J}>]</VP5DRl3#V5!#6TSZ1jW<2n#,J9SC5DJGJ}>]5=>yV55D>`@kF=7=#,KyRT4s5DF+2nH(VX8|4s#,Rl3#V56T=a9{R2;a2zIzVXR*S0Ol!#6TSZMAD#,IYA?5DRl3#V52;;aS0OlWnVU5D=a9{!#
M=cor(log2(exprSet+1),method = "spearman")#correlation matrix
#M=cor((exprSet+1),method = "spearman")#correlation matrix


tmp=data.frame(group_list[,c(3,11,14,18)])
#tmp=data.frame(group_list)
#rownames(tmp)=colnames(M)

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

#################
#PCA
#################
library("wesanderson")
##Input: prepare the count matix that has been normalized and log transformed
#Vs normalized then pca
vsdata <- varianceStabilizingTransformation(dds,blind = F)#nsub the number of genes to subset (1000)
#plotPCA(vsdata, intgroup="group_list") #using the DESEQ2 plotPCA fxn we can
pcaData <- plotPCA(vsdata, intgroup=c("Tissue"), returnData=TRUE)#intgroup= a character vector of names in colData(x) to use for grouping

# {# also possible to perform custom transformation:
# dds <- estimateSizeFactors(dds)
# # shifted log of normalized counts
# se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
#                            colData=colData(dds))
# # the call to DESeqTransform() is needed to
# # trigger our plotPCA method.
# plotPCA( DESeqTransform( se ),intgroup=c("Tissue") )
# }
percentVar <- round(100 * attr(pcaData, "percentVar"))###PCs explaned variance
head(pcaData)
all(rownames(sampleTable) == rownames(pcaData))
df <- cbind(pcaData,sampleTable)
df <- df[,-3]
head(df)
#df$Stage_Sex<- paste0("E",df$Stage,"_",df$Sex) #Giving another group for plotting
#-----------
# Visualize #PCA
PC1_stage <- ggplot(df, aes(x = Tissue, y = PC1, fill = Tissue)) +
  geom_boxplot(fill = c("#81A88D" ,"#FDD262","#D3DDDC","#C7B19C"),coef=1e30) +#ignore the outlier by changing the threshold
  geom_jitter(shape = 21,
              fill = "gray60",alpha=0.9,
              position = position_jitter(0.21)) +
  theme_classic()
PC1_stage
ggsave("Anova_PC1_stage.jpg",PC1_stage,width = 6,height = 4)
#names(df)[9] <- "Breeds"
##PCA
{pca <- 
    ggplot(df, aes(PC1, PC2,color=Tissue,shape=Conditions,label=ID)) +
    #geom_point(size=10, aes(fill=factor(breed_name), alpha=as.character(Sex))) +
    geom_point(size=7,alpha=0.8)+
    #scale_shape_manual(values=c(21,24)) +
    scale_shape_manual(values=c(15,16,17,18))+
    #scale_alpha_manual(values=c("ZZ"=0, "ZW"=1))+
    scale_color_manual(values=c(wes_palette("Darjeeling1",n=5,type ="discrete"),wes_palette("Darjeeling2",n=3,type ="discrete"))) +
    #scale_colour_manual(values=rep(brewer.pal(8,"Set1"),times=1))+
    #scale_x_continuous(limits = c(-50,50),expand = c(0, 0))+
    scale_y_continuous(limits = c(-50,50),expand = c(0, 0))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()+  
    geom_text(vjust =0 , nudge_y =-4,check_overlap = F,size=3)+#,check_overlap = T
    # geom_text(hjust =0.5, nudge_x =5,check_overlap =T,size=3)+#,check_overlap = T
    #guides(colour = guide_legend(title ="Breed"),keyheight = 0.8,label.vjust=1,ncol = 5)+  
    theme_light()+
    theme(legend.position = "bottom")+
    guides(colour = guide_legend(label.hjust=0.01,label.vjust=0.5,ncol =4,override.aes = list(size =2)),
           shape = guide_legend(label.hjust=0.01,label.vjust=0.5,ncol =2),override.aes = list(size =2.5))
}
pca
ggsave("d:/RNA_Seq_LALO/PCA_gene_expression.JPG",width = 7,height =7,pca)

