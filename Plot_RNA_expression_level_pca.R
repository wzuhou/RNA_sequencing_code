{library(dplyr)
  library(ggplot2)
  library("wesanderson")
  library(DESeq2)
  library(edgeR)
  library(limma)
  library("RColorBrewer")
}


GeneOfInterest<-c("FKBP5","COL","PER1","PER2","PER3","SOUL","TBC1D2","CRY1","FMN2","ZP3","TBC1D8",'AMBP','PCK1','UCN3','FBXO32','CPS1','NR1D1','NPC1')

GeneOfInterest_scaffold<-c("FKBP5","PER2","PER3","SOUL","TBC1D2","CRY1","FMN2.1","ZP3","TBC1D8.1",'AMBP','PCK1','UCN3','FBXO32','CPS1','NR1D1','NPC1.1') #remove two gene that do not exist
setwd("M:/ROSLIN/RNA_Seq_LALO/GENE_expression_plot/")
#Visualize 
#########################
#log2 expression level
#########################
{GENE="FKBP5"
gene_name="FKBP5"
TS="Hypothalamus"
}
for (i in GeneOfInterest_scaffold) {
  GENE=i
  gene_name=i
  TS="Gonad"
  # TS="Heart"
  #TS="Hypothalamus"
  # TS="Liver"

{
  try(expr=GENE_all <- plotCounts(dds, gene=GENE, intgroup=c("ID","Tissue","Conditions"),cex=1.2,normalized=F,returnData=T),silent =F)
#GENE_all <- cbind(GENE_all,str_split(row.names(GENE_all),"_",simplify = T)[,1:3])
head(GENE_all)
#names(GENE_all)[3:5] <- c("breed_name","Stage","Numbering")
GENE_TS<- GENE_all[GENE_all$Tissue==paste0(TS),]
# GENE_breed<- GENE_all[GENE_all$breed_name==BD |GENE_all$breed_name==paste0(BD,"B"),]
# GENE_E5<- GENE_E5[GENE_E5$breed_name==BD |GENE_E5$breed_name==paste0(BD,"B"),]
# head(GENE_breed)
}

# # ##ALL Tissues:
# GENE_all$count <- log2(GENE_all$count)#Log2 transform?#GENE_E5| GENE_ALL |tpm_gene_E5
# {w1 <- ggplot(GENE_all,aes(x=Conditions,y=count,fill=Tissue,shape=Tissue))+
#     geom_boxplot(alpha=0.7,color="black")+#,position = "jitter"
#     #geom_point(size=2,alpha=0.7,color="black",position = position_jitter(width =0.1,height = 0))+#
#     scale_fill_manual(values =c("#00AFBB", "#E7B800", "#FC4E07",'darkblue'))+
#     labs(title = paste0(gene_name),subtitle =paste0("All conditions"))+# subtitle = paste0("Stage E",S,"_",BD,"_",gene_name)+
#     #geom_smooth(count~ Conditions, span = 0.8)+
#     facet_wrap(~Tissue,drop=F)+
#     theme_bw()+
#     xlab(NULL) +
#     scale_y_continuous(expand = c(0,1))+
#     scale_shape_manual(values=c(24,23,22,21))+#E5: 24,23, E13: 22,21
#     #scale_x_discrete(breaks=NULL)+
#     ylab("log2 expression count" ) + #"log2 expression count" "TPM"
#     theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
#           plot.subtitle = element_text(hjust = 0.5,size = 15,face = "bold"),
#           axis.text=element_text(size=15),
#           axis.title=element_text(size=15),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank())
#   w1}
#  ggsave(paste0("M:/ROSLIN/RNA_Seq_LALO/GENE_expression_plot/Gene","_",gene_name,".jpg"),w1,width = 9,height = 3.5)
#  ggsave(paste0("M:/ROSLIN/RNA_Seq_LALO/GENE_expression_plot/Gene","_",gene_name,".PDF"),w1,width = 9,height = 3.5)
# #Stage_specific Log2 transform?
###GENE_Tissue
GENE_TS$count <- log2(GENE_TS$count) 
{w2 <- ggplot(GENE_TS,aes(x=Conditions,y=count,colour=Conditions))+
  geom_point(size=7,alpha=0.7)+
  scale_colour_manual(values =c("#46ACC8","#B40F20","#EBCC2A" ,"grey50"),guide="none")+ 
    #scale_shape_manual(values=c(15,19))+#E5: 17,18, E13: 15,19
  labs(title = paste0(gene_name),subtitle =paste("LALO",TS,sep='_'))+# subtitle = paste0("Stage E",S,"_",BD,"_",gene_name)
  theme_bw()+
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0.5))+
  #scale_x_discrete(breaks=NULL)+
  ylab("log2 expression count") +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 15,face = "bold"),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
w2}
ggsave(paste0("M:/ROSLIN//RNA_Seq_LALO/GENE_expression_plot/Gene_",TS,"_",gene_name,".jpg"),w2,width = 3.5,height = 2.5)
ggsave(paste0("M:/ROSLIN/RNA_Seq_LALO/GENE_expression_plot/Gene_",TS,"_",gene_name,".PDF"),w2,width = 3.5,height = 2.5)

}



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
ggsave("Anova_PC1_stage.pdf",PC1_stage,width = 6,height = 4)
names(df)[9] <- "Breeds"
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
  #scale_x_continuous(limits = c(-80,80),expand = c(0, 0))+
  #scale_y_continuous(limits = c(-50,50),expand = c(0, 0))+
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
