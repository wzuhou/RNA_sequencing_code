{library(dplyr)
library(ggplot2)
library("wesanderson")}
#Visualize 
#########################
#log2 expression level
#########################
{GENE="PER2"
gene_name="PER2"
TS="Hypothalamus"
#BD="DPBd"
}


GENE_all <- plotCounts(dds, gene=GENE, intgroup=c("ID","Tissue","Conditions"),cex=1.2,normalized=F,returnData=T)
#GENE_all <- cbind(GENE_all,str_split(row.names(GENE_all),"_",simplify = T)[,1:3])
head(GENE_all)
#names(GENE_all)[3:5] <- c("breed_name","Stage","Numbering")
GENE_TS<- GENE_all[GENE_all$Tissue==paste0(TS),]



# ##ALL:
GENE_all$count <- log2(GENE_all$count)#Log2 transform?#GENE_E5| GENE_ALL |tpm_gene_E5
{w1 <- ggplot(GENE_all,aes(x=Tissue,y=count,fill=Tissue,shape=Conditions))+
    geom_point(size=7,alpha=0.4,color="black")+#,position = "jitter"
    scale_fill_manual(values =c(wes_palette("FantasticFox1",n=8,type="continuous")),guide="none")+
    labs(title = paste0(gene_name),subtitle =paste0(TS))+# subtitle = paste0("Stage E",S,"_",BD,"_",gene_name)
    theme_bw()+
    xlab(NULL) +
    scale_y_continuous(expand = c(0,1))+
    scale_shape_manual(values=c(24,23,22,21))+#E5: 24,23, E13: 22,21
    #scale_x_discrete(breaks=NULL)+
    ylab("log2 expression count" ) + #"log2 expression count" "TPM"
    theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5,size = 15,face = "bold"),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  w1}
 ggsave(paste0("D:/RNA_Seq_LALO/GENE_expression_plot/Gene","_",gene_name,".jpg"),w1,width = 9,height = 3.5)
 ggsave(paste0("D:/RNA_Seq_LALO/GENE_expression_plot/Gene","_",gene_name,".PDF"),w1,width = 9,height = 3.5)
#Stage_specific Log2 transform?
###GENE_TS
GENE_TS$count <- log2(GENE_TS$count) 
w2 <- ggplot(GENE_TS,aes(x=Conditions,y=count,colour=Conditions))+
  geom_point(size=7,alpha=0.7)+
  scale_colour_manual(values =c("#46ACC8","#B40F20","#EBCC2A" ,"#E1AF00"),guide="none")+ 
    #scale_shape_manual(values=c(15,19))+#E5: 17,18, E13: 15,19
  labs(title = paste0(gene_name),subtitle =paste0("E",S,"_",BD))+# subtitle = paste0("Stage E",S,"_",BD,"_",gene_name)
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
w2

ggsave(paste0("D:/RNA_Seq_LALO/GENE_expression_plot/Gene",TS,"_",gene_name,".jpg"),w2,width = 3.5,height = 3.5)

#################
#PCA
#################
library("wesanderson")
##Input: prepare the count matix that has been normalized and log transformed
#Vs normalized then pca
vsdata <- varianceStabilizingTransformation(dds,blind = F)#nsub the number of genes to subset (1000)
#plotPCA(vsdata, intgroup="group_list") #using the DESEQ2 plotPCA fxn we can
pcaData <- plotPCA(vsdata, intgroup=c("Tissue"), returnData=TRUE)#intgroup= a character vector of names in colData(x) to use for grouping

percentVar <- round(100 * attr(pcaData, "percentVar"))###PCs explaned variance
head(pcaData)
all(rownames(sampleTable) == rownames(pcaData))
df <- cbind(pcaData,sampleTable)
df <- df[,-3]
head(df)
#df$Stage_Sex<- paste0("E",df$Stage,"_",df$Sex) #Giving another group for plotting
#-----------
# Visualize PC1
PC1_stage <- ggplot(df, aes(x = Tissue, y = PC1, fill = Tissue)) +
  geom_boxplot(fill = c("#81A88D" ,"#FDD262","#D3DDDC","#C7B19C"),coef=1e30) +#ignore the outlier by changing the threshold
  geom_jitter(shape = 21,
              fill = "gray60",alpha=0.9,
              position = position_jitter(0.21)) +
  theme_classic()
PC1_stage
#ggsave("Anova_PC1_stage.pdf",PC1_stage,width = 6,height = 4)
names(df)[9] <- "Breeds"

##PCA
{pca <- 
  ggplot(df, aes(PC1, PC2,color=Tissue,shape=Conditions,label=ID)) +
  geom_point(size=7,alpha=0.8)+
  scale_shape_manual(values=c(15,16,17,18))+
  scale_color_manual(values=c(wes_palette("Darjeeling1",n=5,type ="discrete"),wes_palette("Darjeeling2",n=3,type ="discrete"))) +
  #scale_colour_manual(values=rep(brewer.pal(8,"Set1"),times=1))+
  #scale_x_continuous(limits = c(-80,80),expand = c(0, 0))+
  #scale_y_continuous(limits = c(-50,50),expand = c(0, 0))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+  
  geom_text(vjust =0 , nudge_y =-4,check_overlap = F,size=3)+#,check_overlap = T
  theme_light()+
  theme(legend.position = "bottom")+
  guides(colour = guide_legend(label.hjust=0.01,label.vjust=0.5,ncol =4,override.aes = list(size =2)),
         shape = guide_legend(label.hjust=0.01,label.vjust=0.5,ncol =2),override.aes = list(size =2.5))
}
pca
ggsave("D:/RNA_Seq_LALO/PCA_gene_expression.JPG",width = 7,height =7,pca)
