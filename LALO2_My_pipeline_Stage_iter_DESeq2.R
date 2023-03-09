################
#RNA-seq Work Flow
#Stage
################
{  library(data.table)
  library(stringr)
  suppressMessages(library(DESeq2))
  library(dplyr)
  library(tidyr)
}

TS="ADRE"
TS="PIT"
# TS="FAT"

TISSUE=c("ADRE","PIT","FAT")#

# #Making the dataframe used for interation between conditions
# dat <- crossing(x=1:4,y=1:4) #Case=y
# dat2 <- dat %>% 
#   filter(!row_number()==2)%>%
#   group_by(grp = paste0(pmin(x, y), pmax(x, y))) %>%
#   filter(!x==y) %>%
#   filter(row_number() == min(row_number())) %>%
#   ungroup() %>%
#   filter(!grp%in%c("13","24"))%>%
#   select(-grp) 
# 
# dat2
#----------------------------------------------
# ifelse: within the folder of tissue, by condition
#----------------------------------------------
source('M:/ROSLIN/RNA_2LALO/Codes/My_RNA_plot_functions.R')
#source('M:/ROSLIN/RNA_2LALO/Codes/My_RNA_plot_functions_simple.R')

{ 
  if (exists('TISSUE')) {
    ITERA=3
  } else{
    ITERA=1  
  }
  
  for (j in seq(ITERA)){  #seq(ITERA)
    #j=4
    try(expr=eval(parse(text=paste0("TS <- TISSUE[",j,"]"))),silent = T)
    try(expr=print(TS),silent = T)
    #-----------------
    #1.Read in Coldata
    #-----------------
    { #Analyses in container
      samples <- read.table("M:/ROSLIN/RNA_2LALO/Sample_LALO2_ID_Tissue.txt", header = TRUE,stringsAsFactors = T)
      #Filtering for Tissue
      if (exists("TS")){
        samples=samples[samples$Tissue==TS,] #select Tissue
        wd=paste0("M:/ROSLIN/RNA_2LALO/",TS)
        setwd(wd)
        samples=subset(samples,Tissue==paste0(TS)) #select Tissue
      } else {
        wd=paste0("M:/ROSLIN/RNA_2LALO/")
        setwd(wd)
      }
      #samples=samples[samples$ID!="193_Heart",]
      ########Filter for specific tissues: heart and hyphothalamus
      #samples=subset(samples,Tissue=="Heart"|Tissue=="Hypothalamus")
      
      # Prepare the class  of columns
      {
        sampleTable <- samples
        sampleTable$Storm <- as.factor(sampleTable$Storm )
        sampleTable$Tissue <- as.factor(sampleTable$Tissue )
        sampleTable$LH_Sub_Stage <- as.factor(sampleTable$LH_Sub_Stage)
        sampleTable$Conditions <- as.factor(sampleTable$Conditions)
        sampleTable$Time_Code <- gsub("\\:", "", sampleTable$Time)
        sampleTable$Time_Code <- as.factor(sampleTable$Time_Code)
        head(sampleTable)
        row.names(sampleTable) <- sampleTable$ID #make the rownames= ID
      }
    }
    #Analyses in container
    #-----------------
    #2.Read in featurecounts output 
    #-----------------
    #a <- fread("M:/ROSLIN/RNA_2LALO/psiclass08022021_vote_Modified_Gene_name.count",header=T)
    #a <- fread("M:/ROSLIN/RNA_2LALO/psiclass08022021_vote_Modified_Gene_name.count_Galgal7",header=T)
    
    #dont need to read in everytime
    a <- fread("M:/ROSLIN/RNA_2LALO/Output_samples.txt_Galgal7_unmapped",header=T)
    meta <- a[,1:6]
    cts <- a[,7:ncol(a)]
    #head(cts)
    #change the ID od columns in count table
    {id<- names(cts) #change the ID od columns
      id <- str_split(id, "/",simplify = T)[,14]
      id <- str_split(id, "Al",simplify = T)[,1]
      names(cts) <- id}
    cts <- as.matrix(cts)
    row.names(cts)=meta$Geneid
    head(cts)
    rm(a)
    #Make sure the order in colsata is the same as in count file
    cts<- cts[,rownames(sampleTable)]#cts has to be matrix
    all(rownames(sampleTable) == colnames(cts))
    dim(cts)
    #-----------------
    #3. Run DESeq2 
    #-----------------
    dds <- DESeqDataSetFromMatrix(countData=cts,
                                  colData=sampleTable,
                                  design =~ Conditions) # ~ Stage+Sex+pheno|For stage DEGs:Sex+Stage
    #dds is now ready for DESeq() see DESeq2 vignette
    keep <- rowSums(counts(dds)) >= 36#a minimal pre-filtering to keep only rows that have at least 96 reads total.
    table(keep)
    dds <- dds[keep,]
    dim(dds)
    
    exprSet=assay(dds)
    #Summary a table with the number of gene expressed (non zero) in each sample
    SAMPLEgene <- as.data.frame(apply(exprSet,2, function(x) sum(x>0)))
    names(SAMPLEgene) <- paste0("expressed_gene_count_out_of_",dim(dds)[1])
    SAMPLEgene$Per <- round(SAMPLEgene[,1]/dim(dds)[1],digits=5)
    write.csv(SAMPLEgene,paste0("M:/ROSLIN/RNA_2LALO/Filter_",TS,"_Ind_expressed_gene_count.csv"),row.names = T)
    
    #  { #Optional Filter for genelist e.g. GWAS genes/DEGs
    #      geneSet <- read.table("M:/ROSLIN/RNA-seq/F_quantification/E5/Breeds/DPTS//GENE_DEseq2_Vol_E5_DPTS_Volcano_ALL.txt",header=F,col.names="gene_id")
    #    head(geneSet)
    #    keep_geneset <- (rownames(exprSet)%in%geneSet$gene_id)
    #    table(keep_geneset)
    #    dds <- dds[keep_geneset,]
    #  }
    dim(dds)
    dds <- DESeq(dds,minReplicatesForReplace=Inf)###for schijd
    
    #  ###DEGs between PHENOTYPE
    #source('M:/ROSLIN/RNA_2LALO/Codes/My_RNA_plot_functions.R')
    ###########Gene of interest for highlight
    
    GeneOfInterest<-c("FKBP5","COL","PER1","PER2","PER3","SOUL","TBC1D2","CRY1","FMN2","ZP3","TBC1D8",'AMBP','PCK1','UCN3','FBXO32','CPS1','NR1D1','NPC1')
    #GeneOfInterest<-read.table('M:/ROSLIN/Klebsiella_variicola_RNA/Gene_of_interest.txt',header=F)$V1
    
    
    ####Interate between conditions
    # for (i in 1:length(dat2$x)){
    #   i=1
    # Base=dat2$x[i]
    # Case=dat2$y[i]
    Base=3
    Case=4
     print(paste('Control:', Base,
                 ';Case:', Case)) 

    
    Pair=paste(Case,Base,sep="_")
    dds$condition <- factor(dds$Conditions, levels = c(Case,Base))
    dds$condition <- relevel(dds$condition, ref = as.character(Base)) #as.char, otherwise give error in this command
    res <- results(dds,contrast=c("Conditions",Case,Base)) #Contrast: The level given last is the base level for the comparison.
    resultsNames(dds)
    res <- lfcShrink(dds,contrast=c("Conditions",Case,Base), res=res, quiet=T,type='normal')####Shrink log2 fold changes
    # resNorm <- lfcShrink(dds, coef=2, type="normal")
    
    # xlim <- c(1,1e5); ylim <- c(-2,2)
    # plotMA(res, xlim=xlim, ylim=ylim, main="lfs") Plot the MA,shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet
    # plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
    
    # # ###  DEGs between STAGES
    # dds$condition <- factor(dds$Stage, levels = c("13","5")) ###DEGs between stages
    # dds$condition <- relevel(dds$condition, ref = "5")###DEGs between stages
    # res <- results(dds,contrast=c("Stage","13","5"))###DEGs between stages
    # res <- lfcShrink(dds,contrast=c("Stage","13","5"), res=res, quiet=T) ###DEGs between stages
    #######  E5 is used as the base for DEG
    # #plotMA(res, ylim=c(-5,5))
    head(res)
    
    resOrdered <- res[order(res$padj),]#Sort summary list by p-value
    head(resOrdered)
    DEG =as.data.frame(resOrdered)
    DEG = na.omit(DEG)
    
    # 
    # ####Output 1: all colums
    {DEseq_DEG <- DEG
      nrDEG=DEseq_DEG #[,c(2,6)]
      nrDEG$Gene <- rownames(nrDEG)
      #colnames(nrDEG)=c('log2FoldChange','padj')
      #head(nrDEG)
      prefix=paste0(TS,"_",Pair)
      print(prefix)
      write.table(nrDEG,paste0("1DESeq2_output_",prefix,".txt"),row.names = F,quote = F,sep="\t")
    }
    
    # #OUtput 2: Plot volcano and UP/Down regulated
    #source('M:/ROSLIN/RNA_2LALO/Codes/My_RNA_plot_functions.R')
    #source('M:/ROSLIN/RNA_2LALO/Codes/My_RNA_plot_functions_simple.R')
    {
      DEseq_DEG <- DEG
      DEseq_DEG=DEseq_DEG[,c(2,6)] #5: unadjusted P-value 6:Padj
      colnames(DEseq_DEG)=c('log2FoldChange','padj')
      #row.names(DEseq_DEG) <- str_split(row.names(DEseq_DEG), "[.]",simplify = T)[,1] #Change gene ID
      prefix=paste(TS,"Conditions",Case,"VS",Base,sep="_")
      draw_v_U_D(exprSet,DEseq_DEG,paste0(prefix))
    }
      
      # #rm(list=ls())
    }#Interate between combination
  }#Ieration between Tissue
#}#First if




