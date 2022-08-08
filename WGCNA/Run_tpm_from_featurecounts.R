# Convert featurecounts to TPM
a <- fread("M:/ROSLIN/RNA_Seq_LALO/Output_samples.txt_Galgal7_unmapped",header=T)
meta <- a[,1:6]
cts <- a[,7:ncol(a)]
#rm(a)
#head(cts)

#change the ID od columns in count table
{id<- names(cts) #change the ID od columns
  id <- str_split(id, "/",simplify = T)[,14]
  id <- str_split(id, "A",simplify = T)[,1]
  names(cts) <- id}
cts <- as.matrix(cts)
row.names(cts)=meta$Geneid

#Make sure the order in colsata is the same as in count file
sampleTable<- read.table("M:/ROSLIN/RNA_Seq_LALO/WGCNA/LALO_trait.txt", header = TRUE,stringsAsFactors = T)
# Prepare the class  of columns
{  sampleTable$Storm <- as.factor(sampleTable$Storm )
  sampleTable$Tissue <- as.factor(sampleTable$Tissue )
  sampleTable$LH_Sub_Stage <- as.factor(sampleTable$LH_Sub_Stage)
  sampleTable$Conditions <- as.factor(sampleTable$Conditions)
  sampleTable$Time_Code <- as.factor(sampleTable$Time_Code)
  sampleTable$Time_Code <- gsub("\\:", "", sampleTable$Time)
  head(sampleTable)
  row.names(sampleTable) <- sampleTable$ID #make the rownames= ID
}
cts<- cts[,rownames(sampleTable)]#cts has to be matrix
all(rownames(sampleTable) == colnames(cts))
dim(cts)
head(cts)
head(meta)
Correct_featurecounts <- cbind(meta,cts)
head(Correct_featurecounts)
source("M:/ROSLIN/RNA_Seq_LALO/Codes/Cal_TPM.R")
tpm_featurecounts <- calc_tpm_from_featurecounts(Correct_featurecounts)
test <- head(tpm_featurecounts)
ncol(test)
tpm_featurecounts_strip <- tpm_featurecounts[,c(1:6,55:102)]
head(tpm_featurecounts_strip)
{id<- names(tpm_featurecounts_strip) 
id_1 <- id[1:6]
id_2 <- str_split(id[7:length(id)], "m_",simplify = T)[,2]
id <- c(id_1,id_2)}
all(names(Correct_featurecounts) == id)
all(rownames(sampleTable) == id[7:length(id)])

names(tpm_featurecounts_strip) <- id
write.table(tpm_featurecounts_strip,"M:/ROSLIN/RNA_Seq_LALO/WGCNA/LALO_TPM_genecounts.txt",row.names = F,col.names = T,sep="\t",quote = F)
