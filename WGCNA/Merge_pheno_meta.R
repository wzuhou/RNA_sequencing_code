library(plyr)
df1<-read.table('M:/ROSLIN/RNA_Seq_LALO/Sample_strip_ID_Tissue(Time).txt',header=T)
df <- read.table('M:/ROSLIN/RNA_Seq_LALO/1Phenotype/Phenotype_filtered.txt',header = T)

head(df)
df<-as.data.frame(df)

df$"LH_Sub_Stage"<-as.factor(df$"LH_Sub_Stage")
df$"Testicular_Volume"<-1*df$"Testicular_Volume"
df$CP_Height <- as.numeric(df$CP_Height)
df$CP_Volume=(4/3)*pi*((0.5*df$CP_Width)^2)*(0.5*df$CP_Height)
df2<-df[,c(1,11:29,34)]
names(df2)[1]<-'Num'

df3 <- join(df1,df2,by='Num',type='left')

for (i in seq(ncol(df3))){
print(paste(i,all(!is.na(df3[,i]))))
}
names(df3)
df3<-  df3[,c(1:3,8,10:12,14,16,20:27,34:37,40)]
write.table(df3,'M:/ROSLIN/RNA_Seq_LALO/WGCNA/LALO_trait.txt',quote = F, sep='\t',row.names = F)
