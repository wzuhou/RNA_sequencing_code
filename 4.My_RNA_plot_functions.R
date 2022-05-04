#draw_v_U_D(exprSet,DEseq_DEG,paste0(prefix))
# need_DEG=DEseq_DEG
# n=paste0(prefix)
###########
#FUNCTION##
###########
draw_v_U_D <- function(exprSet,need_DEG,n='DEseq2'){
  ## we only need two columns of DEG, which are log2FoldChange and pvalue
    
  #logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
  logFC_cutoff=0.58
  
  need_DEG$change = as.factor(ifelse(need_DEG$padj < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  need_DEG$Gene <- row.names(need_DEG)
  write.table(subset(need_DEG,change=="UP"),
              paste0(n,"_Volcano_UP.txt"),row.names=F,quote=F,sep="\t")
  write.table(subset(need_DEG,change=="DOWN"),
              paste0(n,"_Volcano_DOWN.txt"),row.names=F,quote=F,sep="\t")
 write.table(row.names(subset(need_DEG,change=="UP")),
              paste0("GENE_",n,"_Volcano_UP.txt"),row.names=F,quote=F,sep="\t",col.names=F)
  write.table(row.names(subset(need_DEG,change=="DOWN")),
              paste0("GENE_",n,"_Volcano_DOWN.txt"),row.names=F,quote=F,sep="\t",col.names=F)
  write.table(row.names(subset(need_DEG,change=="UP" |change=="DOWN")),
              paste0("GENE_",n,"_Volcano_ALL.txt"),row.names=F,quote=F,sep="\t",col.names=F)

  this_tile <- paste0(n,
                      ### '\nCutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
  )
  library(ggplot2)
  library(ggrepel)
  #########Optional: highlight & label
  need_DEG<-need_DEG %>%
    # Add highlight and annotation information
     mutate( is_highlight=ifelse( -log10(padj)>max(-log10(padj))*0.95 ,"yes", "no")) %>%
    # mutate( is_annotate=ifelse(Gene %in% GeneOfInterest|-log10(padj)>max(-log10(padj))*0.95 , "yes", "no"))
    mutate( is_annotate=ifelse(row_number()<=8 , "yes", "no"))
  #mutate( is_annotate=ifelse(-log10(padj)>max(-log10(padj))*0.95 , "yes", "no"))
  #########OptionaL#####################
  g = ggplot(data=need_DEG, 
             aes(x=log2FoldChange, y=-log10(padj), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    ####################OPTIONAL
    # Add highlighted points
    geom_point(data=subset(need_DEG, is_highlight=="yes"), color="firebrick4", size=2,alpha=0.8) +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=subset(need_DEG, is_annotate=="yes"), aes(label=Gene), size=2) +
    ####################OPTIONAL
    theme_set(theme_set(theme_classic(base_size=15)))+
    xlab("log2 fold change") + ylab("-log10 FDR P-value") +
    ggtitle( this_tile ) +
    theme(plot.title = element_text(size=13,hjust = 0.5))+
    #scale_colour_manual(values = c('dodgerblue4','black','firebrick4')) ## corresponding to the levels(res$change)+
    scale_colour_manual(values = c('#3B9AB2','black','#EBCC2A'))+
    theme_classic()
  #print(g)
  ggsave(g,filename = paste0("DEseq2_Vol_",n,'_volcano_Label_Highlight.png'),width = 4.5,height = 4.5)
}
