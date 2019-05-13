library(DESeq2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pasilla)
library(ggplot2)
library(dplyr)
library(topGO)
library(tidyverse)
library(reshape2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(grid)
require(gridExtra)
library(RColorBrewer)
#library(ggz


setwd("~/dataOS/CS_RNA/rePloting")

#----------------Data Preperaion-----------------------

FPKM_info_all <- read.csv("~/dataOS/CS_RNA/surface_click/NK92/20190408/20190408_FPKMs.txt",
                          sep="\t") %>% as.data.frame(.)
#annotation info
?select

#biotypes
Biotype_df <- FPKM_info_all[,c("gene_id","Biotype")]

#names
cols <- c("GENENAME","SYMBOL")
anno_info <-AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(Biotype_df$gene_id), columns=cols, keytype="ENSEMBL")

#columns(org.Hs.eg.db)
############## Surface Clieck data
setwd("~/dataOS/CS_RNA/surface_click/NK92/20190408/")
RNAclick_counts <- cbind ((read.csv("htseqCountTotal.txt",sep = "\t")),
                          (read.csv("~/dataOS/CS_RNA/surface_click/NK92/20190409/htseqCountTotal.txt",sep = "\t")[,-1]))
                    
                
colnames(RNAclick_counts) <- c("gene_id","NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2",
                               "NK92_Tot_Click_1","NK92_Tot_Click_3","NK92_Sur_Click_2","NK92_Tot_Click_3_1")
rownames(RNAclick_counts) <- RNAclick_counts$gene_id 
RNAclick_counts$gene_id <- NULL

#Process RNAclick_count data 
tt <- RNAclick_counts %>% as.matrix(.) %>% rowMaxs(.) 
  (rowMaxs(RNAclick_counts)) %>% as.data.frame(.)
summary(tt)

hist(log10(tt))

## subset by filterng the rows with row sum smaller than 7
RNAclick_select_counts <- subset(RNAclick_counts,rowSums(RNAclick_counts)>=7)

#----------overall Analysis-------------
setwd("~/dataOS/CS_RNA/surface_click/NK92/")
# Build the dds object 
RNAclick_select_counts_temp <- RNAclick_select_counts[,-7]
dds_tot <- DESeqDataSetFromMatrix(countData = as.matrix(RNAclick_select_counts_temp),
                                  colData = data.frame(row.names = colnames(as.matrix(RNAclick_select_counts_temp)),
                                                       condition = c("Surface","Control","Control","Total","Total","Surface")),
                                  design=~condition)
# Plot PCA

require("ggrepel")
z <- plotPCA(rlog(dds_tot))+geom_label_repel(aes(label = name))+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18,face="bold"),
        plot.title = element_text(size=15,face="bold"))+
  #guides(colour = guide_legend(override.aes = list(alpha = 1/4,size = 5,shape = 1)),size = FALSE)+
  ggtitle("PCA of all the Surface_click Samples")

z  

ggsave("PCA_20190408.png",z,height = 6 , width = 10)
plotDispEsts(dds_tot)
#-------------- NK92_Sur_Click_1 VS Negative controls ---------------
#set wd
setwd("~/dataOS/CS_RNA/surface_click/NK92/Pair_wise_comp/NK92_1_V_NCs")
SC_df <- format_samples(RNAclick_counts,c("NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2"))
DEres_SC <- perform_DEseq(SC_df,1,2)
#cand_SC <- subset(DEres_SC,(log2FoldChange>1&padj<0.05)|(log2FoldChange< -1&padj<0.05))%>% 
#  as.data.frame(.) 
cand_SC <- subset(DEres_SC,padj<0.05)

#-------------NK92_Tot_3 VS negative control---------
SC_df <- format_samples(RNAclick_counts,c("NK92_Tot_Click_3_1","NK92_Tot_Click_3"))
DEres_SC <- perform_DEseq(SC_df,1,1)
cand_SC <- subset(DEres_SC,(log2FoldChange>2&padj<0.05)|(log2FoldChange< -2&padj<0.05))%>% 
  as.data.frame(.) 
#test
cand_SC <- subset(DEres_SC,log2FoldChange>6) %>% as.data.frame() 
#temp_geneList <- cbind(cand_SC$Gene_name,cand_SC$log2FoldChange) %>% as.data.frame(.)

# 
# gene_ID <- cand_SC$Gene_name %>% as.character(.)
# cols <- c("GENENAME")
# anno_info <- select(org.Hs.eg.db, keys=gene_ID, columns=cols, keytype="ENSEMBL")
# 
# #format output dataframe
# res <- merge (cand_SC,anno_info, by.x = "Gene_name", by.y = "ENSEMBL")
# res_final <- merge (res,Biotype,by.x = "Gene_name",by.y = "gene_id")
# write.csv(res_final,"temp_gene.csv")
  #as.data.frame(DEres_SC) 
#length(which(cand_SC$padj== 0))
length(which(is.na(cand_SC$padj)&is.na(cand_SC$log2FoldChange)))
dim(cand_SC)

# -------------- NK92_Sur_Click_1+Two Negative control VS NK92_Sur_Click_2 -------------
#set wd
setwd("~/dataOS/CS_RNA/surface_click/NK92/Pair_wise_comp/NK92_Sur_1_V_NK92_Sur_2/")
SC_df <- format_samples(RNAclick_counts,c("NK92_Sur_Click_2","NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2"))
DEres_SC <- perform_DEseq(SC_df,1,3)

cand_2 <- subset(DEres_SC,(log2FoldChange>2&padj<0.05)) %>% as.data.frame(.)
cand_1 <- subset(DEres_SC,(log2FoldChange< -2&padj<0.05)) %>% as.data.frame(.)

gene_ID_2 <- cand_2$Gene_name %>% as.character(.)
gene_ID_1 <- cand_1$Gene_name %>% as.character(.)

cols <- c("GENENAME")
anno_info <- select(org.Hs.eg.db, keys=gene_ID, columns=cols, keytype="ENSEMBL")

#format output dataframe
res_1 <- merge (cand_1,anno_info, by.x = "Gene_name", by.y = "ENSEMBL")
res_1_final <- merge(res_1,Biotype_df,by.x  = "Gene_name",by.y = "gene_id")

res_2 <- merge (cand_2,anno_info, by.x = "Gene_name", by.y = "ENSEMBL")
res_2_final <- merge(res_2,Biotype_df,by.x  = "Gene_name",by.y = "gene_id")

tt <- table(res_1_final$Biotype) %>% as.data.frame(.)

write.csv(res_2_final,"temp_gene_2.csv")
write.csv(res_1_final,"temp_gene_1.csv")

#----------------- Sur_1&2 Vs Tot-----------
setwd("~/dataOS/CS_RNA/surface_click/NK92/Pair_wise_comp/NK92_Sur_all_V_NK92_Tot/")
SC_df <- format_samples(RNAclick_counts,c("NK92_Sur_Click_2","NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2",
                                          "NK92_Tot_Click_1","NK92_Tot_Click_3"))
DEres_SC <- perform_DEseq(SC_df,4,2)
res_temp<- get_DEgenes_info(DEres_SC,2,0.05,"all_sur_V_tot.csv")


SC_df <- format_samples(RNAclick_counts,c("NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2",
                                          "NK92_Tot_Click_1","NK92_Tot_Click_3"))
DEres_SC <- perform_DEseq(SC_df,3,2)
res_temp<- get_DEgenes_info(DEres_SC,2,0.05,"all_sur_1_V_tot.csv")


#-----------Sur_2 Vs Tot---------------
SC_df <- format_samples(RNAclick_counts,c("NK92_Sur_Click_2",
                                          "NK92_Tot_Click_1","NK92_Tot_Click_3"))
DEres_SC <- perform_DEseq(SC_df,1,2)

res_temp<- get_DEgenes_info(DEres_SC,2,0.05,"all_sur_2_V_tot.csv")

#get the dataframe of the normalized count 
DEres_SC_dds <- perform_DEseq_dds(SC_df,1,2)
Norm_counts <- (counts(DEres_SC_dds)/sizeFactors(DEres_SC_dds)) %>% as.data.frame(.)
#with psuedo counts
Norm_counts$NK92_Tot <- ((rowSums(Norm_counts[,2:3]))/2)+1
Norm_counts$NK92_Tot_Click_1 <- NULL
Norm_counts$NK92_Tot_Click_3 <- NULL
Norm_counts <-rownames_to_column(Norm_counts,"Gene_name")

#plot
setwd("~/dataOS/CS_RNA/surface_click/NK92/Pair_wise_comp/")
plot_cands <- merge (res_temp,Norm_counts,by = "Gene_name")%>% .[,c(1,3,9,11,12)]
write.csv(plot_cands,"cand_genes.csv",sep = "\t",row.names = FALSE)

plot_cands_withFunc <- read.csv("Func_annotation.csv",sep = "\t")
#plot_cands_withFunc <- cbind(plot_cands_withFunc,plot_cands[,2])
#pp <- ggplot(data = temp_df,aes(x = log2(NK92_Sur_Click_2),y = log2(NK92_Tot),label = SYMBOL,size = log2FoldChange))

#------plot with NonFunc Genes-------------
remove(total_df)
for (i in (6:ncol(plot_cands_withFunc))){
  #print (i)
  temp_df <- plot_cands_withFunc[,c(1,2,3,4,5,i)] %>% subset(.,!is.na(.[,6]))
  colnames(temp_df)[6]<-"func_tag" 
  if (!exists("total_df")){
    total_df <- temp_df
  }
  else{
    total_df <- rbind(total_df,temp_df)
  }
}
# NA_df <- subset(plot_cands_withFunc,rowSums(is.na(plot_cands_withFunc[,6:17])) == 12) %>% .[,c(1,2,3,4,5)]
# NA_df$func_tag <- "None"
# total_df <- rbind(total_df,NA_df)
#rowSums(is.na(plot_cands_withFunc[,6:17]))

total_df$func_tag <- factor(total_df$func_tag,levels =
                                   c("Transmembrane helix","Transmembrane","Cell membrane","Membrane",
                                         "Hydrogen ion transport","Signal","Glycoprotein","Disulfide bond",
                                         "Retinitis pigmentosa","CF(0)","ATP synthesis","Leigh syndrome"))
                                        
#colorCount = ncol(plot_cands_withFunc)-5
getPalette = colorRampPalette(brewer.pal(11, "RdYlGn"))
colors = getPalette(12)

#c(colors,"black")
pp <- ggplot(data = total_df,aes(x = log2(NK92_Sur_Click_2),y = log2(NK92_Tot),
                           label = SYMBOL,size = log2FoldChange,color = func_tag))+
  geom_point(position="jitter",alpha = 1/3)+
  scale_color_manual(values = colors)+theme_classic2()+
  scale_size_continuous(range = c(3, 12))+xlab("log2(K562 Surface Sample Gene Count)")+
  ylab("log2(K92 Total Sample Gene Count)")+ggtitle("Function Distribution")+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18,face="bold"),
        plot.title = element_text(size=18,face="bold"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size = 5)),size = FALSE)
  
pp

ggsave("Function_dotPlot.png",pp,height = 6 , width = 10)


ppp <- ggplot(data = plot_cands_withFunc,aes(x = log2(NK92_Sur_Click_2),y = log2(NK92_Tot),
                                  label = SYMBOL,size = log2FoldChange))+
  geom_point(position="jitter",alpha = 1/4)+
  scale_color_manual(values = colors)+theme_classic2()+
  scale_size_continuous(range = c(3, 12))+xlab("log2(NK92 Surface Sample Gene Count)")+
  ylab("log2(NK92 Total Sample Gene Count)")+ggtitle("Gene Distribution")+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18,face="bold"),
        plot.title = element_text(size=15,face="bold"))

ppp

ggsave("Function_dotPlot_all.png",ppp,height = 6 , width = 10)
# position_jitter(h=0.15,w=0.15)
#pp
# temp_df <- plot_cands_withFunc[,c(1,2,3,4,5,6)] %>% subset(.,!is.na(.[,6]))
# total_df <- rbind(total_df,temp_df)
# +
#   geom_point(color = colors[1],alpha = 1/4)





+geom_label_repel(aes(label = SYMBOL))
  


#getPalette = colorRampPalette(brewer.pal(12, "Set3"))
#geom_text(aes(label=SYMBOL),hjust=0, vjust=0)



# DEres_SC <- perform_DEseq(SC_df,1,1)
# condition <- factor(c(rep("Surface",1), rep("Total", 1)),levels = c("Total","Surface"))
# col_data <- data.frame(row.names = colnames(as.matrix(SC_df)),condition)
# dds <- DESeqDataSetFromMatrix(countData = as.matrix(SC_df),colData = col_data,design=~condition)
# 
# #perform DEseq without pre-filtering
# dds <- DESeq(dds,betaPrior=FALSE)
# res_dds <- results(dds) %>% as.data.frame(.) %>% rownames_to_column(.,"Gene_name")




#z <- plotPCA(rlog(dds))
#nudge <- position_nudge(x=2,y = 4)
#z + geom_label(aes(label = name), position = nudge)



#----------------------------Functions-----------------------------
# Method for formating df into required format
# format: surface samples then total samples, accordingly as vector
# input: DataFrame
# output: formated Dataframe
format_samples <- function (df,order_vec){
  df <- df[order_vec]
  return (df)
}

# Method for performing DESeq analysis
# Input: Method df
#        surface sample number 
#        total sample numMembers vectors, in the form of c(x,y), x = number of surface samples, y = number of total samples
# Output: Method result df
perform_DEseq <- function (df,ss,ts){
  #setting conditions for DEseq of two dfs
  condition <- factor(c(rep("Surface",ss), rep("Total", ts)),levels = c("Total","Surface"))
  col_data <- data.frame(row.names = colnames(as.matrix(df)),condition)
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(df),colData = col_data,design=~condition)
  
  #perform DEseq without pre-filtering
  dds <- DESeq(dds,betaPrior=FALSE)
  res_dds <- results(dds) %>% as.data.frame(.) %>% rownames_to_column(.,"Gene_name")
  
  return(res_dds)
}

perform_DEseq_dds <- function(df,ss,ts){
  #setting conditions for DEseq of two dfs
  condition <- factor(c(rep("Surface",ss), rep("Total", ts)),levels = c("Total","Surface"))
  col_data <- data.frame(row.names = colnames(as.matrix(df)),condition)
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(df),colData = col_data,design=~condition)
  
  #perform DEseq without pre-filtering
  dds <- DESeq(dds,betaPrior=FALSE)
  #res_dds <- results(dds) %>% as.data.frame(.) %>% rownames_to_column(.,"Gene_name")
  
  return(dds)
}




get_DEgenes_info <- function(DEres_df,log2FoldChaneg_thresh,Pval_thresh,fileName){
  cand_names<- subset(DEres_df,(log2FoldChange>log2FoldChaneg_thresh&
                                  padj<Pval_thresh)) %>% as.data.frame(.)
  gene_ID <- cand_names$Gene_name %>% as.character(.)
  
  res_df <- merge (cand_names,anno_info, by.x = "Gene_name", by.y = "ENSEMBL")
  res_df_final <- merge(res_df,Biotype_df,by.x  = "Gene_name",by.y = "gene_id")
  
  write.csv(res_df_final,fileName)
  return (res_df_final)
}
