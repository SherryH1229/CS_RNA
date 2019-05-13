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
#library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(grid)
require(gridExtra)
library(ggfortify)


setwd("~/dataOS/CS_RNA/rePloting")

#----------------Data Preperaion-----------------------

sample_info_1 <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_seq/K562/20190412/20190412_FPKMs.txt",sep="\t")
sample_info_2 <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_seq/K562/20190421/20190421_FPKMs.txt",sep  = "\t")
sample_info_1$X <- NULL
sample_info_1$K562_TotalRNA_FPKM <- NULL
FPKM_info_all <- cbind(sample_info_1,sample_info_2$K562_TotalRNA_FPKM) %>% .[!duplicated(.$gene_id),]
colnames(FPKM_info_all)[8:14] <- c("K562_Seq_Sur_2","K562_Seq_Sur_4","K562_Seq_Sur_1","K562_Seq_Sur_3","K562_Seq_Sur_5",
                                 "K562_Seq_Tot_2","K562_Seq_Tot_1")


#annotation info

#biotypes
Biotype_df <- FPKM_info_all[,c("gene_id","Biotype")]
gene_ID <- as.vector(Biotype_df$gene_id)
#names
cols <- c("GENENAME","SYMBOL")
anno_info <- select(org.Hs.eg.db, keys=gene_ID, columns=cols, keytype="ENSEMBL")
############## Surface seq data
setwd("~/dataOS/CS_RNA/surface_seq/K562/20190412/")
RNAseq_counts <- cbind ((read.csv("htseqCountTotal.txt",sep = "\t")),
                          (read.csv("~/dataOS/CS_RNA/surface_seq/K562/20190421/htseqCountTotal.txt",sep = "\t")[,-1]))
RNAseq_counts$K562_TotalRNA <- NULL             
                
colnames(RNAseq_counts) <- c("gene_id","K562_Seq_Sur_4","K562_Seq_Sur_2","K562_Seq_Sur_3",
                               "K562_Seq_Sur_1","K562_Seq_Sur_5","K562_Seq_Tot_1","K562_Seq_Tot_2")

rownames(RNAseq_counts) <- RNAseq_counts$gene_id 
RNAseq_counts$gene_id <- NULL

#Process RNAclick_count data 
# tt <- RNAclick_counts %>% as.matrix(.) %>% rowMaxs(.) 
#   (rowMaxs(RNAclick_counts)) %>% as.data.frame(.)
# summary(tt)
# 
# hist(log10(tt))
# 
# ## subset by filterng the rows with row sum smaller than 7
# RNAclick_select_counts <- subset(RNAclick_counts,rowSums(RNAclick_counts)>=7)


############### SurfaceSeq data
# RNAseq_counts <- read.csv("count_matrix.txt",sep = "\t", head = TRUE,skip=1,row.names = "Geneid")

#----------overall Analysis-------------
setwd("~/dataOS/CS_RNA/surface_seq/K562/Deseq/Deseq_cands/")
# Build the dds object 
dds_tot <- DESeqDataSetFromMatrix(countData = as.matrix(RNAseq_counts),
                                  colData = data.frame(row.names = colnames(as.matrix(RNAseq_counts)),
                                                       condition = c("Surface","Surface","Surface","Surface","Surface","Total","Total")),
                                  design=~condition)

# Plot PCA

require("ggrepel")
z <- plotPCA(rlog(dds_tot))+geom_label_repel(aes(label = name))+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18,face="bold"),
        plot.title = element_text(size=15,face="bold"))+
  ggtitle("PCA of all the Surface_seq Samples")

z  

ggsave("PCA_20190408.png",z,height = 6 , width = 10)
#plotDispEsts(dds_tot)
#-------------- DEseq Surface_seq ---------------
#set wd
#setwd("~/dataOS/CS_RNA/surface_click/NK92/Pair_wise_comp/NK92_1_V_NCs")

SS_df <- format_samples(RNAseq_counts,c("K562_Seq_Sur_4","K562_Seq_Sur_2","K562_Seq_Sur_3",
                                          "K562_Seq_Sur_1","K562_Seq_Sur_5","K562_Seq_Tot_1","K562_Seq_Tot_2"))
DEres_SS <- perform_DEseq(SS_df,5,2)
res_temp<- get_DEgenes_info(DEres_SS,2,0.05,"DEseq_cands_SurfaceSeq.csv")

####Tot_1
SS_df <- format_samples(RNAseq_counts,c("K562_Seq_Sur_4","K562_Seq_Sur_2","K562_Seq_Sur_3",
                                        "K562_Seq_Sur_1","K562_Seq_Sur_5","K562_Seq_Tot_1"))
DEres_SS <- perform_DEseq(SS_df,5,1)
res_temp<- get_DEgenes_info(DEres_SS,2,0.05,"DEseq_cands_SurfaceSeq_NoRibDep.csv")

####Tot_2
SS_df <- format_samples(RNAseq_counts,c("K562_Seq_Sur_4","K562_Seq_Sur_2","K562_Seq_Sur_3",
                                        "K562_Seq_Sur_1","K562_Seq_Sur_5","K562_Seq_Tot_2"))
DEres_SS <- perform_DEseq(SS_df,5,1)
res_temp<- get_DEgenes_info(DEres_SS,2,0.05,"DEseq_cands_RibDep.csv")


#get the dataframe of the normalized count 
DEres_SS_dds <- perform_DEseq_dds(SS_df,5,2)
Norm_counts <- (counts(DEres_SS_dds)/sizeFactors(DEres_SS_dds)) %>% as.data.frame(.)

#with psuedo counts
Norm_counts$K562_Tot <- ((rowSums(Norm_counts[,6:7]))/2)+1
Norm_counts$K562_Sur <- ((rowSums(Norm_counts[,1:5]))/2)+1
Norm_counts[,1:7]<- NULL
Norm_counts <-rownames_to_column(Norm_counts,"Gene_name")

plot_cands <- merge (res_temp,Norm_counts,by = "Gene_name")%>% .[,c(1,3,9,11,12)]
write.csv(plot_cands,"cand_genes.csv",sep = "\t",row.names = FALSE)

plot_cands_withFunc <- read.csv("Func_annotation.csv",sep = "\t")


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


cols <-c("all_genes" = "black","cellSurfaceFunc_genes" = "red")
ppp <- ggplot(data = plot_cands_withFunc,aes(x = log2(K562_Sur),y = log2(K562_Tot),
                                             label = SYMBOL,size = log2FoldChange))+
  geom_point(aes(color = "all_genes"),alpha = 1/4)+
  geom_point(data = total_df,aes(color = "cellSurfaceFunc_genes"),alpha = 1/2)+
  geom_label_repel(data = total_df,aes(label = SYMBOL,size = 2))+
  #geom_point(total_df,cl)
  scale_color_manual(values = colors)+theme_classic2()+
  scale_size_continuous(range = c(3, 12))+xlab("log2(K562 Surface Sample Gene Count)")+
  ylab("log2(K562 Total Sample Gene Count)")+ggtitle("Gene Distribution")+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18,face="bold"),
        plot.title = element_text(size=15,face="bold"))+
  scale_color_manual("Gene_type",values = cols)+
  #xlim(0,15)+ylim(0,15)+
  guides(colour = guide_legend(override.aes = list(alpha = 1/4,size = 5)),size = FALSE)
  #guides(colour = guide_legend(size = FALSE))

ppp
setwd("~/dataOS/CS_RNA/surface_seq/K562/Deseq/Deseq_cands/")
ggsave("Function_dotPlot.png",ppp,height = 6 , width = 10)
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

res_dds <- results(dds) %>% as.data.frame(.) %>% rownames_to_column(.,"Gene_name")

#total of 21 genes
crossed_genes <- get_crossed_geneList(DEres_SC,DEres_SS,FPKM_info_all,2,0.05)

g_foldChange <- plot_FoldChange(DEres_SC,DEres_SS,crossed_genes,fish_string,
                                fish_sig_string,"log2FoldChange_SurfaceClick",
                                "log2FoldChange_SurfaceSeq")
ggsave("Scatter_plot_FoldChange.png",g_foldChange,height = 6 , width = 10)

Volcano_plots <- plot_Volcanos(DEres_SC,"Surface_click_genes","Candidate_SS_genes_inSC",DEres_SS,"Surface_seq_genes","Candidate_SC_genes_inSS")
ggsave ("Volcano_plots.png",Volcano_plots,height = 6.8, width = 20)

#----------------------------Functions-----------------------------
# Method for formating df into required format
# format: surface samples then total samples, accordingly as vector
# input: DataFrame
# output: formated Dataframe
format_samples <- function (df,order_vec){
  df <- df[order_vec]
  return (df)
}

# Method for plotting two methods global Surface_mean over Total Mean
# Input format: FirstMethod dataframe
#               FirstMethod surface sample colnums as vector of nums
#               FirstMethod total sample Colnums as vector of nums
#               FirstMethod gene tag, tag that shows on the legend 
#               SecondMethod DataFrame
#               SecondMethod surface sample colnums as vector of nums
#               SecondMethod total sample Colnums as vector of nums
#               SecondMethod gene tag, tag that shows on the legend
# Output: ggplot plot scatter plot object with linear regression line and Pearson's correlation
plot_Mean <- function (SC_df,surf_vec_SC,tot_vec_SC,SC_tag,SS_df,surf_vec_SS,tot_vec_SS,SS_tag){
  #prepare FC tables 
  FC_df_SC <- as.data.frame(rowMeans(SC_df[,surf_vec_SC]+1)) %>% cbind(.,rowMeans(SC_df[,tot_vec_SC]+1))
  FC_df_SC$Tag <- SC_tag
  colnames(FC_df_SC) <- c("SurfaceMean","TotalMean","Tag")
  
  #surface seq data
  FC_df_SS <- as.data.frame(rowMeans(SS_df[,surf_vec_SS]+1)) %>% cbind(.,rowMeans(SS_df[,tot_vec_SS]+1))
  FC_df_SS$Tag <- SS_tag
  colnames(FC_df_SS) <- c("SurfaceMean","TotalMean","Tag")
  
  totol_mean_df <- rbind(FC_df_SC,FC_df_SS)
  #get pearson's correlation and format tags
  correltion_SS <- cor.test(FC_df_SS$SurfaceMean, FC_df_SS$TotalMean, method = "pearson", conf.level = 0.95)
  correltion_SC <- cor.test(FC_df_SC$SurfaceMean, FC_df_SC$TotalMean, method = "pearson", conf.level = 0.95)
  
  cor_tag_SS <- paste (paste(SS_tag," Pearson's Correlation: ",sep = ""),
                       format(round(correltion_SS$estimate, 2), nsmall = 2),sep = "")
  cor_tag_SC <- paste (paste(SC_tag," Pearson's Correlation: ",sep = ""),
                       format(round(correltion_SC$estimate, 2), nsmall = 2),sep = "")
  
  #plot scatter plot with log2 transformation
  pp <- ggplot(totol_mean_df, aes(x=TotalMean, y= SurfaceMean,color = factor(Tag)))+
    geom_point(alpha = 1/10)+scale_x_continuous(trans="log2")+scale_y_continuous(trans="log2")+
    geom_smooth(method = lm,se = TRUE)+
    stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))+
    annotate("text",x = 0.02*max(totol_mean_df$TotalMean),
             y = 9*min(totol_mean_df$SurfaceMean),color = "#F8766D",label = cor_tag_SC)+
    annotate("text",x = 0.02*max(totol_mean_df$TotalMean),
             y = 1.5*min(totol_mean_df$SurfaceMean),color = "#00BFC4",label = cor_tag_SS)+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 18,face="bold"),
          axis.text.y = element_text(size = 18,face="bold"),
          plot.title = element_text(size=15,face="bold"))+
    ggtitle("Distribution of log2(SurfaceMean) Vs log2(Total Mean) of genes per method")
  pp$labels$colour <- "Gene_type"
  
  return (pp)
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
  
  #####With_filter######
  # dds_SS_keep <-dds_SS[rowSums(counts(dds_SS)) >10,]
  # dds_SS_keep <- DESeq(dds_SS_keep,betaPrior=FALSE)
  # res_dds_SS_keep <- results(dds_SS_keep)
  # 
  # dds_SC_keep <-dds_SC[rowSums(counts(dds_SC)) >10,]
  # dds_SC_keep<-DESeq(dds_SC_keep,betaPrior=FALSE)
  # res_dds_SC_keep<-results(dds_SC_keep)
  # 
  # res_dds_SC_df <- as.data.frame(res_dds_SC_keep) %>% rownames_to_column(.,"Gene_name")
  # res_dds_SS_df <- as.data.frame(res_dds_SS_keep) %>% rownames_to_column(.,"Gene_name")
  # count_merge_10 <- merge(res_dds_SC_df,res_dds_SS_df,by = "Gene_name")
  
  return(res_dds)
}

# Function for getting and annotating the crossedTech genes
# Threshold: P value smaller than 0.05 and log2foldchange larger than 2
# Input: Two DEseq result dfs, SC first and SS (has to be in order)
#        The raw gene info with ensembleID and RNA type
#        threshold for log2FoldChange
#        Theshold for p value
# Output: DataFrame with gene name, ENSEMBL, Official GeneID, log2FoldChange for both methods
get_crossed_geneList <- function (DE_df_1,DE_df_2,FPKM_info,log2Thresh,pvalThresh){
  #find cross tech genes and filter genes by input threshold
  count_merge_NF <- merge(DE_df_1,DE_df_2,by = "Gene_name")
  cross_tech_genes <- subset(count_merge_NF,log2FoldChange.x>log2Thresh&
                               log2FoldChange.y>log2Thresh&padj.x<pvalThresh&padj.y<pvalThresh)
  crossed_gene_info <- merge(cross_tech_genes,FPKM_info, by.x = "Gene_name",
                             by.y = "gene_short_name") %>% .[,c(1,3,9,14,15)]
  
  #annotate cross tech genes
  gene_ID <- crossed_gene_info$tracking_id %>% as.character(.)
  cols <- c("GENENAME")
  anno_info <- select(org.Mm.eg.db, keys=gene_ID, columns=cols, keytype="ENSEMBL")
  
  #format output dataframe
  res <- merge (crossed_gene_info,anno_info, by.x = "tracking_id", by.y = "ENSEMBL")
  res <- res[,c(1,2,5,6,3,4)]
  
  #save files to current directory
  write.csv(res,"crossed_gene_specific_info.csv",row.names=FALSE)
  
  return (res)
}

# Function for plot the global scater plot of log2FoldChange with linear regression and Pearson's correlation
# Input: Deres_SC, DEres_SS (has to be in order)
#        crossed_genes df (could be generated by using function get_crossed_geneList())
#        fish_str: a string of gene name - all the fish genes for detection
#        fish_sig_str: a string of gene name - genes that showed positibe signals for detection
#        The x axis label (input as string)
#        The y axis label (input as string)
# Output: Return the plot object
#         Save a CSV file for count NAs
plot_FoldChange <- function (df_1,df_2,cross_df,fish_str,fish_sig_str,xlab_tag,ylab_tag){
  #merge the dataframe for later use
  count_merge_NF <- merge(df_1,df_2,by = "Gene_name")
  
  #Get the data from the count_merge_NF by using FISH data
  fish_df <- count_merge_NF[count_merge_NF$Gene_name %in% fish_str, ]
  fish_sig_df <- count_merge_NF[count_merge_NF$Gene_name %in% fish_sig_str, ]
  
  #preform Pearson's correlation test
  correlation <- cor.test(count_merge_NF$log2FoldChange.x, count_merge_NF$log2FoldChange.y, method = "pearson", conf.level = 0.95)
  cor_tag <- paste ("Pearson's Correlation: ",format(round(correlation$estimate, 2), nsmall = 2),sep = "")
  
  #plot graphs
  line_types <- c("Linear_regression"=1,"Loess_0.05"=6)
  p <- ggplot(count_merge_NF, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
    geom_hline(aes(yintercept = 0),color = "gray") +
    geom_vline(aes(xintercept = 0),color = "gray")+
    geom_point(aes(color = "Genes"),alpha = 1/20)+
    xlab(xlab_tag)+ylab(ylab_tag)+
    geom_smooth(method = "lm", aes(linetype = "Linear_regression"),color = "cornflowerblue")+
    geom_point(data = cross_df,aes(color = "Genes_Crossed_Techniques"),alpha = 1/4)+
    stat_regline_equation(label.x = 0,label.y = -5,
                          aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))+
    annotate("text",x = 4,y = 12,label = cor_tag)+
    geom_point(data = fish_df,aes(color = "FISH_genes"))+
    geom_point(data = fish_sig_df,aes(color = "FISH_pos_genes"))+
    #annotate("text",x = -5,y = 10,label = cor_tag)
    scale_color_manual(name = "Gene_type",values = c("green","#e59a3d","black","red"))+
    scale_linetype_manual(name = "Regression_method",values = c("solid"))+
    ylim(-10,15)+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 18,face="bold"),
          axis.text.y = element_text(size = 18,face="bold"),
          plot.title = element_text(size=15,face="bold"))+
    ggtitle("All Genes log2FoldChange per Technique with Linear Regression Line")
  
  p <- p + geom_text(data = fish_df,aes(label=Gene_name),color = "#e59a3d",hjust=0, vjust=0)
  #Generating the count table 
  #for both techniques
  v1 <- length(which(is.na(count_merge_NF$log2FoldChange.y)&is.na(count_merge_NF$log2FoldChange.x)))
  #only Surface Click NA
  v2 <- length(which(!is.na(count_merge_NF$log2FoldChange.y)&is.na(count_merge_NF$log2FoldChange.x)))
  #only Surface Seq NA
  v3 <- length(which(is.na(count_merge_NF$log2FoldChange.y)&!is.na(count_merge_NF$log2FoldChange.x)))
  
  tt <- cbind(as.data.frame(c("Both Techs","Surface_Click_only","Surface_seq_only","Total_NA")),
              as.data.frame(c(v1,v2,v3,sum(v1,v2,v3))))
  colnames(tt) <-c("Type","Count")
  
  #output files 
  write.csv(tt,"NA_count_summary.csv",row.names=FALSE)
  return(p)
}

# Function for plot the global scater plot with no intersection of genes
# of log2FoldChange with linear regression and Pearson's correlation
# Input: Deres_SC, DEres_SS (has to be in order)
#        crossed_genes df (could be generated by using function get_crossed_geneList())
#        The x axis label (input as string)
#        The y axis label (input as string)
# Output: Return the plot object
#         Save a CSV file for count NAs
plot_FoldChange_noCross <- function (df_1,df_2,xlab_tag,ylab_tag){
  #merge the dataframe for later use
  count_merge_NF <- merge(df_1,df_2,by = "Gene_name")
  
  #preform Pearson's correlation test
  correlation <- cor.test(count_merge_NF$log2FoldChange.x, count_merge_NF$log2FoldChange.y, method = "pearson", conf.level = 0.95)
  cor_tag <- paste ("Pearson's Correlation: ",format(round(correlation$estimate, 2), nsmall = 2),sep = "")
  
  #plot graphs
  line_types <- c("Linear_regression"=1,"Loess_0.05"=6)
  p <- ggplot(count_merge_NF, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
    geom_hline(aes(yintercept = 0),color = "gray") +
    geom_vline(aes(xintercept = 0),color = "gray")+
    geom_point(aes(color = "Genes"),alpha = 1/20)+
    xlab(xlab_tag)+ylab(ylab_tag)+
    geom_smooth(method = "lm", aes(linetype = "Linear_regression"),color = "cornflowerblue")+
    #geom_point(data = cross_df,aes(color = "Crossed_Technique_genes"),alpha = 1/4)+
    stat_regline_equation(label.x = 0,label.y = -5,
                          aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))+
    annotate("text",x = 4,y = 12,label = cor_tag)+
    #annotate("text",x = -5,y = 10,label = cor_tag)
    scale_color_manual(name = "Gene_type",values = c("black"))+
    scale_linetype_manual(name = "Regression_method",values = c("solid"))+
    ylim(-10,15)+
    ggtitle("All Genes log2FoldChange per Technique with Linear Regression Line")
  
  #Generating the count table 
  #for both techniques
  v1 <- length(which(is.na(count_merge_NF$log2FoldChange.y)&is.na(count_merge_NF$log2FoldChange.x)))
  #only Surface Click NA
  v2 <- length(which(!is.na(count_merge_NF$log2FoldChange.y)&is.na(count_merge_NF$log2FoldChange.x)))
  #only Surface Seq NA
  v3 <- length(which(is.na(count_merge_NF$log2FoldChange.y)&!is.na(count_merge_NF$log2FoldChange.x)))
  
  tt <- cbind(as.data.frame(c("Both Techs","Surface_Click_only","Surface_seq_only","Total_NA")),
              as.data.frame(c(v1,v2,v3,sum(v1,v2,v3))))
  colnames(tt) <-c("Type","Count")
  
  #output files 
  write.csv(tt,"NA_count_summary.csv",row.names=FALSE)
  return(p)
}

# Functions for plottong volcano plot
# input: DEresSC: the DEres df for the first method
#        SC_tag: gene tag for the first method, whill be shown in the legend, color red
#        SS_inSC_tag: gene tag for the second method's genes in first method, whill be shown in the legend, color blue
#        DEres_SS: DEres df for the second method
#        SS_tag: gene tag for the second method, whill be shown in the legend, color blue
#        SC_inSC_tag: gene tag for the first method's genes in second method, whill be shown in the legend, color red      
# output: Three volcano plots as list
#         1: total SC, SS
#         2: SS and SS selected
#         3: SC and SC selected
plot_Volcanos<- function(df_SC,SC_tag,SS_inSC_tag,df_SS,SS_tag,SC_inSS_tag){
  cand_SC <- subset(df_SC,log2FoldChange>2&padj<0.05) %>% as.data.frame() %>% .$Gene_name
  select_df_SS <- df_SS[df_SS$Gene_name %in% cand_SC, ]
  cand_SS <- subset(df_SS,log2FoldChange>2&padj<0.05) %>% as.data.frame() %>% .$Gene_name
  select_df_SC <- df_SC[df_SC$Gene_name %in% cand_SS, ]
  
  #two method total genes
  cols <- c("#00BFC4","#F8766D")
  names(cols) <- c(SS_tag,SC_tag)
  #cols <-c(SS_tag = "#00BFC4",SC_tag = "#F8766D")
  g_total <- ggplot(data=df_SS, aes(x=log2FoldChange, y=padj)) +
    geom_rect(aes(xmin = 2,xmax = Inf,ymin = -Inf,ymax = 0.05,fill = "Significant Area"),alpha = 0.2)+
    geom_point(aes(color = SS_tag),alpha=0.2, size=1.75)+
    geom_point(aes(color = SC_tag),data = df_SC,alpha=0.4)+
    geom_hline(aes(yintercept = 0.05),color = "gray") +
    scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.05)),trans = 'reverse')+
    geom_vline(aes(xintercept = 2),color = "gray")+
    scale_x_continuous(breaks = sort(c(seq(-10, 10, length.out=5), 2)))+
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1))+
    scale_fill_manual('Highlight',values = 'lightblue',
                      guide = guide_legend(override.aes = list(alpha = 1)))+
    scale_color_manual(name = "Gene Type",values = cols)+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 18,face="bold"),
          axis.text.y = element_text(size = 18,face="bold"),
          plot.title = element_text(size=15,face="bold"),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text=element_text(size=12))+
    guides(color = guide_legend(order = 1),
           fill = guide_legend(order = 2))+
    ggtitle("All genes in Two Techniques")
  
  #Tech SS with cands in SC
  cols <- c("#00BFC4","#d14238")
  names(cols) <- c(SS_tag,SC_inSS_tag)
  #cols <-c(SS_tag = "#00BFC4",SC_inSS_tag = "#F8766D")
  #title_tag <- paste("VolcanoPlot: ",SC_inSS_tag,sep = "")
  g_SS <- ggplot(data=df_SS, aes(x=log2FoldChange, y=padj)) +
    geom_rect(aes(xmin = 2,xmax = Inf,ymin = -Inf,ymax = 0.05,fill = "Significant Area"),alpha = 0.2)+
    geom_point(aes(color = SS_tag),alpha=0.2, size=1.75)+
    geom_point(aes(color = SC_inSS_tag),data = select_df_SS,alpha=0.4)+
    geom_hline(aes(yintercept = 0.05),color = "gray") +
    scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.05)),trans = 'reverse')+
    geom_vline(aes(xintercept = 2),color = "gray")+
    scale_x_continuous(breaks = sort(c(seq(-10, 10, length.out=5), 2)))+
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1))+
    scale_color_manual(name = "Gene Type",values = cols)+
    scale_fill_manual('Highlight',values = 'lightblue',
                      guide = guide_legend(override.aes = list(alpha = 1)))+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 18,face="bold"),
          axis.text.y = element_text(size = 18,face="bold"),
          plot.title = element_text(size=15,face="bold"),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text=element_text(size=12))+
    guides(color = guide_legend(order = 1),
           fill = guide_legend(order = 2))+
    ggtitle(SC_inSS_tag)
  
  #Tech SC with cands in SS
  cols <- c("#F8766D","#00979b")
  names(cols)<- c(SC_tag,SS_inSC_tag)
  #cols <-c(SC_tag = "#F8766D",SS_inSC_tag = "#00BFC4")
  #title_tag <- paste("VolcanoPlot: ",SS_inSC_tag,sep = "")
  g_SC <- ggplot(data=df_SC, aes(x=log2FoldChange, y=padj)) +
    geom_rect(aes(xmin = 2,xmax = Inf,ymin = -Inf,ymax = 0.05,fill = "Significant Area"),alpha = 0.2)+
    geom_point(aes(color = SC_tag),alpha=0.2, size=1.75)+
    geom_point(aes(color = SS_inSC_tag),data = select_df_SC,alpha=0.4)+
    geom_hline(aes(yintercept = 0.05),color = "gray") +
    scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.05)),trans = 'reverse')+
    geom_vline(aes(xintercept = 2),color = "gray")+
    scale_x_continuous(breaks = sort(c(seq(-10, 10, length.out=5), 2)))+
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1))+
    scale_color_manual(name = "Gene Type",values = cols)+
    scale_fill_manual('Highlight',values = 'lightblue',
                      guide = guide_legend(override.aes = list(alpha = 1)))+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 18,face="bold"),
          axis.text.y = element_text(size = 18,face="bold"),
          plot.title = element_text(size=15,face="bold"),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text=element_text(size=12))+
    guides(color = guide_legend(order = 1),
           fill = guide_legend(order = 2))+
    ggtitle(SS_inSC_tag)
  
  p_com <- grid.arrange(g_total, g_SS, g_SC,nrow = 1,
                        top = textGrob("Volcano Plots", gp=gpar(fontsize=20,fontface="bold")))
  return (p_com)
  # g
  # ggsave("VolcanoPlot_SCgenes_in_SC.pdf")
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
