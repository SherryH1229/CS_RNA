library(DESeq2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pasilla)
library(ggplot2)
library(dplyr)
library(topGO)
library(tidyverse)

#----------------Prepare Original SurfaceClick data-----------------------
setwd("~/dataOS/CS_RNA/surface_click/featureCount/")
RNAclick_counts <- read.csv("count_matrix.txt",sep = "\t", head = TRUE,skip=1,row.names = "Geneid")

#format the data, remove unneccessary cells and trnasform the colnames
RNAclick_counts <- RNAclick_counts[ ,6:ncol(RNAclick_counts)]
sampleName <- str_split_fixed(colnames(RNAclick_counts),"\\.",7) 
sampleName <- gsub(".bam","",sampleName)
sampleName <- gsub("datapool.","",sampleName)[,7]
colnames(RNAclick_counts) <- sampleName
#------------Perpare Original SurfaceSeq data -----------------
setwd("~/dataOS/CS_RNA/surface_seq/featureCount")
RNAseq_counts <- read.csv("count_matrix.txt",sep = "\t", head = TRUE,skip=1,row.names = "Geneid")

#format the data, remove unneccessary cells and trnasform the colnames
RNAseq_counts <- RNAseq_counts[ ,6:ncol(RNAseq_counts)]
sampleName <- str_split_fixed(colnames(RNAseq_counts),"\\.",7) 
sampleName <- gsub(".accepted_hits.bam","",sampleName)
sampleName <- gsub("datapool.","",sampleName)[,7]

colnames(RNAseq_counts) <- sampleName
RNAseq_counts$NP.ER.Son.1 <- RNAseq_counts$NP.ER.Son.1+RNAseq_counts$NP.ER.Son.2
RNAseq_counts$NP.ER.Son.2 <- NULL
colnames(RNAseq_counts) <- c("NP.ER.Son","NP.SAL.Son","NP.SAL.Ext","Tot","Tot.ER","IC.SAL.1","IC.SAL.2")
#--------------Plot meanSurface Vs meanTotal from original data (non filtering) ---------------
#calculated Mean surface and mean total for SC and SS
#Surface click Data
FC_df_SC <- as.data.frame(rowMeans(RNAclick_counts[,1:2]+1)) %>% cbind(.,rowMeans(RNAclick_counts[,4:5])+1)
FC_df_SC$Tag <- "SC_genes"
colnames(FC_df_SC) <- c("SurfaceMean","TotalMean","Tag")
#surface seq data
FC_df_SS <- as.data.frame(rowMeans(RNAseq_counts[,c(1,2,3,6,7)]+1)) %>% cbind(.,rowMeans(RNAseq_counts[,4:5])+1)
FC_df_SS$Tag <- "SS_gene"
colnames(FC_df_SS) <- c("SurfaceMean","TotalMean","Tag")

#combine data together
totol_mean_methodOnly_df <- rbind(FC_df_SC,FC_df_SS)
cor.test(FC_df_SS$SurfaceMean, FC_df_SS$TotalMean, method = "pearson", conf.level = 0.95)

pp <- ggplot(totol_mean_methodOnly_df, aes(x=TotalMean, y= SurfaceMean,color = factor(Tag)))+
  geom_point(alpha = 1/10)+scale_x_continuous(trans="log2")+scale_y_continuous(trans="log2")+
  geom_smooth(method = lm,se = TRUE)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))+
  annotate("text",x = 17,y = 17,color = "red4",label = "SC Pearson's Correlation: 0.9729862\nP-value < 2.2e-16")+
  annotate("text",x = 1024,y = 1024,color = "blue4",label = "SS Pearson's Correlation: 0.4374715\nP-value < 2.2e-16")+
  ggtitle("Distribution of log2(SurfaceMean) Vs log2(Total Mean) of genes per method")
pp$labels$colour <- "Gene_type"
pp
ggsave("Gene_Method_Dist_log2.pdf")

#-----------Deseq-------------
condition <- factor(c(rep("Surface",3), rep("Total", 2),rep("Surface",2)),levels = c("Total","Surface"))
col_data_SS <- data.frame(row.names = colnames(as.matrix(RNAseq_counts)),condition)
dds_SS <- DESeqDataSetFromMatrix(countData = as.matrix(RNAseq_counts),colData = col_data_SS,design=~condition)

condition <- factor(c(rep("Surface",2), rep("Total", 3)),levels = c("Total","Surface"))
col_data_SC <- data.frame(row.names = colnames(as.matrix(RNAclick_counts)),condition)
dds_SC <- DESeqDataSetFromMatrix(countData = as.matrix(RNAclick_counts),colData = col_data_SC,design=~condition)

###---------------Non_filter--------------
dds_SS <- DESeq(dds_SS,betaPrior=FALSE)
res_dds_SS <- results(dds_SS)
dds_SC<-DESeq(dds_SC,betaPrior=FALSE)
res_dds_SC<- results(dds_SC)

res_dds_SC_df <- as.data.frame(res_dds_SC) %>% rownames_to_column(.,"Gene_name")
res_dds_SS_df <- as.data.frame(res_dds_SS) %>% rownames_to_column(.,"Gene_name")

##Get cross tech gens with filter of log2FoldChange > 2 and padj < 0.05
count_merge_NF <- merge(res_dds_SC_df,res_dds_SS_df,by = "Gene_name")
cross_tech_genes <- subset(count_merge_NF,log2FoldChange.x>2&
                             log2FoldChange.y>2&padj.x<0.05&padj.y<0.05)
#get candidate gene info
setwd("~/dataOS/CS_RNA/surface_seq/FPKM_plot")
sample_info_SS <- read.csv("FPKM_merged_selected_20181008.txt",sep="\t") 

crossed_gene_info <- merge(cross_tech_genes,sample_info_SS, by.x = "Gene_name",by.y = "gene_short_name") %>% 
  .[,c(1,3,9,14,15)]
setwd("~/dataOS/CS_RNA/crossed_validation")
write.csv(crossed_gene_info,"crossed_gene_info.csv")
DAVID_info_table <- read.csv("output_info.csv",sep = "|",head = TRUE)

merge_info <- merge(crossed_gene_info,DAVID_info_table,
                    by.x = "Gene_name",by.y = "OFFICIAL_GENE_SYMBOL")
write.csv(merge_info,"crossed_gene_specific_info.csv")

#####With_filter#####
dds_SS_keep <-dds_SS[rowSums(counts(dds_SS)) >10,]
dds_SS_keep <- DESeq(dds_SS_keep,betaPrior=FALSE)
res_dds_SS_keep <- results(dds_SS_keep)

dds_SC_keep <-dds_SC[rowSums(counts(dds_SC)) >10,]
dds_SC_keep<-DESeq(dds_SC_keep,betaPrior=FALSE)
res_dds_SC_keep<-results(dds_SC_keep)

res_dds_SC_df <- as.data.frame(res_dds_SC_keep) %>% rownames_to_column(.,"Gene_name")
res_dds_SS_df <- as.data.frame(res_dds_SS_keep) %>% rownames_to_column(.,"Gene_name")
count_merge_10 <- merge(res_dds_SC_df,res_dds_SS_df,by = "Gene_name")

#Count the number of NA
length(which(is.na(test4$log2FoldChange.y)&is.na(test4$log2FoldChange.x)))
length(which(!is.na(test4$log2FoldChange.y)&is.na(test4$log2FoldChange.x)))
rownames(test4,which(is.na(test4$log2FoldChange.y)|is.na(test4$log2FoldChange.x)))

#####scatter plot of log2FoldChange#####
#Pierseen correlation
cor.test(count_merge_NF$log2FoldChange.x, count_merge_NF$log2FoldChange.y, method = "pearson", conf.level = 0.95)
linearMod <- lm(log2FoldChange.y ~ log2FoldChange.x,data = count_merge_NF)
summary(linearMod)

line_types <- c("Linear_regression"=1,"Loess_0.05"=6)
p <- ggplot(count_merge_NF, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_hline(aes(yintercept = 0),color = "gray") +
  geom_vline(aes(xintercept = 0),color = "gray")+
  geom_point(aes(color = "Genes"),alpha = 1/20)+
  xlab("log2FoldChange_SurfaceClick")+ylab("log2FoldChange_SurfaceSeq")+
  geom_smooth(method = "lm", aes(linetype = "Linear_regression"),color = "cornflowerblue")+
  geom_point(data = cross_tech_genes,aes(color = "Corssed_Technique_genes"),alpha = 1/4)+
  stat_regline_equation(label.x = -7,label.y = 3,
                        aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))+
  annotate("text",x = 4,y = 12,label = "Pearson's Correlation: 0.4370606\nP-value < 2.2e-16")+
  scale_color_manual(name = "Gene_type",values = c("red","black"))+
  scale_linetype_manual(name = "Regression_method",values = c("solid"))+ylim(-10,15)+
  ggtitle("All Genes log2FoldChange per Technique with Linear Regression Line")
p
ggsave("Nonfil_allGene_foldChange_NF.pdf")

#----------------Candidate genes--------------
#########Surface seq############
setwd("~/dataOS/CS_RNA/surface_seq/FPKM_plot")
sample_info_SS <- read.csv("FPKM_merged_selected_20181008.txt",sep="\t") 

crossed_gene_info <- merge(cross_tech_genes,sample_info_SS, by.x = "Gene_name",by.y = "gene_short_name") %>% 
  .[,c(1,3,9,14,15)]
setwd("~/dataOS/CS_RNA/crossed_validation")
write.csv(crossed_gene_info,"crossed_gene_info.csv")
DAVID_info_table <- read.csv("output_info.csv",sep = "|",head = TRUE)

merge_info <- merge(crossed_gene_info,DAVID_info_table,
                    by.x = "Gene_name",by.y = "OFFICIAL_GENE_SYMBOL")
write.csv(merge_info,"crossed_gene_specific_info.csv")
##########Surface click#############
#sample_info_SC <- read.csv("/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt",sep="\t") 
#Gene_type_info_SC <- as.data.frame(sample_info_SC[,c(1,2,3)])
all_gene_matrix_SC <- as.data.frame(res_dds_SC_keep) %>% rownames_to_column(.,"gene_Name")
#%>% .[,c("log2FoldChange","padj"),drop=FALSE] %>% 
#  rownames_to_column(.,"gene_Name") %>% merge(Gene_type_info_SC,.,by.x="gene_short_name",by.y = "gene_Name") %>% 
#  .[order(-.$log2FoldChange),]


cand_list_SC <- (subset(all_gene_matrix_SC,log2FoldChange>=2 &!(is.na(padj)) & padj < 0.05))$gene_Name %>% 
  as.character(.)
candiates_SC <- all_gene_matrix_SC[all_gene_matrix_SC$gene_Name %in% cand_list_SC, ] %>% .[order(-.$log2FoldChange),]



ttt <- merge (candiates_SC,candiates_SS,by = "gene_Name")

p <- ggplot(ttt, aes(x=log2FoldChange.x, y=log2FoldChange.y)) + geom_point()+
  xlab("padj_SurfaceSeq")+ylab("padj_SurfaceClick")+
  ggtitle("All Genes padj per Technique")
p
#-------------Gene count table---------
RNA_annotation <- as.character(unique(sample_info_SC$Biotype))
RNA_mainType <- grep("RNA",RNA_annotation,value=TRUE) %>% c(.,"protein_coding")
RNA_otherType <- RNA_annotation[!(RNA_annotation %in% RNA_mainType)]
RNA_main_df <- as.data.frame(RNA_mainType)
RNA_other_df <- as.data.frame(RNA_otherType)

SC_df <- counts_matrix_SC %>% as.data.frame(.)
SC_df$gene_short_name <-rownames(SC_df)
SC_df <- merge(SC_df,Gene_type_info_SC,by = "gene_short_name")
SC_geneCount_1000<- get_gene_count(1000,SC_df,2,6,8)
write.csv(SC_geneCount_1000,"SC_geneCount_1000.csv")

SS_df <- counts_matrix_SS %>% as.data.frame(.)
SS_df$gene_short_name <-rownames(SS_df)
SS_df <- merge(SS_df,Gene_type_info_SS,by = "gene_short_name")
SS_geneCount_1000 <- get_gene_count(1000,SS_df,2,8,10) 
write.csv(SS_geneCount_1000,"SS_geneCount_1000.csv")

#all genes
test <- filter_all(counts_matrix_SC, all_vars(. > 10))
test <- counts_matrix_SC[apply(counts_matrix_SC[,c(1,2,3)], MARGIN = 1, function(x) all(x > 100)), ] %>% as.data.frame(.)

gene_list_100 <- as.vector(rownames(test))
test_res <- gene_list_100[gene_list_100 %in% cross_candidate_list]
#-------------Mereged candidate-------------
#crossed candiated
candidate_merge <- merge(candiates_SC,candiates_SS,by = "tracking_id")
write.csv(candidate_merge,"cross_candidate_genes_info.csv")

#Candidates in SS not in SC, and merge back to SC
cross_candidate_list <- as.vector(candidate_merge$gene_short_name.x)
candiate_inSS_notSC <-cand_list_SS[!(cand_list_SS %in% cross_candidate_list)]
df_SC <- all_gene_matrix_SC[all_gene_matrix_SC$gene_short_name %in% candiate_inSS_notSC, ]
df_SS <- all_gene_matrix_SS[all_gene_matrix_SS$gene_short_name %in% candiate_inSS_notSC, ]
df_SC_SS <- merge(df_SC,df_SS,by = "tracking_id")
write.csv(df_SC_SS,"candidate_inSS_notSC_genes_info.cvs")

#all genes
L10_counts_matrix_SC<- counts_matrix_SC[apply(counts_matrix_SC, MARGIN = 1, function(x) all(x > 10)), ] %>% as.data.frame(.)
L10_counts_matrix_SS<- counts_matrix_SS[apply(counts_matrix_SS, MARGIN = 1, function(x) all(x > 10)), ] %>% as.data.frame(.)

RM_SC <-data.frame(as.vector(rownames(L10_counts_matrix_SC)),as.vector(rowMeans(L10_counts_matrix_SC)))
colnames(RM_SC) <- c("gene_name","mean_count")
RM_SS <-data.frame(as.vector(rownames(L10_counts_matrix_SS)),as.vector(rowMeans(L10_counts_matrix_SS)))
colnames(RM_SS) <- c("gene_name","mean_count")

RM_total <- merge(RM_SC,RM_SS,by = "gene_name")
#---------Merge all genes----------------

test <- merge (all_gene_matrix_SC,all_gene_matrix_SS,by = "tracking_id")

p <- ggplot(test, aes(x=log2FoldChange.x, y=log2FoldChange.y)) + geom_point()+
  xlab("log2FoldChange_SurfaceClick")+ylab("log2FoldChange_SurfaceSeq")
p
p <- ggplot(test, aes(x=padj.x, y=padj.y)) + geom_point()+
  xlab("padj_SurfaceClick")+ylab("padj_SurfaceSeq")+ggtitle()
p
#----------plot scatter plots------------------
setwd("~/dataOS/CS_RNA/")
p <- ggplot(candidate_merge, aes(x=log2FoldChange.x, y=log2FoldChange.y)) + geom_point()+
  xlab("log2FoldChange_SurfaceClick")+ylab("log2FoldChange_SurfaceSeq")
p + geom_point(data=df_SC_SS,colour="red")
ggsave("Scatter_plot_log2FoldChange.pdf")

p <- ggplot(candidate_merge, aes(x=padj.x, y=padj.y)) + geom_point()+
  xlab("padj_SurfaceClick")+ylab("padj_SurfaceSeq")
p + geom_point(data=df_SC_SS,colour="red")
ggsave("Scatter_plot_padj.pdf")

p <- ggplot(RM_total, aes(x=mean_count.x, y=mean_count.y)) + geom_point()+
  xlab("mean_count_SurfaceClick")+ylab("mean_count_SurfaceSeq")+scale_x_log10()+scale_y_log10()
p
ggsave("Scatter_plot_mean_count.pdf")

t1 <- RM_total[,c(1,2)]
t1$method <- "Surface_click"
colnames(t1)[2]<-"mean_count" 

t2 <- RM_total[,c(1,3)]
t2$method <- "Surface_seq"
colnames(t2)[2]<-"mean_count" 

t3 <- rbind(t1,t2)
ggplot(t3, aes(x=mean_count, color=method)) +geom_density()+xlim(0, 10000)
ggsave("mean_count_distribution.pdf")
#---------DAVID----------
DAVID_info <- read.csv("DAVID_indo_procc.csv",sep = "\t",row.names=NULL)
colnames(DAVID_info)[1] <- c("tracking_id")
all_info <- merge(DAVID_info,candidate_merge,by  = "tracking_id")
write.csv(all_info,"call_info.csv")

#------------Functions------------
count_gene_num <- function(thresh,colNum,info_col,ori_df){
  sub_df <- ori_df[c(colNum,info_col)]
  sample_name <- colnames(sub_df)[1]
  sub_df_keep <- sub_df[which(sub_df[1] > thresh),]
  gene_count_table <- table(sub_df_keep[2]) %>% as.data.frame(.)
  gene_count_df <- merge(RNA_main_df,gene_count_table,by.x = "RNA_mainType", by.y  = "Var1")
  gene_count_other <- merge(RNA_other_df,gene_count_table,by.x = "RNA_otherType", by.y  = "Var1")
  colnames(gene_count_df) <- c("RNA_type",sample_name)
  gene_count_other_list <- c("Other",sum(gene_count_other[2])) %>% t(.) %>% as.data.frame(.)
  colnames(gene_count_other_list) <- c("RNA_type",sample_name)
  gene_count_df_res <- rbind(gene_count_df,gene_count_other_list)
  rownames(gene_count_df_res) <- gene_count_df_res$RNA_type
  gene_count_df_res$RNA_type <- NULL
  return (gene_count_df_res)
}

get_gene_count <- function (threshold_num,raw_df,start_pos,end_pos,info_col){
  for (index in start_pos:end_pos){
    temp_df <- count_gene_num(threshold_num,index,info_col,raw_df)
    if (!exists("count_total_df")){
      count_total_df <- temp_df
    }
    else{
      count_total_df <- cbind(count_total_df,temp_df)
    }
  }
  return (count_total_df)
}

getCands <- function(res_dds){
  temp_df <- as.data.frame(res_dds) %>% rownames_to_column(.,"Gene_name")
  temp_cand_list <- (subset(temp_df,log2FoldChange>=2&!(is.na(padj)) & padj < 0.05))$Gene_name %>% as.character(.)
  temp_cand_select <- temp_df[temp_df$Gene_name %in% temp_cand_list, ] 
  
  return (temp_cand_select)
}