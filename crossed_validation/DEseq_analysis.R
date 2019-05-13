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

#----------------Data Preperaion-----------------------
#FPKM, ensembleID,and biotype sample info for all samples
FPKM_info_all <- read.csv("/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt",
                           sep="\t") %>% as.data.frame(.)
############## Surface Clieck data
setwd("~/dataOS/CS_RNA/surface_click/featureCount/")
RNAclick_counts <- read.csv("count_matrix.txt",sep = "\t", head = TRUE,skip=1,row.names = "Geneid")

#format the data, remove unneccessary cells and trnasform the colnames
RNAclick_counts <- RNAclick_counts[ ,6:ncol(RNAclick_counts)]
sampleName <- str_split_fixed(colnames(RNAclick_counts),"\\.",7) 
sampleName <- gsub(".bam","",sampleName)
sampleName <- gsub("datapool.","",sampleName)[,7]
colnames(RNAclick_counts) <- sampleName

#Get FPKM info
FPKM_info_SC <- FPKM_info_all[,c(1,2,3,14:16,19,20)]
############### SurfaceSeq data
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

#FPKM, ensembleID,and biotype sample info 
FPKM_info_SS <- read.csv("~/dataOS/CS_RNA/surface_seq/FPKM_plot/FPKM_merged_selected_20181008.txt",
                           sep="\t") %>% as.data.frame(.) %>% .[,c(1:5,7:11)]

#-------------------FPKM_bar_plot:Diff FPKM Threshd: 0, 5, 10----------
#remove unqualified data
FPKM_selected_info <- FPKM_info_all[,-which(
  names(FPKM_info_all) %in% c("EL4.Mem.ER","EL4.NP.ER.Ext","EL4.Click.3","EL4.Click.woCu"))] %>% .[!duplicated(.$tracking_id),]

#Get type info

RNA_annotation <- unique(FPKM_selected_info$Biotype) %>% as.character(.)
RNA_mainType <- grep("RNA",RNA_annotation,value=TRUE) %>% c(.,"protein_coding")
#RNA_mainType
#Get Count of FPKM with three threshold 0,5,10, and save as csv file
#Thresh = 0
setwd("~/dataOS/CS_RNA/crossed_validation")
FPKM_0_SS <- geneCount_FPKM(FPKM_selected_info,4,11,0)
FPKM_0_SC <- geneCount_FPKM(FPKM_selected_info,12,16,0)
FPKM_0_selected_SS <- as.data.frame(t(geneCount_FPKM_select(FPKM_0_SS)))%>% rownames_to_column(.,"RNA_type")
FPKM_0_selected_SC <- as.data.frame(t(geneCount_FPKM_select(FPKM_0_SC)))%>% rownames_to_column(.,"RNA_type")
write.csv(FPKM_0_selected_SS,"SS_FPKM_0_GeneCount.csv")
write.csv(FPKM_0_selected_SC,"SC_FPKM_0_GeneCount.csv")
FPKM_0_selected_SS$Condition <- "FPKM_0"
FPKM_0_selected_SC$Condition <- "FPKM_0"

#Thresh  = 5
FPKM_5_SS <- geneCount_FPKM(FPKM_selected_info,4,11,5)
FPKM_5_SC <- geneCount_FPKM(FPKM_selected_info,12,16,5)
FPKM_5_selected_SS <- as.data.frame(t(geneCount_FPKM_select(FPKM_5_SS)))%>% rownames_to_column(.,"RNA_type")
FPKM_5_selected_SC <- as.data.frame(t(geneCount_FPKM_select(FPKM_5_SC)))%>% rownames_to_column(.,"RNA_type")
write.csv(FPKM_5_selected_SS,"SS_FPKM_5_GeneCount.csv")
write.csv(FPKM_5_selected_SC,"SC_FPKM_5_GeneCount.csv")
FPKM_5_selected_SS$Condition <- "FPKM_5"
FPKM_5_selected_SC$Condition <- "FPKM_5"

#Thresh = 10
FPKM_10_SS <- geneCount_FPKM(FPKM_selected_info,4,11,10)
FPKM_10_SC <- geneCount_FPKM(FPKM_selected_info,12,16,10)
FPKM_10_selected_SS <- as.data.frame(t(geneCount_FPKM_select(FPKM_10_SS)))%>% rownames_to_column(.,"RNA_type")
FPKM_10_selected_SC <- as.data.frame(t(geneCount_FPKM_select(FPKM_10_SC)))%>% rownames_to_column(.,"RNA_type")
write.csv(FPKM_10_selected_SS,"SS_FPKM_10_GeneCount.csv")
write.csv(FPKM_10_selected_SC,"SC_FPKM_10_GeneCount.csv")
FPKM_10_selected_SS$Condition <- "FPKM_10"
FPKM_10_selected_SC$Condition <- "FPKM_10"

#Format plotting data
FPKM_all_SS<- rbind(FPKM_0_selected_SS,FPKM_5_selected_SS,FPKM_10_selected_SS)
FPKM_all_SC<- rbind(FPKM_0_selected_SC,FPKM_5_selected_SC,FPKM_10_selected_SC)

FPKM_all_melt_SS <- melt(FPKM_all_SS,id.vars = c("RNA_type","Condition"))
FPKM_all_melt_SS$Condition <- factor(FPKM_all_melt_SS$Condition,levels = c("FPKM_0","FPKM_5","FPKM_10"))
FPKM_all_melt_SS$RNA_type <- factor(FPKM_all_melt_SS$RNA_type,levels = FPKM_0_selected_SS$RNA_type)

#plot distribution of reads
P <- ggplot(FPKM_all_melt_SS, aes(x=variable, y= value, fill = RNA_type)) +
  geom_bar(stat = "identity") + coord_flip()+facet_grid(Condition ~.) +
  theme(axis.text.x = element_text(angle=-40, hjust=.1)) +
  scale_fill_discrete(name = "Assigned_status") + ylab("Count")+ xlab("Sample_ID")+
  scale_fill_brewer(palette="Set3")+ggtitle("FeatureCount_Assign_Distribution")
P
#--------------Scatter Plot --- meanSurface Vs meanTotal ---------------
#calculated Mean surface and mean total for SC and SS
#Surface click Data
FC_df_SC <- as.data.frame(rowMeans(RNAclick_counts[,1:2]+1)) %>% cbind(.,rowMeans(RNAclick_counts[,4:5])+1)
FC_df_SC$Tag <- "SC_genes"
colnames(FC_df_SC) <- c("SurfaceMean","TotalMean","Tag")
#surface seq data
FC_df_SS <- as.data.frame(rowMeans(RNAseq_counts[,c(1,2,3,6,7)]+1)) %>% cbind(.,rowMeans(RNAseq_counts[,4:5])+1)
FC_df_SS$Tag <- "SS_gene"
colnames(FC_df_SS) <- c("SurfaceMean","TotalMean","Tag")

#combination of data
totol_mean_methodOnly_df <- rbind(FC_df_SC,FC_df_SS)
#calculate peasron correlation
correltion_SS <- cor.test(FC_df_SS$SurfaceMean, FC_df_SS$TotalMean, method = "pearson", conf.level = 0.95)
correltion_SC <- cor.test(FC_df_SC$SurfaceMean, FC_df_SC$TotalMean, method = "pearson", conf.level = 0.95)

cor_tag_SS <- paste ("SS Pearson's Correlation: ",format(round(correltion_SS$estimate, 2), nsmall = 2),
                     "\nP-value < ",correltion_SS$p.value,sep = "")
cor_tag_SC <- paste ("SC Pearson's Correlation: ",format(round(correltion_SC$estimate, 2), nsmall = 2),
                     "\nP-value < ",correltion_SC$p.value,sep = "")

#plot scatter plot with log2 transformation
pp <- ggplot(totol_mean_methodOnly_df, aes(x=TotalMean, y= SurfaceMean,color = factor(Tag)))+
  geom_point(alpha = 1/10)+scale_x_continuous(trans="log2")+scale_y_continuous(trans="log2")+
  geom_smooth(method = lm,se = TRUE)+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))+
  annotate("text",x = 0.05*max(totol_mean_methodOnly_df$TotalMean),
           y = 9*min(totol_mean_methodOnly_df$SurfaceMean),color = "#F8766D",label = cor_tag_SC)+
  annotate("text",x = 0.05*max(totol_mean_methodOnly_df$TotalMean),
           y = 1.5*min(totol_mean_methodOnly_df$SurfaceMean),color = "#00BFC4",label = cor_tag_SS)+
  ggtitle("Distribution of log2(SurfaceMean) Vs log2(Total Mean) of genes per method")
pp$labels$colour <- "Gene_type"
pp
#ggsave("Gene_Methods_Dist_log2.pdf")
#-----------Preform DEseq & Get Crossed Tech Genes -------------
#data structure, factor arrangement
condition <- factor(c(rep("Surface",3), rep("Total", 2),rep("Surface",2)),levels = c("Total","Surface"))
col_data_SS <- data.frame(row.names = colnames(as.matrix(RNAseq_counts)),condition)
dds_SS <- DESeqDataSetFromMatrix(countData = as.matrix(RNAseq_counts),colData = col_data_SS,design=~condition)

condition <- factor(c(rep("Surface",2), rep("Total", 3)),levels = c("Total","Surface"))
col_data_SC <- data.frame(row.names = colnames(as.matrix(RNAclick_counts)),condition)
dds_SC <- DESeqDataSetFromMatrix(countData = as.matrix(RNAclick_counts),colData = col_data_SC,design=~condition)

#perform DEseq without pre-filtering
dds_SS <- DESeq(dds_SS,betaPrior=FALSE)
res_dds_SS <- results(dds_SS)
dds_SC<-DESeq(dds_SC,betaPrior=FALSE)
res_dds_SC<- results(dds_SC)

#####With_filter
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

res_dds_SC_df <- as.data.frame(res_dds_SC) %>% rownames_to_column(.,"Gene_name")
res_dds_SS_df <- as.data.frame(res_dds_SS) %>% rownames_to_column(.,"Gene_name")

#Get cross tech gens with filter of log2FoldChange > 2 and padj < 0.05
count_merge_NF <- merge(res_dds_SC_df,res_dds_SS_df,by = "Gene_name")
cross_tech_genes <- subset(count_merge_NF,log2FoldChange.x>2&
                             log2FoldChange.y>2&padj.x<0.05&padj.y<0.05)

crossed_gene_info <- merge(cross_tech_genes,FPKM_info_SS, by.x = "Gene_name",by.y = "gene_short_name") %>% 
  .[,c(1,3,9,14,15)]
#setwd("~/dataOS/CS_RNA/crossed_validation")
write.csv(crossed_gene_info,"crossed_gene_info.csv")
DAVID_info_table <- read.csv("output_info.csv",sep = "|",head = TRUE)
merge_info <- merge(crossed_gene_info,DAVID_info_table,
                    by.x = "Gene_name",by.y = "OFFICIAL_GENE_SYMBOL")
write.csv(merge_info,"crossed_gene_specific_info.csv")
#-------------Scatter plot: log2FoldChange & Count NA --------------
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

###Count the number of NA for each technique
#both techniques
length(which(is.na(count_merge_NF$log2FoldChange.y)&is.na(count_merge_NF$log2FoldChange.x)))
#only Surface Click NA
length(which(!is.na(count_merge_NF$log2FoldChange.y)&is.na(count_merge_NF$log2FoldChange.x)))
#only Surface Seq NA
length(which(is.na(count_merge_NF$log2FoldChange.y)|is.na(count_merge_NF$log2FoldChange.x)))

#-------------Volcano plot--------------------
#get significant genes in each method: Threshold is log2FoldChange larger than 2 and padj smaller than 0,05
cand_SC <- subset(res_dds_SC,log2FoldChange>2&padj<0.05) %>% as.data.frame()
cand_SS <- subset(res_dds_SS,log2FoldChange>2&padj<0.05) %>% as.data.frame()

#two method total genes
cols <-c("Surface_seq_genes" = "#00BFC4","Surface_click_genes" = "#F8766D")
g <- ggplot(data=res_dds_SS_df, aes(x=log2FoldChange, y=padj)) +
  geom_rect(aes(xmin = 2,xmax = Inf,ymin = -Inf,ymax = 0.05,fill = "Significant Area"),alpha = 0.2)+
  geom_point(aes(color = "Surface_seq_genes"),alpha=0.2, size=1.75)+
  geom_point(aes(color = "Surface_click_genes"),data = res_dds_SC_df,alpha=0.4)+
  geom_hline(aes(yintercept = 0.05),color = "gray") +
  scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.05)),trans = 'reverse')+
  geom_vline(aes(xintercept = 2),color = "gray")+
  scale_x_continuous(breaks = sort(c(seq(-10, 10, length.out=5), 2)))+
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1))+
  scale_fill_manual('Highlight',values = 'lightblue', 
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  scale_color_manual(name = "Gene Type",values = cols)+
  ggtitle("VolcanoPlot: All genes in Two Techniques")
g
ggsave("VolcanoPlot_allGene_twoTeches.pdf")

#Tech SS with cands in SC
cols <-c("Surface_seq_genes" = "#00BFC4","Candidate_SC_genes" = "#F8766D")
g <- ggplot(data=res_dds_SS_df, aes(x=log2FoldChange, y=padj)) +
  geom_rect(aes(xmin = 2,xmax = Inf,ymin = -Inf,ymax = 0.05,fill = "Significant Genes"),alpha = 0.2)+
  geom_point(aes(color = "Surface_seq_genes"),alpha=0.2, size=1.75)+
  geom_point(aes(color = "Candidate_SC_genes"),data = cand_SC,alpha=0.4)+
  geom_hline(aes(yintercept = 0.05),color = "gray") +
  scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.05)),trans = 'reverse')+
  geom_vline(aes(xintercept = 2),color = "gray")+
  scale_x_continuous(breaks = sort(c(seq(-10, 10, length.out=5), 2)))+
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 1))+
  scale_color_manual(name = "Gene Type",values = cols)+
  scale_fill_manual('Highlight',values = 'lightblue', 
                    guide = guide_legend(override.aes = list(alpha = 1))) 
g
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
#------------Functions------------
count_gene_num <- function(thresh,colNum,info_col,ori_df){
#   sub_df <- ori_df[c(colNum,info_col)]
#   sample_name <- colnames(sub_df)[1]
#   sub_df_keep <- sub_df[which(sub_df[1] > thresh),]
#   gene_count_table <- table(sub_df_keep[2]) %>% as.data.frame(.)
#   gene_count_df <- merge(RNA_main_df,gene_count_table,by.x = "RNA_mainType", by.y  = "Var1")
#   gene_count_other <- merge(RNA_other_df,gene_count_table,by.x = "RNA_otherType", by.y  = "Var1")
#   colnames(gene_count_df) <- c("RNA_type",sample_name)
#   gene_count_other_list <- c("Other",sum(gene_count_other[2])) %>% t(.) %>% as.data.frame(.)
#   colnames(gene_count_other_list) <- c("RNA_type",sample_name)
#   gene_count_df_res <- rbind(gene_count_df,gene_count_other_list)
#   rownames(gene_count_df_res) <- gene_count_df_res$RNA_type
#   gene_count_df_res$RNA_type <- NULL
#   return (gene_count_df_res)
# }
# 
# get_gene_count <- function (threshold_num,raw_df,start_pos,end_pos,info_col){
#   for (index in start_pos:end_pos){
#     temp_df <- count_gene_num(threshold_num,index,info_col,raw_df)
#     if (!exists("count_total_df")){
#       count_total_df <- temp_df
#     }
#     else{
#       count_total_df <- cbind(count_total_df,temp_df)
#     }
#   }
#   return (count_total_df)
# }
# 
# getCands <- function(res_dds){
#   temp_df <- as.data.frame(res_dds) %>% rownames_to_column(.,"Gene_name")
#   temp_cand_list <- (subset(temp_df,log2FoldChange>=2&!(is.na(padj)) & padj < 0.05))$Gene_name %>% as.character(.)
#   temp_cand_select <- temp_df[temp_df$Gene_name %in% temp_cand_list, ] 
#   
#   return (temp_cand_select)
}

geneCount_FPKM <- function(sample_info_table,sCol,eCol,Thresh){
  for (type in RNA_annotation){
    info <- as.data.frame(colSums(sample_info_table[sample_info_table$Biotype==type,sCol:eCol]>Thresh))
    colnames(info)<-type
    if (!exists("FPKM_thresh")){
      FPKM_thresh <- info
    }
    else{
      FPKM_thresh <-cbind(FPKM_thresh,info)
    }
  }
  FPKM_thresh<-as.data.frame(t(FPKM_thresh))
  FPKM_thresh<-as.data.frame(t(FPKM_thresh))
  return(FPKM_thresh)
}

geneCount_FPKM_select<-function(FPKM_dataframe){
  for (i in names(FPKM_dataframe)){
    if (i %in% RNA_mainType){
      select <- FPKM_dataframe[,i,drop=FALSE]
      if (!exists("FPKM_select")){
        FPKM_select <- select
      }
      else{
        FPKM_select <-cbind(FPKM_select,select)
      }
    }
    else{
      other <- FPKM_dataframe[,i,drop=FALSE]
      #print (other)
      if (!exists("other_sum")){
        other_sum <-other
      }
      else{
        other_sum <- other_sum+other
        colnames(other_sum) <- "other"
      }
    }
  }
  #protein_coding <- FPKM_dataframe[,"protein_coding",drop=FALSE]
  #FPKM_select <-cbind(FPKM_select,protein_coding)
  FPKM_select <-cbind(FPKM_select,other_sum)
  #FPKM_select$sum <-colSums(FPKM_select)
  FPKM_select$miRNA<-NULL
  return (FPKM_select)
}

plot_Barplot <- function(dataFrame,axis_x,axis_y){
  dataFrame_trans <- melt(dataFrame,axis_x)
  dataFrame_trans$value <- as.numeric(gsub(",","",dataFrame_trans$value))
  colnames(dataFrame_trans)[2] <- axis_y
  
  dataFrame_barplot <- ggplot(dataFrame_trans,aes(x = dataFrame_trans[,axis_x], y = value)) + 
    geom_bar(aes(fill = dataFrame_trans[,axis_y]),stat = "identity",position = "dodge") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle=-40, hjust=.1))
  
  
  dataFrame_barplot <- dataFrame_barplot+labs(fill=axis_y)+xlab(axis_x)
  
  #P <- ggplot(F_melt, aes(x=variable, y= value, fill = status)) +geom_bar(stat = "identity")
  return (dataFrame_barplot)
}
