library(magrittr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(reshape2)
library(RColorBrewer)
#--------------------process original FPKM file Surface Click----------------------
#setwd("~/dataOS/CS_RNA")
sample_info_1 <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_click/NK92/20190408/20190408_FPKMs.txt",sep = "\t")
sample_info_2 <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_click/NK92/20190409/201904010_FPKMs.txt",sep="\t")

temp_sample_info <- data.frame(sample_info_2$NK92.Sur.Click.2_FPKM)
colnames(temp_sample_info) <- c("NK92.Sur.Click.2_FPKM")
sample_info <- cbind(sample_info_1,temp_sample_info)
#colnames(sample_info)[8:13]<- gsub("(\\w+)\\.(\\w+)\\.(\\w+)\\.([0-9]+)_FPKM","\\1_\\2_\\3_\\4",colnames(sample_info)[8:13])
#colnames(sample_info)[9:11]<- c("NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2")
colnames(sample_info)[8:13]<- c("NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2",
                                "NK92_Tot_Click_1","NK92_Tot_Click_3","NK92_Sur_Click_2")

#sub(".*\\s([0-9]+)\\snomination.*$", "\\1", awards)
RNA_annotation <- as.character(unique(sample_info$Biotype))
RNA_annotation
RNA_mainType <- grep("RNA",RNA_annotation,value=TRUE)
RNA_mainType

source("/dataOS/sherry/CS_RNA/Functions.R")
test <- get_FPKMcount_table(8,13,colnames(sample_info)[8:13],RNA_annotation)

FPKM_0 <- geneCount_FPKM(sample_info,8,13,0)
FPKM_0_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_0)))%>% rownames_to_column(.,"RNA_type")
#write.csv(FPKM_0_selected,"FPKM_0_GeneCount.csv")
FPKM_0_selected$Condition <- "FPKM_0"

FPKM_5 <- geneCount_FPKM(sample_info,8,13,5)
FPKM_5_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_5)))%>% rownames_to_column(.,"RNA_type")
#write.csv(FPKM_5_selected,"FPKM_5_GeneCount.csv")
FPKM_5_selected$Condition <- "FPKM_5"

FPKM_all<- rbind(FPKM_0_selected,FPKM_5_selected)

FPKM_10 <- geneCount_FPKM(sample_info,8,13,10)
FPKM_10_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_10)))%>% rownames_to_column(.,"RNA_type")
#write.csv(FPKM_10_selected,"FPKM_10_GeneCount.csv")
FPKM_10_selected$Condition <- "FPKM_10"

FPKM_all <- rbind(FPKM_all,FPKM_10_selected) 

FPKM_all_melt <- melt(FPKM_all,id.vars = c("RNA_type","Condition"))
FPKM_all_melt$Condition <- factor(FPKM_all_melt$Condition,levels = c("FPKM_0","FPKM_5","FPKM_10"))
FPKM_all_melt$RNA_type <- factor(FPKM_all_melt$RNA_type,levels = FPKM_0_selected$RNA_type)
#FPKM_all_melt$variable <- gsub("_FPKM","",FPKM_all_melt$variable)
FPKM_all_melt$variable <- factor(FPKM_all_melt$variable,levels =
                                   rev(c("NK92_Sur_Click_1","NK92_Sur_Click_1_1","NK92_Sur_Click_1_2",
                                         "NK92_Sur_Click_2","NK92_Tot_Click_1","NK92_Tot_Click_3")))
                                   # rev(c("K562_Seq_Sur_1","K562_Seq_Sur_2","K562_Seq_Sur_3","K562_Seq_Sur_4","K562_Seq_Sur_5",
                                   #       "K562_Seq_Tot_1","K562_Seq_Tot_2")))



P <- ggplot(FPKM_all_melt, aes(x=variable, y= value, fill = RNA_type)) +
  geom_bar(stat = "identity") + coord_flip()

colourCount = length(unique(FPKM_all_melt$RNA_type))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
P <- P + facet_grid(Condition ~.) +theme() +
  scale_fill_discrete(name = "Assigned_status") + ylab("Count")+ xlab("Sample_ID")+
  #scale_fill_brewer(palette="Set3")+
  ggtitle("FeatureCount_Assign_Distribution")+
  theme(text = element_text(size = 18,face="bold"),
        #axis.text.x = element_text(angle=-40, hjust=.1),
        #axis.text.y = element_text(angle=-40, hjust=.1),
        axis.text.x = element_text(size = 18,face="bold",angle=-40, hjust=.1),
        axis.text.y = element_text(size = 18,face="bold"),
        plot.title = element_text(size=15,face="bold"))+
  scale_fill_manual(values = getPalette(colourCount))
P

#setwd ("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_seq/K562/")
setwd ("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_click/NK92/")
ggsave("FPKM_NK92_190421.png",P,height= 6 , width = 10)






#test <- data.frame(sample_info_2$NK92.Sur.Click.2_FPKM,sample_info_2$NK92.Tot.Click.Ctrl.3_FPKM)
#"/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt"
#sample_info <- read.csv("/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt",sep="\t")

# sample_info_1 <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_seq/K562/20190412/20190412_FPKMs.txt",sep="\t")
# sample_info_2 <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/sherry/CS_RNA/surface_seq/K562/20190421/20190421_FPKMs.txt",sep  = "\t")
# sample_info_1$X <- NULL
# sample_info_1$K562_TotalRNA_FPKM <- NULL
# sample_info <- cbind(sample_info_1,sample_info_2$K562_TotalRNA_FPKM) %>% .[!duplicated(.$gene_id),]
# colnames(sample_info)[8:14] <- c("K562_Seq_Sur_2","K562_Seq_Sur_4","K562_Seq_Sur_1","K562_Seq_Sur_3","K562_Seq_Sur_5",
#                                  "K562_Seq_Tot_2","K562_Seq_Tot_1")



# FPKM_table_selected<-sample_info %>% .[!duplicated(.$gene_id),]
# rownames(FPKM_table_selected)<-FPKM_table_selected$gene_id
# FPKM_table_selected$gene_id<-NULL
#FPKM_table_selected$X <- NULL
#FPKM_table_selected$K562_TotalRNA_FPKM <- NULL
#FPKM_table_selected<- FPKM_table_selected[,c(1:3,14:18)] 
#test <- sample_info[,c(14:18)]
#-------------------FPKM_bar_plot----------

#setwd("~/dataOS/CS_RNA/surface_click/plots/barplot")
FPKM_all_melt$variable <- factor(FPKM_all_melt$variable,levels =
                                   rev(c("K562_Seq_Sur_1","K562_Seq_Sur_2","K562_Seq_Sur_3","K562_Seq_Sur_4","K562_Seq_Sur_5",
                                         "K562_Seq_Tot_1","K562_Seq_Tot_2")))


#------------------------------function-------------------------------
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
  protein_coding <- FPKM_dataframe[,"protein_coding",drop=FALSE]
  FPKM_select <-cbind(FPKM_select,protein_coding)
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
 