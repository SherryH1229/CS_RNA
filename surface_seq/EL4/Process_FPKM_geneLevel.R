library(magrittr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(reshape2)
#--------------------process original FPKM file----------------------
#setwd("~/dataOS/CS_RNA")
sample_info <- read.csv("/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt",sep="\t")
sample_info$EL4.NP.ER.Son.1 <- sample_info$EL4.NP.ER.Son.1+sample_info$EL4.NP.ER.Son.2
sample_info$EL4.NP.ER.Son.2 <- NULL
colnames(sample_info)[6]<-"EL4.NP.ER.Son"


#FPKM_table_selected<-sample_info %>% .[!duplicated(.$tracking_id),]
#rownames(FPKM_table_selected)<-FPKM_table_selected$tracking_id
#FPKM_table_selected$tracking_id<-NULL
#FPKM_table_selected<- FPKM_table_selected[,c(5:12)] 

#FPKM_table_selected$EL4.NP.ER.Son.1<- FPKM_table_selected$EL4.NP.ER.Son.2+FPKM_table_selected$EL4.NP.ER.Son.1
#FPKM_table_selected$EL4.NP.ER.Son.2 <- NULL
#colnames(FPKM_table_selected)[1]<-"EL4.NP.ER.Son"
#-------------------Plot heatmap----------------------------
FPKM_matrix<-as.matrix(FPKM_table_selected)
head(FPKM_matrix)
#png("FPKM_heatmap_ploy.png")
#pdf("FPKM_heatmap_ploy.pdf",width=10,height = 10)
#png(filename="FPKM_heatmap_plot.png", width=900, height = 900,bg="white")
#par(mar=c(5,6,4,1)+.1)
#par(mar=c(8, 4, 2, 2) + 0.1)
png(filename="FPKM_heatmap_plot.png")
op <- par(mar = c(10,4,4,2) + 0.1)
#par(mar=c(5.1,4.1,4.1,2.1))
heatmap(FPKM_matrix)
dev.off()

#-------------------FPKM_count-----------------------------
RNA_annotation <- as.character(unique(sample_info$Biotype))
RNA_mainType <- grep("RNA",RNA_annotation,value=TRUE)
RNA_mainType
setwd("~/dataOS/CS_RNA/surface_seq/barplot")


FPKM_0 <- geneCount_FPKM(sample_info,6,12,0)
FPKM_0_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_0)))%>% rownames_to_column(.,"RNA_type")
write.csv(FPKM_0_selected,"FPKM_0_GeneCount.csv")
FPKM_0_selected$Condition <- "FPKM_0"


FPKM_5 <- geneCount_FPKM(sample_info,6,12,5)
FPKM_5_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_5)))%>% rownames_to_column(.,"RNA_type")
write.csv(FPKM_5_selected,"FPKM_5_GeneCount.csv")
FPKM_5_selected$Condition <- "FPKM_5"

FPKM_all<- rbind(FPKM_0_selected,FPKM_5_selected)

FPKM_10 <- geneCount_FPKM(sample_info,6,12,10)
FPKM_10_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_10)))%>% rownames_to_column(.,"RNA_type")
write.csv(FPKM_10_selected,"FPKM_10_GeneCount.csv")
FPKM_10_selected$Condition <- "FPKM_10"

FPKM_all <- rbind(FPKM_all,FPKM_10_selected) 

FPKM_all_melt <- melt(FPKM_all,id.vars = c("RNA_type","Condition"))
FPKM_all_melt$Condition <- factor(FPKM_all_melt$Condition,levels = c("FPKM_0","FPKM_5","FPKM_10"))
FPKM_all_melt$RNA_type <- factor(FPKM_all_melt$RNA_type,levels = FPKM_0_selected$RNA_type)
FPKM_all_melt$variable <- factor(FPKM_all_melt$variable,levels = c("EL4.NP.ER.Son","EL4.NP.SAL.Son","EL4.NP.SAL.Ext","EL4.IC.SAL.1","EL4.IC.SAL.2","EL4.Tot","EL4.Tot.ER"))

P <- ggplot(FPKM_all_melt, aes(x=variable, y= value, fill = RNA_type)) +
  geom_bar(stat = "identity") + coord_flip()
P <- P + facet_grid(Condition ~.) +theme(axis.text.x = element_text(angle=-40, hjust=.1)) +
  scale_fill_discrete(name = "Assigned_status") + ylab("Count")+ xlab("Sample_ID")+
  scale_fill_brewer(palette="Set3")+ggtitle("FeatureCount_Assign_Distribution")
P
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
