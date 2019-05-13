#File that contains all functions 

#-----------Plot for Qvalue and Pvalue-------------
# Plot P_vale for every type of bread
plot_Pvalue <- function(stats_df){
  pp_Pval <- ggplot(stats_df, aes(x = Pvalue))+
    geom_histogram(bins = 100,fill = "steelblue") +facet_wrap(~ lable,nrow = 2)+
    ggtitle("Pvalue Stats for 8 different beads")+theme_classic()+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 12,face="bold",angle = 90),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size=18,face="bold"),
          legend.position = "top",plot.margin=unit(c(1,1,1,1),"cm"))
  return (pp_Pval)
}

# Plot Q_vale for every type of bread
plot_Qvalue <- function(stats_df){
  pp_Qval <- ggplot(stats_df, aes(x = Qvalue))+
    geom_histogram(bins = 100,fill = "steelblue") +facet_wrap(~ lable,nrow = 2)+
    ggtitle("Qvalue Stats for 8 different beads")+theme_classic()+
    theme(text = element_text(size = 18,face="bold"),
          axis.text.x = element_text(size = 12,face="bold",angle = 90),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size=18,face="bold"),
          legend.position = "top",plot.margin=unit(c(1,1,1,1),"cm"))
  return (pp_Qval)
}

#------------Plot FPKM distribution: Get necessay data for plotting -------------------
#subfuncion for getting count with certain threshold for each gene
#count the FPKM with different FPKM threshold
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

#subfuncion for getting count with certain threshold for each gene for a subset of type
#count the FPKM with different FPKM threshold
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

# Main function
#get the merged FPKM table for plotting the stacked Bar plot
get_FPKMcount_table <- function(sample_info,startCol,endCol,levelVec,RNA_annotation){
  FPKM_0 <- geneCount_FPKM(sample_info,startCol,endCol,0)
  FPKM_0_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_0)))%>% rownames_to_column(.,"RNA_type")
  write.csv(FPKM_0_selected,"FPKM_0_GeneCount.csv",col.names = FALSE)
  FPKM_0_selected$Condition <- "FPKM_0"
  
  FPKM_5 <- geneCount_FPKM(sample_info,startCol,endCol,5)
  FPKM_5_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_5)))%>% rownames_to_column(.,"RNA_type")
  write.csv(FPKM_5_selected,"FPKM_5_GeneCount.csv")
  FPKM_5_selected$Condition <- "FPKM_5"
  
  FPKM_all<- rbind(FPKM_0_selected,FPKM_5_selected)
  
  FPKM_10 <- geneCount_FPKM(sample_info,startCol,endCol,10)
  FPKM_10_selected <- as.data.frame(t(geneCount_FPKM_select(FPKM_10)))%>% rownames_to_column(.,"RNA_type")
  write.csv(FPKM_10_selected,"FPKM_10_GeneCount.csv",col.names = FALSE)
  FPKM_10_selected$Condition <- "FPKM_10"
  
  FPKM_all <- rbind(FPKM_all,FPKM_10_selected) 
  
  FPKM_all_melt <- melt(FPKM_all,id.vars = c("RNA_type","Condition"))
  FPKM_all_melt$Condition <- factor(FPKM_all_melt$Condition,levels = c("FPKM_0","FPKM_5","FPKM_10"))
  FPKM_all_melt$RNA_type <- factor(FPKM_all_melt$RNA_type,levels = FPKM_0_selected$RNA_type)
  FPKM_all_melt$variable <- factor(FPKM_all_melt$variable,levels =rev(levelVec))
  
  return (FPKM_all_melt)
}

#-------------------DEseq---------------------
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
  countMatrix <- sapply(df,as.numeric) %>% as.matrix(.)
  condition <- factor(c(rep("Surface",ss), rep("Total", ts)),levels = c("Total","Surface"))
  col_data <- data.frame(row.names = colnames(countMatrix),condition)
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,colData = col_data,design=~condition)
  
  #perform DEseq without pre-filtering
  dds <- DESeq(dds,betaPrior=FALSE)
  res_dds <- results(dds) %>% as.data.frame(.) 
  res_dds$Gene_name <- rownames(df)
  #%>% rownames_to_column(.,"Gene_name")
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

#?write_csv
get_DEgenes_info <- function(DEres_df,log2FoldChaneg_thresh,Pval_thresh,anno_info,fileName){
  cand_names<- subset(DEres_df,(log2FoldChange>log2FoldChaneg_thresh)&
                                  (padj<Pval_thresh)) %>% as.data.frame(.)
  gene_ID <- cand_names$Gene_name %>% as.character(.)
  
  res_df <- merge (cand_names,anno_info, by.x = "Gene_name", by.y = "ENSEMBL")
  #res_df_final <- merge(res_df,Biotype_df,by.x  = "Gene_name",by.y = "gene_id")
  if (nrow(res_df) != 0){
    write.table(res_df,fileName,row.names = FALSE)
  }
  return (res_df)
}