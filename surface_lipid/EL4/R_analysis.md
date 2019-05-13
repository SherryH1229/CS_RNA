EL4\_lipid\_pullDown
================
Xuerui Huang
5/5/2019

-   [Load requried package](#load-requried-package)
-   [Load and Process Data](#load-and-process-data)
-   [DEseq](#deseq)
    -   [Overall filtering](#overall-filtering)
    -   [Independent Filtering](#independent-filtering)

Load requried package
=====================

``` r
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(grid)

source("~/dataOS/CS_RNA/Annotation.R")
source("~/dataOS/CS_RNA/Functions.R")
```

Load and Process Data
=====================

Load count data

``` r
#batch 1
EL4_count_1 <- read.csv("/dataOS/frankyan/OTHERS/lipid_RNA/Results/pipeOutput/mm_lipidRNA-L1_2/countsTable/mm_lipidRNA-L1_2.EX.counts.table",sep = "\t") %>% .[,c(1,7:ncol(.))] %>% column_to_rownames(.,"Geneid")

#batch 2
EL4_count_2 <- read.csv("/dataOS/frankyan/OTHERS/lipid_RNA/Results/pipeOutput/mm_lipidRNA-L2/countsTable/mm_lipidRNA-L2.EX.counts.table",sep = "\t")  %>% .[,c(1,7:ncol(.))] %>% column_to_rownames(.,"Geneid") %>% .[,-c(7,8)]

# Reform data
data_list <- c()
for (i in seq(1, ncol(EL4_count_1), by=2)){
  data_list <- c(data_list,EL4_count_1[,c(i,i+1)],EL4_count_2[,c(i,i+1)])
}
count_df <- as.data.frame(data_list) %>% set_rownames(.,gsub("\\..*","", rownames(EL4_count_1)))
```

Process count data: Stats

``` r
#Count occurance of 0
zero_count <- rowSums(count_df == 0)
perc_count <- table(zero_count) %>% as.data.frame(.)
perc_count$Percentage <- round((perc_count$Freq/sum(perc_count$Freq))*100 , 1) %>% paste0(.,"%")

# plot histogram of zero count
pp_bar <- ggplot(perc_count,aes(x = zero_count,y = Freq))+
  geom_bar(stat="identity", fill="steelblue")+theme_classic()+
  geom_text(aes(label=Percentage),hjust = -0.05, color="black", size=2.5,angle=90)+
  ggtitle("Histogram of Zero Count")+
  scale_y_continuous(expand = expand_scale(mult = c(0, .13)))+
  #xlab("Occur_Freq_thresh")+ggtitle("Count of Genes Across Different Num of Samples: Batch1")+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 12,face="bold",angle=90),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size=18,face="bold"),
        legend.position = "top",plot.margin=unit(c(1,1,1,1),"cm"))
pp_bar
```

![](R_analysis_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
ggsave("Hist_0_occur.png",pp_bar,height = 6 , width = 10)
```

DEseq
=====

Perform DEseq and generate data required for analysis

Overall filtering
-----------------

``` r
#filter 
fil_count_df <- count_df[rowSums(count_df == 0)<=34,]
temp_df <- fil_count_df[,c(1:4,5:8)]
temp_DEres_df <- perform_DEseq(temp_df,4,4)

OutputVec_ENSEMBL <- c()
OutputVec_SYMBOL <- c()
File_Name_vec <- c()

remove(stats_df)
#OutputVec_value <- c()
for (i in seq(5, ncol(fil_count_df), by=4)){
  fileName <- gsub("\\..*","",colnames(fil_count_df)[i])
  
  temp_df <- fil_count_df[,c(1:4,i:(i+3))]
  temp_DEres_df <- perform_DEseq(temp_df,4,4)
  temp_stats <- cbind(temp_DEres_df$pvalue,temp_DEres_df$padj) %>% as.data.frame(.)
  temp_stats$Tag <- fileName
  
  if (!exists("stats_df")){
    stats_df <- temp_stats
  }
  else{
    stats_df <- rbind(stats_df,temp_stats)
  }
  #break
  outputFileName <-  paste0(fileName,"_cands.csv")
  File_Name_vec <- c(File_Name_vec,fileName)
  
  temp_cand <- get_DEgenes_info(temp_DEres_df,2,0.5,anno_info_Mm,outputFileName)
  OutputVec_ENSEMBL <- c(OutputVec_ENSEMBL, temp_cand$Gene_name %>% as.data.frame(.))
  OutputVec_SYMBOL <- c(OutputVec_SYMBOL,temp_cand$SYMBOL %>% as.data.frame(.))
}
#File_title <- gsub("\\..*","",colnames(EL4_count))[seq(3,ncol(EL4_count), by=2)]
#File_title

#Format output res into dataframe
# Op_ENSEMBL <- sapply(OutputVec_ENSEMBL, '[', seq(max(sapply(OutputVec_ENSEMBL, length)))) %>% 
#   as.data.frame(.) %>% set_colnames(.,File_title)
# Op_SYMBOL <- sapply(OutputVec_SYMBOL, '[', seq(max(sapply(OutputVec_SYMBOL, length)))) %>% 
#   as.data.frame(.) %>% set_colnames(.,File_title)
```

Plot

``` r
fil_stats_df <- na.omit(stats_df) %>% set_colnames(.,c("Pvalue","Qvalue","lable")) %>% as.data.frame(.)

pp_P <- plot_Pvalue(fil_stats_df)
pp_P
```

![](R_analysis_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
ggsave("Hist_Pvalue.png",pp_P,height = 6 , width = 10)

pp_Q <- plot_Qvalue(fil_stats_df)
pp_Q
```

![](R_analysis_files/figure-markdown_github/unnamed-chunk-5-2.png)

``` r
ggsave("Hist_Qvalue.png",pp_Q,height = 6 , width = 10)
```

Independent Filtering
---------------------

``` r
#filter 
OutputVec_ENSEMBL <- c()
OutputVec_SYMBOL <- c()
File_Name_vec <- c()
remove(stats_df)

for (i in seq(5, ncol(fil_count_df), by=4)){
  fileName <- gsub("\\..*","",colnames(fil_count_df)[i])
  
  temp_df <- fil_count_df[,c(1:4,i:(i+3))]
  temp_df <- temp_df[rowSums(temp_df == 0)<=6,]
  temp_DEres_df <- perform_DEseq(temp_df,4,4)
  temp_stats <- cbind(temp_DEres_df$pvalue,temp_DEres_df$padj) %>% as.data.frame(.)
  temp_stats$Tag <- fileName
  
  if (!exists("stats_df")){
    stats_df <- temp_stats
  }
  else{
    stats_df <- rbind(stats_df,temp_stats)
  }

  outputFileName <-  paste0(fileName,"_cands.csv")
  File_Name_vec <- c(File_Name_vec,fileName)
  
  temp_cand <- get_DEgenes_info(temp_DEres_df,2,0.5,anno_info_Mm,outputFileName)
  OutputVec_ENSEMBL <- c(OutputVec_ENSEMBL, temp_cand$Gene_name %>% as.data.frame(.))
  OutputVec_SYMBOL <- c(OutputVec_SYMBOL,temp_cand$SYMBOL %>% as.data.frame(.))
}
```

Plot

``` r
fil_stats_df <- na.omit(stats_df) %>% set_colnames(.,c("Pvalue","Qvalue","lable")) %>% as.data.frame(.)

pp_P <- plot_Pvalue(fil_stats_df)
pp_P
```

![](R_analysis_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
ggsave("Hist_Pvalue_ind.png",pp_P,height = 6 , width = 10)

pp_Q <- plot_Qvalue(fil_stats_df)
pp_Q
```

![](R_analysis_files/figure-markdown_github/unnamed-chunk-7-2.png)

``` r
ggsave("Hist_Qvalue_ind.png",pp_Q,height = 6 , width = 10)
```
