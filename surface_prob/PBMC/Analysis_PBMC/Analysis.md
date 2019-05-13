Analysis
================
Xuerui Huang
5/9/2019

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
library(RColorBrewer)
source("~/dataOS/CS_RNA/Annotation.R")
source("~/dataOS/CS_RNA/Functions.R")
```

Load data and format
====================

``` r
# laod FPKM data
PBMC_FPKM <- read.csv("~/dataOS/CS_RNA/surface_lipid/PBMC/20190502_FPKMs.txt",sep = "\t")
colnames(PBMC_FPKM)[8:11] <- gsub("*_FPKM","",colnames(PBMC_FPKM)[8:11])
#load gene Count data
PBMC_count <- read.csv("~/dataOS/CS_RNA/surface_lipid/PBMC/htseqCountTotal.txt",sep = "\t")
```

``` r
PBMC_FPKM_count <- get_FPKMcount_table(PBMC_FPKM,8,11,colnames(PBMC_FPKM)[8:11])
```

    ## Warning in write.csv(FPKM_0_selected, "FPKM_0_GeneCount.csv", col.names =
    ## FALSE): attempt to set 'col.names' ignored

    ## Warning in write.csv(FPKM_10_selected, "FPKM_10_GeneCount.csv", col.names =
    ## FALSE): attempt to set 'col.names' ignored

``` r
colourCount = length(unique(PBMC_FPKM_count$RNA_type))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

P <- ggplot(PBMC_FPKM_count, aes(x=variable, y= value, fill = RNA_type)) +
  geom_bar(stat = "identity") + coord_flip()+facet_grid(Condition ~.) +theme() +
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
```

    ## Scale for 'fill' is already present. Adding another scale for 'fill',
    ## which will replace the existing scale.

``` r
P
```

![](Analysis_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
ggsave("FPKM_PBMC.png",P,height= 6 , width = 10)
```
