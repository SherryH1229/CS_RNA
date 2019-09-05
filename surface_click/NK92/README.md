Untitled
================
Xuerui Huang
7/16/2019

laod package
============

``` r
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(grid)
library(RColorBrewer)
library(gage)
library(gageData)
library(clusterProfiler)
library(org.Hs.eg.db)

# libaray for plot PCA
require("ggrepel")

# libaray for plot heatmap
library(ComplexHeatmap)
library(circlize)

source("~/dataOS/CS_RNA/Annotation.R")
source("~/dataOS/CS_RNA/Functions.R")
```

Basic info
==========

FeatureCount is default to be non-strandSpecific, HTseq is default to be strandSpecific. Use featureCount with Strand Specificity for future analysis

Load data

``` r
# Starnd specific featureCount res
STAR_SS_counts <- read.csv("~/dataOS/CS_RNA/surface_click/NK92/Data_pool/datapool/countsMatrix_SS.txt",
                            sep = "\t",skip = 1) %>% .[,c(1,7:ncol(.))] %>%
                            set_colnames(.,gsub(".bam","",colnames(.))) %>%
                            set_colnames(.,gsub("\\.","_",colnames(.)))
```

Analysis
========

PCA
---

Format dataframe

``` r
# Format 
RNAclick_count_df <- column_to_rownames(STAR_SS_counts,"Geneid") %>% .[rowSums(.)> 0,]
RNAclick_count_df$NK92_Tot_Click_Ctrl_3 <- NULL
```

plot PCA by using plotPCA() function

``` r
# setup DEseq dataframe with condition

# with negative control
# dds_tot <- DESeqDataSetFromMatrix(countData = as.matrix(RNAclick_count_df),
#                                   colData = data.frame(row.names =colnames(as.matrix(RNAclick_count_df)),            
#                                                        condition=c(rep("Surface",5),rep("Neg_ctrl",5),rep("Total",3))),
#                                   design=~condition)

# setup dataframe with only surface samples and total samples
ST_df <- format_samples(RNAclick_count_df,c("NK92_Sur_Click_1","NK92_Sur_Click_2","NK92_Sur_Click_3",
                                            "NK92_Sur_Click_4","NK92_Sur_Click_5",
                                            "NK92_Tot_Click_1","NK92_Tot_Click_3",
                                            "NK92_Tot_Click_4"))
# no negative controls, output to plot samples only
dds_tot <- DESeqDataSetFromMatrix(countData = as.matrix(ST_df),
                                  colData = data.frame(row.names =colnames(as.matrix(ST_df)),            
                                                       condition=c(rep("Surface",5),rep("Total",3))),
                                  design=~condition)
# Plot PCA
z <- plotPCA(rlog(dds_tot))+geom_label_repel(aes(label = name),size = 2.5)+
  theme_classic2()+
  theme(text = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size = 12,face="bold"),
        plot.title = element_text(size=12,face="bold"))+
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  ggtitle("PCA of all the Surface_click Samples")

z
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
ggsave("./plots/PCA_NK92_samplesOnly.png",z,height= 4, width = 5)
```

DEseq
-----

All Surface VS All Total

``` r
# filter out genes with 0 signals across all sampples 
RNAclick_count_df <- column_to_rownames(STAR_SS_counts,"Geneid") %>% .[rowSums(.)> 0,]

# remove thee negative control for total sample 
RNAclick_count_df$NK92_Tot_Click_Ctrl_3 <- NULL


# setup dataframe with only surface samples and total samples
ST_df <- format_samples(RNAclick_count_df,c("NK92_Sur_Click_1","NK92_Sur_Click_2","NK92_Sur_Click_3",
                                            "NK92_Sur_Click_4","NK92_Sur_Click_5",
                                            "NK92_Tot_Click_1","NK92_Tot_Click_3",
                                            "NK92_Tot_Click_4"))
# perform DEseq by using function perform_DEseq
ST_DEseq <- perform_DEseq(ST_df,5,3)
```

    ## converting counts to integer mode

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
# get the information of differentially expressed gene by applying the threshold of log2FoldChange larger than 2 
# and qvalue smaller than 0.01. the output will be in dataframe and csv format
# two different threshold 
res.ST.2_0.01 <- get_DEgenes_info(ST_DEseq,2,0.01,anno_info_Hs,"ST_2_0.01.csv")
res.ST.1_0.05 <- get_DEgenes_info(ST_DEseq,1,0.05,anno_info_Hs,"ST_1_0.05.csv")
```

Visualization of results
------------------------

plot by using complexHeatmap package

``` r
# remove ribosomal RNA related genes
res.ST_select <- res.ST.2_0.01[- grep("riboso", res.ST.2_0.01$GENENAME),]

# reformat
heatmap.df <- ST_df[rownames(ST_df) %in% res.ST_select$Gene_name, ]
temp<- as.data.frame(res.ST_select$Gene_name,res.ST_select$SYMBOL) %>% rownames_to_column(.,"SYMBOL") %>% set_colnames(c("SYMBOL","ENSEMBL")) %>% na.omit(.)
heatmap.df <- merge(heatmap.df,temp,by.x ="row.names",by.y = "ENSEMBL") %>% column_to_rownames(.,"SYMBOL")
heatmap.df$Row.names <- NULL 
heatmap.df <- t(apply(heatmap.df , 1, function(x) x/sum(x))) # normalized heatmap by row

# setup heatmap colors
col_fun = colorRamp2(c(1, 0), c( "midnightblue","white"))
col_fun(seq(-3, 3))
```

    ## [1] "#FFFFFFFF" "#FFFFFFFF" "#FFFFFFFF" "#FFFFFFFF" "#191970FF" "#191970FF"
    ## [7] "#191970FF"

``` r
ha = columnAnnotation(foo = anno_text(c(rep("Sample",4),rep("Total",3)), location = 0.5, just = "center",
    gp = gpar(fill = rep(2:4, each = 4), col = "white", border = "black"),
    width = max_text_width(month.name)*1))

fa = rowAnnotation(foo = anno_mark(at = c(12,26,32,60), labels = rownames(heatmap.df)[c(12,26,32,60)]))

pp <- Heatmap(heatmap.df, name = "Scaled Counts", #row_order = rownames(heatmap.df),
              column_order = colnames(heatmap.df),
        column_names_side = "bottom",#column_names_rot = 45,
        show_row_names = FALSE,
        column_split = c(rep("Surface RNA",5),
                         rep("Total RNA",3)),column_gap = unit(3, "mm"),
        right_annotation = fa,
        col = col_fun)
pp
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
png("./plots/STHeatmap_2_0.01_select.png",width = 500, height = 500,
    units = "px", pointsize = 12, bg = "white")
draw(pp)
dev.off()
```

    ## png 
    ##   2

``` r
knitr::include_graphics("./plots/STHeatmap_2_0.01_select.png")
```

<img src="./plots/STHeatmap_2_0.01_select.png" width="500" />

``` r
# For sample pool with negative controls 
# ha = columnAnnotation(foo = anno_text(c(rep("Sample",4),rep("Neg_ctrl",4),rep("Total",4)), location = 0.5, just = "center",
#     gp = gpar(fill = rep(2:4, each = 4), col = "white", border = "black"),
#     width = max_text_width(month.name)*1))
# fa = rowAnnotation(foo = anno_mark(at = c(12,26,32,60), labels = rownames(heatmap.df)[c(12,26,32,60)]))
# 
# pp <- Heatmap(heatmap.df, name = "Scaled Counts", #row_order = rownames(heatmap.df),
#               column_order = colnames(heatmap.df),
#         column_names_side = "bottom",#column_names_rot = 45,
#         show_row_names = FALSE,
#         column_split = c(rep("Sample",5),rep("Neg_control",5),
#                          rep("Total",3)),column_gap = unit(3, "mm"),
#         right_annotation = fa,
#         col = col_fun)
```

GO and KEGG Enrichment

``` r
#GO on geneset of 2 and 0.01
res.ST_select <- res.ST.2_0.01[- grep("riboso", res.ST.2_0.01$GENENAME),]
ego <- enrichGO(gene = res.ST_select$Gene_name, OrgDb='org.Hs.eg.db', keyType ="ENSEMBL",ont="all",
                pvalueCutoff = 0.05, 
                pAdjustMethod = "fdr",qvalueCutoff = 0.05)
#ego@result
ggsave("./plots/GO_ST_2_0.01_selected.png",dotplot(ego, split="ONTOLOGY") + 
         facet_grid(ONTOLOGY~., scale="free")+
         theme(text = element_text(size=12),
               plot.title = element_text(size=15,face="bold"))+
         ggtitle("GO Enrichment"),height= 6, width = 10)
```

    ## wrong orderBy parameter; set to default `orderBy = "x"`

``` r
knitr::include_graphics("./plots/GO_ST_2_0.01_selected.png")
```

<img src="./plots/GO_ST_2_0.01_selected.png" width="3000" />

``` r
res.ST_select <- res.ST.1_0.05[- grep("riboso", res.ST.1_0.05$GENENAME),]
#Kegg on gene set of 1 and 0.05
gene.uniport <- bitr(res.ST_select$Gene_name, fromType = "ENSEMBL",
        toType = c("UNIPROT"),
        OrgDb = org.Hs.eg.db)
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## Warning in bitr(res.ST_select$Gene_name, fromType = "ENSEMBL", toType =
    ## c("UNIPROT"), : 11.74% of input gene IDs are fail to map...

``` r
eke <- enrichKEGG(gene = gene.uniport$UNIPROT, organism = "hsa", keyType ="uniprot", pvalueCutoff = 0.01,qvalueCutoff = 0.05)
ggsave("./plots/KEGG_ST_2_0.01_selected.png",dotplot(eke)+
         theme(text = element_text(size=12),
               plot.title = element_text(size=15,face="bold"))+
         ggtitle("KEGG Enrichment"),height= 4, width = 8) 
```

    ## wrong orderBy parameter; set to default `orderBy = "x"`

``` r
knitr::include_graphics("./plots/KEGG_ST_2_0.01_selected.png")
```

<img src="./plots/KEGG_ST_2_0.01_selected.png" width="2400" />
