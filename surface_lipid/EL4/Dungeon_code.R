## Batch 1
```{r, results='hide', message=FALSE, warning=FALSE}
#format table
EL4_count <- EL4_exon_count[,c(1,7:ncol(EL4_exon_count))]
EL4_count <- column_to_rownames(EL4_count,"Geneid")

#Output candidate genes (in files) for each set of comparison
#Output a merged result for ENSELBl and gene_SYMBL
#Used pvalue instead of padj
OutputVec_ENSEMBL <- c()
OutputVec_SYMBOL <- c()
for (i in seq(3, ncol(EL4_count), by=2)){
  temp_df <- EL4_count[c(i,i+1,1,2)]
  temp_DEres_df <- perform_DEseq(temp_df,2,2)
  outputFileName <- gsub("_.*","",colnames(EL4_count)[i]) %>% paste0(.,"_cands.csv")
  temp_cand <- get_DEgenes_info(temp_DEres_df,2,0.5,anno_info_Mm,outputFileName)
  OutputVec_ENSEMBL <- c(OutputVec_ENSEMBL, temp_cand$Gene_name %>% as.data.frame(.))
  OutputVec_SYMBOL <- c(OutputVec_SYMBOL,temp_cand$SYMBOL %>% as.data.frame(.))
}

#get output colnames
File_title <- gsub("_.*","",colnames(EL4_count))[seq(3,ncol(EL4_count), by=2)]

#Format output res into dataframe
Op_ENSEMBL <- sapply(OutputVec_ENSEMBL, '[', seq(max(sapply(OutputVec_ENSEMBL, length)))) %>% 
  as.data.frame(.) %>% set_colnames(.,File_title)
Op_SYMBOL <- sapply(OutputVec_SYMBOL, '[', seq(max(sapply(OutputVec_SYMBOL, length)))) %>% 
  as.data.frame(.) %>% set_colnames(.,File_title)
```
## Batch 2
```{r, results='hide', message=FALSE, warning=FALSE}
#format table
EL4_count_2 <- EL4_exon_count_2[,c(1,7:ncol(EL4_exon_count_2))] %>% column_to_rownames(.,"Geneid")

#Output candidate genes (in files) for each set of comparison
#Output a merged result for ENSELBl and gene_SYMBL
#Used pvalue instead of padj
OutputVec_ENSEMBL_2 <- c()
OutputVec_SYMBOL_2 <- c()
seq_list <- seq(3, ncol(EL4_count_2), by=2) %>% .[-(3)]

for (i in seq_list){
  print (i)
  temp_df <- EL4_count_2[c(i,i+1,1,2)]
  temp_DEres_df <- perform_DEseq(temp_df,2,2)
  outputFileName <- gsub("_.*","",colnames(EL4_count_2)[i]) %>% paste0(.,"_cands.csv")
  print (outputFileName)
  temp_cand <- get_DEgenes_info(temp_DEres_df,2,0.5,anno_info_Mm,outputFileName)
  OutputVec_ENSEMBL_2 <- c(OutputVec_ENSEMBL_2, temp_cand$Gene_name %>% as.data.frame(.))
  OutputVec_SYMBOL_2 <- c(OutputVec_SYMBOL_2,temp_cand$SYMBOL %>% as.data.frame(.))
}


#get output colnames
File_title_2 <- gsub("_.*","",colnames(EL4_count_2))[seq_list]

#Format output res into dataframe
Op_ENSEMBL_2 <- sapply(OutputVec_ENSEMBL_2, '[', seq(max(sapply(OutputVec_ENSEMBL_2, length)))) %>% 
  as.data.frame(.) %>% set_colnames(.,File_title_2)
Op_SYMBOL_2 <- sapply(OutputVec_SYMBOL_2, '[', seq(max(sapply(OutputVec_SYMBOL_2, length)))) %>% 
  as.data.frame(.) %>% set_colnames(.,File_title_2)
```
# Plot result

## Count
```{r}
#Get all the genes in batch 1
Uni_gene <- c()
for (i in c(1:ncol(Op_SYMBOL))){
  Uni_gene <- c(Uni_gene,Op_SYMBOL[,i] %>% .[!is.na(.)] %>% as.vector(.))
}

#Get all the genes in batch 2
# Uni_gene_2 <- c()
# for (i in c(1:ncol(Op_SYMBOL_2))){
#   Uni_gene_2 <- c(Uni_gene_2,Op_SYMBOL_2[,i] %>% .[!is.na(.)] %>% as.vector(.))
# }

#gene cands duplicated occured genes
dup_genes <- table(Uni_gene) %>% as.data.frame(.) %>% subset(.,.$Freq>4) %>% .[order(-.$Freq),]
high_occur_genes <- merge(dup_genes,anno_info_Mm,by.x = "Uni_gene",by.y = "SYMBOL")
write.table(high_occur_genes,"high_occur_genes.csv",row.names = FALSE)
dup_genes_names <- dup_genes$Uni_gene %>% as.vector(.)

# dup_genes_2 <- table(Uni_gene_2) %>% as.data.frame(.) %>% subset(.,.$Freq>4) %>% .[order(-.$Freq),]
# dup_genes_names_2 <- dup_genes_2$Uni_gene %>% as.vector(.)
#show sample list for batch 1
dup_genes_names[1:10]

#show sample list for batch 2
#dup_genes_names_2[1:10]
```



## Heatmap

### Plot heatmap batch
```{r}
#get count df for count occurace of genes in each file, 1==occur, 0==not occur
dup_genes_all <- table(Uni_gene) %>% as.data.frame(.) %>% .[order(-.$Freq),] %>% 
  .$Uni_gene %>% as.vector(.)

#dup_genes_2 <- table(Uni_gene_2) %>% as.data.frame(.) %>% subset(.,.$Freq>4) %>% .[order(-.$Freq),]
#dup_genes_names_2 <- dup_genes_2$Uni_gene %>% as.vector(.)
ct_all <- c()
for (i in c(1:length(File_title))){
  ct_occur <- c()
  for (gene in dup_genes_all){
    if (gene %in% Op_SYMBOL[,i]){
      ct_occur <- c(ct_occur,1)
    } 
    else{
      ct_occur <- c(ct_occur,0)
    }
  }
  ct_all <- c(ct_all,as.data.frame(ct_occur))
}

count_df <- cbind(dup_genes_all,as.data.frame(ct_all)) %>% as.data.frame(.) %>% set_colnames(.,c("SYMBOL",File_title))

# Run clustering
count_matrix <- as.matrix(t(count_df[,c(2:ncol(count_df))]))
rownames(count_matrix) <- File_title
dendro <- hclust(d = dist(x = count_matrix))

File_title
File_title[dendro$order]
#plot dendrogram
require("ggdendro")
count_dendro <- as.dendrogram(dendro)
dendro.plot <- ggdendrogram(data = count_dendro, rotate = TRUE)
# Preview the plot
print(dendro.plot)

#adjust plot order
count_df_melt <- melt(count_df)
count_df_melt$SYMBOL <- factor(count_df_melt$SYMBOL,levels =dup_genes_names)
count_df_melt$variable <- factor(count_df_melt$variable,levels = File_title[dendro$order])

#plot
pp <- ggplot(count_df_melt,aes(x=SYMBOL,y = variable,fill=factor(value)))+
  geom_tile()+ scale_fill_manual(values=c("lightblue","steelblue"),name = "Appear",
                                 label=c("Not Occur","Occur"))+
  theme_minimal()+ggtitle("Heatmap of Ocurrence of Genes")+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 4.5,face="bold",angle=-45,hjust=.1),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size=15,face="bold"),
        legend.position = "top")

#show plot
pp

#save plot with dandogram
png("Heatmap.png", width = 9,height = 5,units = "in",res = 1200,pointsize = 4)
grid.newpage()+
  print(pp, vp = viewport(x = 0.35, y = 0.5, width = 0.65, height = 1.0))+
  print(dendro.plot, vp =  viewport(x = 0.77, y = 0.46, width = 0.18, height = 0.72))
dev.off()
```


### Plot heatmap batch 2
```{r}
#get count df for count occurace of genes in each file, 1==occur, 0==not occur
dup_genes_all_2 <- table(Uni_gene_2) %>% as.data.frame(.) %>% .[order(-.$Freq),] %>% 
  .$Uni_gene_2 %>% as.vector(.)


ct_all <- c()
for (i in c(1:length(File_title_2))){
  ct_occur <- c()
  for (gene in dup_genes_all_2){
    if (gene %in% Op_SYMBOL_2[,i]){
      ct_occur <- c(ct_occur,1)
    } 
    else{
      ct_occur <- c(ct_occur,0)
    }
  }
  ct_all <- c(ct_all,as.data.frame(ct_occur))
}

count_df <- cbind(dup_genes_all_2,as.data.frame(ct_all)) %>% as.data.frame(.) %>%
  set_colnames(.,c("SYMBOL",File_title_2))

# Run clustering
count_matrix <- as.matrix(t(count_df[,c(2:ncol(count_df))]))
rownames(count_matrix) <- File_title_2
dendro <- hclust(d = dist(x = count_matrix))

#plot dendrogram
require("ggdendro")
count_dendro <- as.dendrogram(dendro)
dendro.plot <- ggdendrogram(data = count_dendro, rotate = TRUE)
# Preview the plot
print(dendro.plot)

#adjust plot order
count_df_melt <- melt(count_df)
count_df_melt$SYMBOL <- factor(count_df_melt$SYMBOL,levels =dup_genes_names_2)
count_df_melt$variable <- factor(count_df_melt$variable,levels = File_title_2[dendro$order])

#plot
pp <- ggplot(count_df_melt,aes(x=SYMBOL,y = variable,fill=factor(value)))+
  geom_tile()+ scale_fill_manual(values=c("lightblue","steelblue"),name = "Appear",
                                 label=c("Not Occur","Occur"))+
  theme_minimal()+ggtitle("Heatmap of Ocurrence of Genes: Batch2")+
  theme(text = element_text(size = 18,face="bold"),
        axis.text.x = element_text(size = 4.5,face="bold",angle=-45,hjust=.1),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size=15,face="bold"),
        legend.position = "top")

#show plot
pp

#save plot with dandogram
png("Heatmap_4_2.png", width = 9,height = 5,units = "in",res = 1200,pointsize = 4)
grid.newpage()+
  print(pp, vp = viewport(x = 0.35, y = 0.5, width = 0.65, height = 1.0))+
  print(dendro.plot, vp =  viewport(x = 0.77, y = 0.46, width = 0.18, height = 0.72))
dev.off()
```

## Merge with previous Techs
```{r}
#load data
EL4_prev_df <- read.csv("~/dataOS/CS_RNA/Pair_wise_comp/SS_NP_V_SC/crossed_gene_specific_info.csv")

#try all genes
# 2 for batch 1
cands <- table(Uni_gene) %>% as.data.frame(.)
res <- merge(EL4_prev_df,cands,by.x = "Gene_name",by.y = "Uni_gene")
write.table(res,"res_cand.csv",row.names = FALSE)
#%>% nrow(.)

# 3 for batch 2
cands_2 <- table(Uni_gene_2) %>% as.data.frame(.)
res2 <- merge(EL4_prev_df,cands_2,by.x = "Gene_name",by.y = "Uni_gene_2") 
write.table(res2,"res2.csv",row.names = FALSE)
#%>% nrow(.)
```


