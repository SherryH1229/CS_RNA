library(DESeq2)
library(stringr)
library(ggplot2)
library(pasilla)
library(ggplot2)
library(dplyr)
library(topGO)

#fileNames_CS <- list.files(pattern = "*.bam")
#read Data
setwd("~/dataOS/CS_RNA/surface_clic k/featureCount/")

counts <- read.csv("count_matrix.txt",sep = "\t", head = TRUE,skip=1,row.names = "Geneid")

#format the data, remove unneccessary cells and trnasform the colnames
counts <- counts[ ,6:ncol(counts)]
sampleName <- str_split_fixed(colnames(counts),"\\.",7) 
sampleName <- gsub(".bam","",sampleName)
sampleName <- gsub("datapool.","",sampleName)[,7]

colnames(counts) <- sampleName

# Analysis with DESeq2---------------------------------------
counts_matrix <- as.matrix(counts) 
condition <- factor(c(rep("Surface",2), rep("Total", 3)),levels = c("Total","Surface"))

condition
#(condition <- factor(c(rep("Surface",3), rep("Total", 2))))
col_data <- data.frame(row.names = colnames(counts_matrix),condition)
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,colData = col_data,design=~condition)

#plot PCA
z <- plotPCA(rlog(dds), intgroup="condition")
nudge <- position_nudge(x=2,y = 4)
z + geom_label(aes(label = name), position = nudge)

setwd("~/dataOS/CS_RNA/surface_click/plots")

#pre-filtering out the gene with low count (count < 10)
dds_keep <-dds[rowSums(counts(dds)) >=10,]

#running the DESeq analysis
dds_keep<-DESeq(dds_keep,betaPrior=FALSE)
res_dds_keep<-results(dds_keep)
summary(res_dds_keep)
head(res_dds_keep)


#plot dispersion plot 
plotDispEsts(res_dds_keep, main="Dispersion plot")

#Plot p values
hist(res_dds_keep$pvalue)
p <- ggplot() + aes(res_dds_keep$pvalue)+ geom_histogram(binwidth=0.001, colour="black", fill="white")
p + labs(title="Distribution of P-values",
         x ="P_vals", y = "Count")
ggsave("P_value_Dist.pdf",p)
#filtering with certain threshold
candiates_names <- rownames (subset(res_dds_keep,log2FoldChange > 7,!(is.na(padj)))) %>% as.character(.)
lll <-length(candiates_names)
#plotCounts(dds_keep, gene=which.min(res_dds_keep$log2FoldChange), intgroup="condition")
#res_dds_keep[which.min(res_dds_keep$log2FoldChange),]
#---------output info table----------

sample_info <- read.csv("/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt",sep="\t") 

#RNA_mainType <- grep("RNA",as.character(unique(sample_info$Biotype)),value=TRUE)
Gene_type_info <- as.data.frame(sample_info[,c(1,2,3)])
#Gene_type_info_selected <
#resOrdered <-res.shr[order(-res.shr$log2FoldChange),]
all_gene_matrix <- as.data.frame(res_dds_keep) %>% .[,c("log2FoldChange","padj"),drop=FALSE] %>% 
  rownames_to_column(.,"gene_Name") %>% merge(Gene_type_info,.,by.x="gene_short_name",by.y = "gene_Name") %>% 
  .[order(-.$log2FoldChange),]

#colnames(all_gene_matrix)[1] <- "gene_short_name"
#candidates_matrix <- all_gene_matrix[which(all_gene_matrix$log2FoldChange>=9), ] %>% as.data.frame(.)
cand_list <- (subset(all_gene_matrix,log2FoldChange>=7 &!(is.na(padj)) & padj < 0.0025))$gene_short_name %>% 
  as.character(.)
#Gene_type_lfc_ordered<- merge (Gene_type_info,candidates_matrix,by="gene_short_name") %>% .[order(-.$log2FoldChange),]
length(cand_list)
#candiates_0.00025_geneList <- rownames (subset(Gene_type_pval_ordered, padj < 0.00025,!(is.na(padj)))) %>% as.character(.)
candiates <- all_gene_matrix[all_gene_matrix$gene_short_name %in% cand_list, ] %>% .[order(-.$log2FoldChange),]
setwd("~/dataOS/CS_RNA/surface_click/DEseq")
#write.csv(candidates_matrix,"candiates_l2fc<=-9.csv")
write.csv(candiates,"candiates_7_0.0025.csv")


david_res <- read.csv("DAVID.csv",sep = "|",head = FALSE) %>% as.data.frame(.)
david_res$V1 <- NULL
david_res$V5 <- NULL
colnames(david_res)<- c("tracking_id","name","species")

candiates$tracking_id <- as.character(candiates$tracking_id)
david_res$tracking_id <- gsub(" ", "", david_res$tracking_id)

out_res <- merge(candiates, david_res, by = "tracking_id")%>% .[order(-.$log2FoldChange),]
write.csv(out_res,"combined_info_update_2.csv")
#get the first occurance of each type or RNA and write as output csv file 
RNA.first <- all_gene_matrix[match(unique(all_gene_matrix$Biotype), all_gene_matrix$Biotype),]
write.csv(RNA.first,"RNA_first_Occur_l2fc.csv")

#plotCounts(dds_keep, gene="Olfr881", intgroup="condition")
#plot the top expressed genes 
candiates_names <- out_res$gene_short_name %>% as.character(.) %>% head(.,n = 4)
candiates_names
for (geneN in candiates_names){
  nameStr = paste(geneN, ".png")
  png(nameStr)
  plotCounts(dds_keep, gene=geneN, intgroup="condition")
  dev.off()
}


plotCounts(dds_keep, gene="Snord1c", intgroup="condition",transform = TRUE)

