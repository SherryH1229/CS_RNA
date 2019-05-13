library(dplyr)
library(tidyr)

setwd("~/dataOS/CS_RNA/Pair_wise_comp/SS_NP_V_SC/")
Genes_84 <- read.csv("crossed_gene_specific_info.csv")


setwd("~/dataOS/CS_RNA")
data1 <- read.csv("data1.csv")
data2 <- read.csv("data2.csv")
data3 <- read.csv("data3.csv")

data_all <- rbind(data1,data2,data3) 
data_all$X <- NULL

seqFISH_genes <- merge(data_all,Genes_84,by.x = "Genes", by.y = "Gene_name")
write.csv(seqFISH_genes,"seqFISH_genes_Info.csv")
