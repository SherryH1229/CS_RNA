library(DESeq2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
# open files
L1_2_EX_raw <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/frankyan/OTHERS/lipid_RNA/Results/pipeOutput/mm_lipidRNA-L1_2/countsTable/mm_lipidRNA-L1_2.EX.counts.table",
                        sep = "\t")
L2_1_GB_raw <- read.csv("/mnt/extraids/OceanStor-SysCmn-5/frankyan/OTHERS/lipid_RNA/Results/pipeOutput/mm_lipidRNA-L1_2/countsTable/mm_lipidRNA-L1_2.GB.counts.table",
                        sep = "\t")

source("~/dataOS/CS_RNA/Functions.R")
geneCount_FPKM()

##---------------Norman:PBMC----------------------
PBMC_FPKM <- read.csv("/dataOS/nhhuang/csRNASequencing/PBMC_csRNA_Pulldown_Lysis/csRNAPulldown.EX.FPKM.csv",sep = ",")


PBMC_count <- read.csv("/dataOS/nhhuang/csRNASequencing/PBMC_csRNA_Pulldown_Lysis/PBMC_RNA_Pulldown_Counts.txt",
                       sep = "\t",header = FALSE) %>% .[,c(1,7:10)] 
colnames(PBMC_count) <- c("gene_id","BBL_1","BBL_Neg_1","BBL_2","BBL_Neg_2")
PBMC_count <- column_to_rownames(PBMC_count,"gene_id")
PBMC_count <- PBMC_count[c(3:nrow(PBMC_count)),] %>% format_samples(.,c("BBL_1","BBL_2","BBL_Neg_1","BBL_Neg_2"))
