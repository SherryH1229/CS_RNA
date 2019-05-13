library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)


#DataBases
cols <- c("GENENAME","SYMBOL")

#Hs
Biotype_df_Hs <- read.csv("~/dataOS/CS_RNA/surface_click/NK92/20190408/20190408_FPKMs.txt",sep="\t") %>% 
  as.data.frame(.) %>% .[,c(1,7)]
# RNA type
RNA_annotation <- as.character(unique(Biotype_df_Hs$Biotype))
RNA_mainType <- grep("RNA",RNA_annotation,value=TRUE)

colnames(Biotype_df_Hs)[1] <- "tracking_id"
anno_info_Hs <- AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(Biotype_df_Hs$tracking_id), columns=cols, keytype="ENSEMBL")
anno_info_Hs <- merge(anno_info_Hs,Biotype_df_Hs,by.x = "ENSEMBL",by.y ="tracking_id")
remove(Biotype_df_Hs)

#Mm
Biotype_df_Mm <- read.csv("/home/xcao3/membraneRNA/pooledResult201901/cufflink_out/fpkmsWithClickNoDup.txt",sep="\t") %>% 
  as.data.frame(.) %>% .[,c(1,3)]
anno_info_Mm <- AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(Biotype_df_Mm$tracking_id), columns=cols, keytype="ENSEMBL")
anno_info_Mm <- merge(anno_info_Mm,Biotype_df_Mm,by.x = "ENSEMBL",by.y ="tracking_id")
remove(Biotype_df_Mm)
remove(cols)


