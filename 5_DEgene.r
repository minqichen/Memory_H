library(Seurat)
library(data.table)
library(SingleCellExperiment)
library(DEsingle)
imap_20_fil.umap <- readRDS("imap_20_fil.RDS")


for(i in levels(imap_20_fil.umap@active.ident)){
  A_clust.name <- c(rownames(imap_20_fil.umap@meta.data)[grep('^T.HSC',imap_20_fil.umap$orig.ident)])
  C_clust.name<-c(rownames(imap_20_fil.umap@meta.data)[grep(paste0('^',i),imap_20_fil.umap$orig.ident,
                invert=F)])
  AC_sc_<- subset(imap_20_fil.umap,cells=c(A_clust.name,C_clust.name))
  AC_annda_ <-as.SingleCellExperiment(AC_sc_)

  
  group <- factor(c(rep("A",length(A_clust.name)), rep("C",length(C_clust.name))))
  
  # Detecting the DE genes with SingleCellExperiment input sce
  de_result<- DEsingle(AC_annda_, group)
  saveRDS(de_result,file= paste('T.HSC_',i,'_DEs_DEa.rds',sep=''))
  # Dividing the DE genes into 3 categories at threshold of FDR < 0.001
  results.classified <- DEtype(results = de_result, threshold = 0.01)
  results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.01, ]
  
  # Extract three types of DE genes separately
  results.DEs <- results.sig[results.sig$Type == "DEs", ]
  results.DEa <- results.sig[results.sig$Type == "DEa", ]
  results.DEg <- results.sig[results.sig$Type == "DEg", ]

  write.csv(results.sig,file= paste('T.HSC_',i,'_DEs_sig.csv',sep=''))
  write.csv(results.DEs,file= paste('T.HSC_',i,'_DEs_DEs.csv',sep=''))
  
  write.csv(results.DEa,file= paste('T.HSC_',i,'_DEs_DEa.csv',sep=''))
  write.csv(results.DEg,file= paste('T.HSC_',i,'_DEs_DEg.csv',sep=''))
}