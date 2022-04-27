library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(monocle3)
library(monocle)
library(devtools)
library(reticulate)
library(org.Mm.eg.db)
library(clusterProfiler)
library( MAST)
use_python("/home/anaconda3/bin/python",required = T)
py_config()
py_module_available("umap")
library(PythonInR)
imap_20_fil<-readRDS('imap_20_fil.RDS')

markers1 <- FindAllMarkers(imap_20_fil, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(unique(markers1$gene))
markers <- FindAllMarkers(imap_20_fil, only.pos = TRUE,  min.pct = 0.01, logfc.threshold =0.09)

gg<-c()
CC_a<-c()
for(i in c('CD4T','CD8T','ss','sl','LT.HSC','CLP','CMP')){
  aa<- read.csv(paste0('T.HSC_',i,'_DEs_sig.csv'))
  aaa<-as.character(aa$X)
  AA<-as.data.frame(aaa);AA$group<-i
  colnames(AA)<-c('gene','group')
  gg<-rbind(gg,AA)
  aa$group<-i
  CC_a<-rbind(CC_a,aa)
  
}

length(unique(gg$gene))
TT<-setdiff(markers1$gene,gg$gene)
TT_1<-setdiff(gg$gene,markers1$gene)

aaTT<-intersect(markers1$gene,gg$gene)  #交集

ONLY_mak<-markers1[which(markers1$gene %in% aaTT),]

markers <- FindAllMarkers(imap_20_fil, only.pos = TRUE,  min.pct = 0.01, logfc.threshold =0.09)

aTT<-intersect(markers$gene,TT_1) 
aTT_1<-setdiff(TT_1,markers$gene)

ONLY_mak_2 <-markers[which(markers$gene %in% aTT),]

ONLY_makk<-rbind(ONLY_mak,ONLY_mak_2)

length(unique(ONLY_makk$gene) )
ONLY_mak_s <-ONLY_makk[,c(7,6)]
colnames(ONLY_mak_s)<-c('gene','group')


aTT_1<-as.data.frame(aTT_1)
aTT_1$group<-'z'
colnames(aTT_1)<-c('gene','group')

end<-rbind(ONLY_mak_s,aTT_1)
end<-end[order(end$group),]
pdf('HEATMAP.PDF')
DoHeatmap(imap_20_fil, features = end$gene)   
dev.off()