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

imap_20k<-read.csv('output11.csv',header = F)
imap_20k<-as.data.frame(t(imap_20k))
rownames(imap_20k)<-imap_20k[,1]
imap_20k<-imap_20k[,-1]
nn<-as.matrix( imap_20k[1,])
colnames(imap_20k)<-nn;imap_20k<-imap_20k[-1,]
imap_20k<-imap_20k[,order(colnames(imap_20k))]
colnames(imap_20k)
colnames(imap_20k)[332:444]<-sprintf('Tcm.L_%04d',c(1:113))
colnames(imap_20k)[445:493]<-sprintf('Tcm.S_%04d',c(1:49))

ANNO<-read.csv('cellInfo1112.csv')
rownames(ANNO)<-ANNO[,1]
ANNO1<-as.data.frame(ANNO[,2]) 
rownames(ANNO1)<-rownames(ANNO)
colnames(ANNO1)<-'CellType'

library(pheatmap)
library(ggplot2)

imap_20k_1<-apply(imap_20k, 2, as.numeric)
imap_20k_1<-as.data.frame(imap_20k_1)
rownames(imap_20k_1)<-rownames(imap_20k)



###############################
col<-c(HSC ="#f0c7b2",CMP="#97dba5",CLP="#b4c5fa",CD4T="#595ee6",
       CD8T ="#5992e6",T.HSC ='#faeb1b' ,Tcm.S ='#edc2b2',Tcm.L='#e33232')

imap_20 <- CreateSeuratObject(counts = imap_20k, project = "imap", min.cells = 3)


imap_20[["percent.mt"]] <- PercentageFeatureSet(object = imap_20, pattern = "^mt")
VlnPlot(object = imap_20,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.5)
imap_20_fil<-imap_20
# filter data
imap_20_fil <- subset(imap_20, subset = nFeature_RNA >5000 & nFeature_RNA < 10000 & percent.mt < 5)
imap_20_fil
saveRDS(imap_20,'imap_20.RDS')
saveRDS(imap_20_fil,'imap_20_fil.RDS')


pdf('./Dimplot/VlnPlot_QC.pdf',width = 20,height = 4)
VlnPlot(object = imap_20_fil,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,pt.size = 0.5,group.by = 'orig.ident',cols =col)
dev.off()
# Normalizing the data

imap_20_fil<-NormalizeData(imap_20_fil, normalization.method = "LogNormalize", scale.factor = 10000)


# Identification of highly variable features (feature selection)
imap_20_fil <- FindVariableFeatures(imap_20_fil, selection.method = "vst", nfeatures = 5000)

imap_20_fil <- ScaleData(imap_20_fil, verbose = FALSE)
imap_20_fil <- RunPCA(imap_20_fil, npcs = 30, verbose = FALSE)
  
#################3
imap_20_fil <- JackStraw(imap_20_fil, num.replicate = 100)
imap_20_fil <- ScoreJackStraw(imap_20_fil, dims = 1:20)



r=0.5
d=30

imap_20_fil <- RunUMAP(imap_20_fil, reduction = "pca", dims = 1:d)
imap_20_fil <- FindNeighbors(imap_20_fil, reduction = "pca", dims = 1:d)
imap_20_fil <- FindClusters(imap_20_fil, resolution = r)  

p1 <- DimPlot(imap_20_fil, reduction = "umap", group.by = 'orig.ident', pt.size = 1,label = TRUE,cols =col )
p2 <- DimPlot(imap_20_fil, reduction = "umap", label = TRUE, pt.size = 1)

pdf(paste('./Dimplot/DimPlot_dim_',d,'_res_',r,'.pdf'),width = 4.8,height = 4)
p1;p2
dev.off()


write.csv( data.frame(imap_20_fil@meta.data),file = "cellInfo1112.csv")
saveRDS(imap_20_fil,'imap_20_fil_rn_dim30re0.5.RDS')
write.table(imap_20_fil@reductions[["umap"]]@cell.embeddings,'imap_umap_position_220329.txt',sep = '\t')
write.table(as.matrix(imap_20_fil@assays$RNA@counts),file ='imap_count_220329.txt',sep = '\t')

### seurat for scenic
gg<-c()
for(i in c('CD4T','CD8T','ss','sl','LT.HSC','CLP','CMP')){
  aa<- read.csv(paste0('T.HSC_',i,'_DEs_sig.csv'))
  aaa<-as.character(aa$X)
  gg<-c(gg,aaa)
}
gg<-unique(gg)
gg<-as.data.frame(gg)
colnames(gg)<-'NAME'


imap_20k$NAME<-rownames(imap_20k)
imap.mark_all_U_m<-merge(gg,imap_20k,by='NAME')
dup<-data.frame(table(imap.mark_all_U_m$NAME))
dup<-subset(dup,dup$Freq>1)

rownames(imap.mark_all_U_m)<-imap.mark_all_U_m$NAME
imap.mark_all_U_m<-imap.mark_all_U_m[,-1]
imap.mark_all_U_m_s<-CreateSeuratObject(counts = imap.mark_all_U_m, project = "imap", min.cells = 3)

write.csv(t(as.matrix(imap.mark_all_U_m_s@assays$RNA@counts)),
          file = "DEgene_1118.csv")
write.csv( data.frame(imap.mark_all_U_m_s@meta.data),
           file = "cellInfo_1118.csv")
saveRDS(imap.mark_all_U_m_s,file='scenic/int/DEgene_1118.RDS')
