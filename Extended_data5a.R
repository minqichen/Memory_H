library(ggsci)
library(ggpubr)
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
imap_20_fil<-readRDS('imap_20_fil.rds')
source('violin_gene_exp_function_0.2_jitter.R')


fn<-c('CD4_8','CLP','CMP','Tcm','LT.HSC')

for (i in 1:5) {
   # he<-read.csv(paste0('F_',fn[i],'_T.HSC.HEATMAP1223.csv'))
  he<-read.csv(paste0(fn[i],'_T.HSC.HEATMAP1223.csv'))
  gene <- list( he$x)
  
  imap_20_fil <- AddModuleScore(  #对每个模块的基因集打分
    object = imap_20_fil,
    features = gene,
    ctrl = 100, #默认值是100
    name = fn[i])
}

he<-read.csv('ONLY_T.HSC.csv')
# he<-read.csv('F_ONLY_T.HSC.csv')
gene <- list( he$x)
imap_20_fil <- AddModuleScore(  #对每个模块的基因集打分
  object = imap_20_fil,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'noly_T')

colnames(imap_20_fil@meta.data)
  
    name<-as.data.frame( imap_20_fil@assays$RNA@counts)
    kk_cluster <-(as.matrix( imap_20_fil@meta.data$orig.ident)) # kk_cluster是细胞分组信息
    
    pdf('average_group_vio.pdf',width = 12,height = 5)
    for (i in 5:10){
      
      kk<-as.data.frame(t(imap_20_fil@meta.data[i]) )   # kk是表达矩阵
      colnames(kk)<-colnames(name)
      rownames(kk)<-colnames(imap_20_fil@meta.data)[i]
      violin_gene_exp(colnames(imap_20_fil@meta.data)[i],kk,kk_cluster)
    }
    
    dev.off()
    