library(ggplot2)
imap_20_fil<-readRDS('DEgene_1118.RDS')
imap_20_fil <- ScaleData(imap_20_fil, verbose = FALSE)

### T
name<-c('LT.HSC',"CD4T"  , "CD8T" ,"CLP","CMP", "Tcm.L",  "Tcm.S" )
thsc_only<-list()
for( i in c(1:7)){
MK_1<-FindMarkers(imap_20_fil,
                  ident.1 = 'T.HSC', 
                  ident.2 = name[i],
                  min.pct = 0.01,
                  only.pos = T)
MK_1$gene<-rownames(MK_1)
write.csv(MK_1,paste0(name[i],'_T.HSC.MK.csv'))
MK_1_list<-list(MK_1)
thsc_only<-c(thsc_only,MK_1_list)
}

merge_1<-thsc_only[[1]]$gene

for ( i in c(1:7)){
merge_1<-intersect(merge_1,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_1)
write.csv(merge_1,'ONLY_T.HSC.csv')
ggsave(filename = 'ONLY_T.HSC.PDF',width = 11,height = 10)


############         Tcm   #############################################################
name<-c('LT.HSC',"CLP","CMP", 
        # "Tcm.L",  "Tcm.S" 
        'CD4T','CD8T'
)
thsc_only<-list()
for( i in c(1:5)){
  MK_1<-FindMarkers(imap_20_fil,
                    ident.1 = c('Tcm.L','Tcm.S'), 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_1$gene<-rownames(MK_1)
  # write.csv(MK_1,paste0(name[i],'_Tcm.L.MK.csv'))
  MK_1_list<-list(MK_1)
  thsc_only<-c(thsc_only,MK_1_list)
}

merge_1<-thsc_only[[1]]$gene

for ( i in c(1:5)){
  merge_1<-intersect(merge_1,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_1)


name<-c('LT.HSC',"CLP","CMP", 
        # "Tcm.L",  "Tcm.S" 
        'CD4T','CD8T'
)
thsc_only2<-list()
for( i in c(1:5)){
  MK_2<-FindMarkers(imap_20_fil,
                    ident.1 = 'T.HSC', 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_2$gene<-rownames(MK_2)
  MK_2_list<-list(MK_2)
  thsc_only2<-c(thsc_only2,MK_2_list)
}

merge_2<-thsc_only2[[1]]$gene

for ( i in c(1:5)){
  merge_2<-intersect(merge_2,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_2)

merge<-intersect(merge_2,merge_1)
DoHeatmap(imap_20_fil, features =merge)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))

ggsave(filename = './1223/Tcm_T.HSC.HEATMAP1223.PDF',width = 11,height = 10)
write.csv(merge,'./1223/Tcm__T.HSC.HEATMAP1223.csv')

merge<-read.csv('./1223/Tcm__T.HSC.HEATMAP1223.csv')
DoHeatmap(imap_20_fil, features =merge$x)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))
ggsave(filename = './1223/Tcm_T.HSC.HEATMAP1223_1.PDF',width = 11,height = 10)

gene.df.trans<-bitr(merge, fromType = "SYMBOL", 
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego_bp <- enrichGO(
  gene  = gene.df.trans$ENTREZID,
  keyType = "ENTREZID", 
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
write.csv(ego_bp@result,file = './1223/Tcm_T.HSC.GO.CSV')

############  CD4/8  ################################################################################
library(clusterProfiler)
library(org.Mm.eg.db)
name<-c('LT.HSC',"CLP","CMP", 
        "Tcm.L",  "Tcm.S" 
        # 'CD4T','CD8T'
)
thsc_only<-list()
for( i in c(1:5)){
  MK_1<-FindMarkers(imap_20_fil,
                    ident.1 = c('CD4T','CD8T'), 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_1$gene<-rownames(MK_1)
  # write.csv(MK_1,paste0(name[i],'_Tcm.L.MK.csv'))
  MK_1_list<-list(MK_1)
  thsc_only<-c(thsc_only,MK_1_list)
}

merge_1<-thsc_only[[1]]$gene

for ( i in c(1:5)){
  merge_1<-intersect(merge_1,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_1)



thsc_only2<-list()
for( i in c(1:5)){
  MK_2<-FindMarkers(imap_20_fil,
                    ident.1 = 'T.HSC', 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    logfc.threshold = 0.3,
                    only.pos = T)
  MK_2$gene<-rownames(MK_2)
  MK_2_list<-list(MK_2)
  thsc_only2<-c(thsc_only2,MK_2_list)
}

merge_2<-thsc_only2[[1]]$gene

for ( i in c(1:5)){
  merge_2<-intersect(merge_2,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_2)

merge<-intersect(merge_2,merge_1)
DoHeatmap(imap_20_fil, features =merge)
ggsave(filename = './1223/CD4_8_T.HSC.HEATMAP1223.PDF',width = 11,height = 10)
write.csv(merge,'./1223/CD4_8_T.HSC.HEATMAP1223.csv')

merge<-read.csv('./1223/CD4_8_T.HSC.HEATMAP1223.csv')
DoHeatmap(imap_20_fil, features =merge$x)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))
ggsave(filename = './1223/CD4_8_T.HSC.HEATMAP1223_1.PDF',width = 11,height = 10)

gene.df.trans<-bitr(merge, fromType = "SYMBOL", 
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego_bp <- enrichGO(
  gene  = gene.df.trans$ENTREZID,
  keyType = "ENTREZID", 
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
write.csv(ego_bp@result,file = './1223/CD4_8_T.HSC.GO.CSV')


############  CLP  ################################################################################
library(clusterProfiler)
library(org.Mm.eg.db)
name<-c('LT.HSC',"CMP", 
        "Tcm.L",  "Tcm.S" ,
         'CD4T','CD8T'
)
thsc_only<-list()
for( i in c(1:6)){
  MK_1<-FindMarkers(imap_20_fil,
                    ident.1 = c('CLP'), 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_1$gene<-rownames(MK_1)
  MK_1_list<-list(MK_1)
  thsc_only<-c(thsc_only,MK_1_list)
}

merge_1<-thsc_only[[1]]$gene

for ( i in c(1:5)){
  merge_1<-intersect(merge_1,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_1)



thsc_only2<-list()
for( i in c(1:5)){
  MK_2<-FindMarkers(imap_20_fil,
                    ident.1 = 'T.HSC', 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_2$gene<-rownames(MK_2)
  MK_2_list<-list(MK_2)
  thsc_only2<-c(thsc_only2,MK_2_list)
}

merge_2<-thsc_only2[[1]]$gene

for ( i in c(1:5)){
  merge_2<-intersect(merge_2,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_2)

merge<-intersect(merge_2,merge_1)
DoHeatmap(imap_20_fil, features =merge)
ggsave(filename = './1223/CLP_T.HSC.HEATMAP1223.PDF',width = 11,height = 10)
write.csv(merge,'./1223/CLP_T.HSC.HEATMAP1223.csv')

merge<-read.csv('./1223/CLP_T.HSC.HEATMAP1223.csv')
DoHeatmap(imap_20_fil, features =merge$x)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))
ggsave(filename = './1223/CLP_T.HSC.HEATMAP1223_1.PDF',width = 11,height = 10)

gene.df.trans<-bitr(merge, fromType = "SYMBOL", 
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego_bp <- enrichGO(
  gene  = gene.df.trans$ENTREZID,
  keyType = "ENTREZID", 
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
write.csv(ego_bp@result,file = './1223/CLP_T.HSC.GO.CSV')


############  CMP  ################################################################################
if(F){
library(clusterProfiler)
library(org.Mm.eg.db)
name<-c('LT.HSC',"CLP", 
        "Tcm.L",  "Tcm.S" ,
        'CD4T','CD8T'
)
thsc_only<-list()
for( i in c(1:6)){
  MK_1<-FindMarkers(imap_20_fil,
                    ident.1 = c('CMP'), 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_1$gene<-rownames(MK_1)
  MK_1_list<-list(MK_1)
  thsc_only<-c(thsc_only,MK_1_list)
}

merge_1<-thsc_only[[1]]$gene

for ( i in c(1:5)){
  merge_1<-intersect(merge_1,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_1)



thsc_only2<-list()
for( i in c(1:5)){
  MK_2<-FindMarkers(imap_20_fil,
                    ident.1 = 'T.HSC', 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_2$gene<-rownames(MK_2)
  MK_2_list<-list(MK_2)
  thsc_only2<-c(thsc_only2,MK_2_list)
}

merge_2<-thsc_only2[[1]]$gene

for ( i in c(1:5)){
  merge_2<-intersect(merge_2,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_2)

merge<-intersect(merge_2,merge_1)
DoHeatmap(imap_20_fil, features =merge)
ggsave(filename = './1223/CMP_T.HSC.HEATMAP1223.PDF',width = 11,height = 10)
write.csv(merge,'./1223/CMP_T.HSC.HEATMAP1223.csv')

merge<-read.csv('./1223/CMP_T.HSC.HEATMAP1223.csv')
DoHeatmap(imap_20_fil, features =merge$x)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))
ggsave(filename = './1223/CMP_T.HSC.HEATMAP1223_1.PDF',width = 11,height = 10)

gene.df.trans<-bitr(merge, fromType = "SYMBOL", 
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego_bp <- enrichGO(
  gene  = gene.df.trans$ENTREZID,
  keyType = "ENTREZID", 
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
write.csv(ego_bp@result,file = './1223/CMP_T.HSC.GO.CSV')

}
############  LT.HSC  ################################################################################
if(F){
library(clusterProfiler)
library(org.Mm.eg.db)
name<-c('CMP',"CLP", 
        "Tcm.L",  "Tcm.S" ,
        'CD4T','CD8T'
)
thsc_only<-list()
for( i in c(1:6)){
  MK_1<-FindMarkers(imap_20_fil,
                    ident.1 = c('LT.HSC'), 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_1$gene<-rownames(MK_1)
  MK_1_list<-list(MK_1)
  thsc_only<-c(thsc_only,MK_1_list)
}

merge_1<-thsc_only[[1]]$gene

for ( i in c(1:5)){
  merge_1<-intersect(merge_1,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_1)



thsc_only2<-list()
for( i in c(1:5)){
  MK_2<-FindMarkers(imap_20_fil,
                    ident.1 = 'T.HSC', 
                    ident.2 = name[i],
                    min.pct = 0.1,
                    only.pos = T)
  MK_2$gene<-rownames(MK_2)
  MK_2_list<-list(MK_2)
  thsc_only2<-c(thsc_only2,MK_2_list)
}

merge_2<-thsc_only2[[1]]$gene

for ( i in c(1:5)){
  merge_2<-intersect(merge_2,thsc_only[[i]]$gene)
}

DoHeatmap(imap_20_fil, features =merge_2)

merge<-intersect(merge_2,merge_1)
DoHeatmap(imap_20_fil, features =merge)
ggsave(filename = './1223/LT.HSC_T.HSC.HEATMAP1223.PDF',width = 11,height = 10)
write.csv(merge,'./1223/LT.HSC_T.HSC.HEATMAP1223.csv')


merge<-read.csv('./1223/LT.HSC_T.HSC.HEATMAP1223.csv')
DoHeatmap(imap_20_fil, features =merge$x)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))
ggsave(filename = './1223/LT.HSC_T.HSC.HEATMAP1223_1.PDF',width = 11,height = 10)


gene.df.trans<-bitr(merge, fromType = "SYMBOL", 
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego_bp <- enrichGO(
  gene  = gene.df.trans$ENTREZID,
  keyType = "ENTREZID", 
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
write.csv(ego_bp@result,file = './1223/LT.HSC_T.HSC.GO.CSV')
}



AA<-read.csv('./1223/ONLY_T.HSC.csv')
gene.df.trans<-bitr(AA$x, fromType = "SYMBOL", 
                    toType = c("ENSEMBL","ENTREZID"),
                    OrgDb = org.Mm.eg.db)
ego_bp <- enrichGO(
  gene  = gene.df.trans$ENTREZID,
  keyType = "ENTREZID", 
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
write.csv(ego_bp@result,file = './1223/only_T.HSC.GO.CSV')



aa<-FindMarkers(imap_20_fil,
                ident.1 = 'T.HSC', 
                ident.2 = 'CMP',
                min.pct = 0.1,
                logfc.threshold = 0.3,
                only.pos = T)
DoHeatmap(imap_20_fil, features =rownames(aa) )


BB<-FindMarkers(imap_20_fil,
                ident.1 = 'T.HSC', 
                ident.2 = 'CLP',
                min.pct = 0.1,
                logfc.threshold = 0.3,
                only.pos = T)
DoHeatmap(imap_20_fil, features =rownames(BB) )

merge_2<-intersect(merge_1,rownames(BB))
DoHeatmap(imap_20_fil, features =merge_2)

merge<-read.csv('./1223/ONLY_T.HSC.csv')
DoHeatmap(imap_20_fil, features =merge$x)+scale_fill_gradientn(colors = c("#f6f6d8", "#ffffff", "#8a0000"))
ggsave(filename = './1223/ONLY_T.HSC.HEATMAP1223_1.PDF',width = 11,height = 10)
