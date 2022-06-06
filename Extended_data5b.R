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
                    only.pos = F)
  MK_1<-subset(MK_1,MK_1$avg_log2FC < -0 )
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
                    only.pos = F)
  MK_2<-subset(MK_2,MK_2$avg_log2FC < - 0)
  MK_2<-subset(MK_2,MK_2$p_val <  0.5)
  
  MK_2$gene<-rownames(MK_2)
  MK_2_list<-list(MK_2)
  thsc_only2<-c(thsc_only2,MK_2_list)
}
