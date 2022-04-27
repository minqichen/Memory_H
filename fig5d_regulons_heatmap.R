library(pheatmap)
BINmatrix<-read.csv('BINmatrix_S_1207.csv')
rownames(BINmatrix)<-BINmatrix[,1];BINmatrix<-BINmatrix[,-1]
for  (i in rownames(BINmatrix)) {
  if (sum(BINmatrix[i,])< 3){
    BINmatrix[i,]<-0
  }
}
BINmatrix1<-BINmatrix[!duplicated(BINmatrix),]
input.mtx<-BINmatrix1

################
input.mtx1<-apply(input.mtx,2,as.numeric)
rownames(input.mtx1)<-rownames(input.mtx)
input.mtx<-input.mtx1
# 
# input.mtx1<-subset(input.mtx,input.mtx[,1]!='NA')
# OO<- input.mtx1[apply(input.mtx1, 1, function(x) sd(x)==0),]

input.mtx1 <- input.mtx1[apply(input.mtx1,1, function(x) sd(x)!=0),]
input.mtx<-input.mtx1

##################33333
input.dist.row <- as.dist(1-cor(t(input.mtx), method = "pearson"))
input.clust.row <- hclust(input.dist.row, method = "ward.D2")

input.dist.col <- as.dist(1-cor(input.mtx, method = "pearson"))
input.clust.col <- hclust(input.dist.col, method = "ward.D2")
col.panel <- colorRampPalette(colors = c("black", "gold"))

input.cutree.row <- cutree(input.clust.row, k = 1 )
input.cutree.col <- cutree(input.clust.col, k = 1)
# input.meta$clusters <- input.cutree.col

ANNO<-read.csv('cellInfo1112.csv')
rownames(ANNO)<-ANNO[,1]
ANNO1<-as.data.frame(ANNO[,2]) 
rownames(ANNO1)<-rownames(ANNO)
colnames(ANNO1)<-'CellType'


input.heatmap <- pheatmap(input.mtx, color=col.panel(100), cluster_rows = T,
                          cluster_cols = F, legend = T, show_colnames = F, show_rownames = F,
                          clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                          clustering_method = "ward.D2", 
                          annotation_col = data.frame( "Cell type" = ANNO1$CellType, row.names= row.names(ANNO1)),
                          annotation_row = data.frame("Gene clusters" = factor(input.cutree.row), row.names = as.character(rownames(input.mtx) )),
                          cutree_rows = max(input.cutree.row), cutree_cols = max(input.cutree.col))
dev.off()

ann_colors = list(
  Cell.type = c(LT.HSC ="#f0c7b2",CMP="#97dba5",CLP="#b4c5fa",CD4T="#595ee6",
                CD8T ="#5992e6",T.HSC ='#fcfc5b' ,Tcm.S ='#db7878',Tcm.L='#e33232'))
head(ann_colors)
pdf("BINmatrix.pdf", onefile = F)
pheatmap(input.mtx, color= colorRampPalette(c("#F6f6d8","#2c285e"))(100),
         cluster_rows = input.heatmap$tree_row,annotation_colors = ann_colors,
         # cluster_cols = input.heatmap$tree_col, legend =T, 
         show_colnames = F, show_rownames = T,cluster_cols = F,
         annotation_col = data.frame( "Cell type" = ANNO1$CellType, row.names= row.names(ANNO1)),
         annotation_row = data.frame("regulon clusters" = factor(input.cutree.row), 
                                     row.names = as.character(rownames(input.mtx) )),
         cutree_rows = max(input.cutree.row), cutree_cols = 8,fontsize_row=3,
         treeheight_row =20,
         gaps_col= c(46, 46+42,46+42+30,46+42+30+22,
                     46+42+30+22+191,46+42+30+22+191+49,
                     46+42+30+22+191+113+49,46+42+30+22+191+113+49+49),
         main = paste("BINmatrix"),scale="none")

dev.off()
