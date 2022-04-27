library(SCENIC)
library(Seurat)
library(ggplot2)
scRNA<-readRDS('imap_20_fil_rn_dim30re0.5.RDS')
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
cellInfo <-  data.frame(scRNA@meta.data)

scenicOptions<-readRDS('int/scenicOptions.Rds')
##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
#调整部分二进制阈值
if(T){
  TZ<-c("Suz12 (15g)","Rcor1 (160g)","Pole3 (136g)","Ctcf_extended (2975g)","Myc (216g)",
        'Suz12_extended (1868g)','Pole3_extended (203g)')
AUCmatrix_pick<-cbind(AUCmatrix[,"Suz12 (15g)"],
                      AUCmatrix[,"Rcor1 (160g)"],
                      AUCmatrix[,"Pole3 (136g)"],
                      AUCmatrix[,"Ctcf_extended (2975g)"],
                      AUCmatrix[,"Myc (216g)"],
                      AUCmatrix[,"Suz12_extended (1868g)"],
                      AUCmatrix[,"Pole3_extended (203g)"])
colnames(AUCmatrix_pick)<-c("Suz12 (15g)","Rcor1 (160g)","Pole3 (136g)",
                            "Ctcf_extended (2975g)","Myc (216g)","Suz12_extended (1868g)",
                            'Pole3_extended (203g)')
rownames(AUCmatrix_pick)<-rownames(AUCmatrix)
AUCmatrix_pick[,"Suz12 (15g)"] [ AUCmatrix_pick[,"Suz12 (15g)"] >= 0.232] <- 1
AUCmatrix_pick[,"Suz12 (15g)"] [ AUCmatrix_pick[,"Suz12 (15g)"] <0.232] <- 0

AUCmatrix_pick[,"Rcor1 (160g)"] [ AUCmatrix_pick[,"Rcor1 (160g)"] >= 0.3] <- 1
AUCmatrix_pick[,"Rcor1 (160g)"] [ AUCmatrix_pick[,"Rcor1 (160g)"] <0.3] <- 0

AUCmatrix_pick[,"Pole3 (136g)"] [ AUCmatrix_pick[,"Pole3 (136g)"] >= 0.3] <- 1
AUCmatrix_pick[,"Pole3 (136g)"] [ AUCmatrix_pick[,"Pole3 (136g)"] <0.3] <- 0

AUCmatrix_pick[,"Ctcf_extended (2975g)"] [ AUCmatrix_pick[,"Ctcf_extended (2975g)"] >= 0.544] <- 1
AUCmatrix_pick[,"Ctcf_extended (2975g)"] [ AUCmatrix_pick[,"Ctcf_extended (2975g)"] < 0.544] <- 0

AUCmatrix_pick[,"Myc (216g)"] [ AUCmatrix_pick[,"Myc (216g)"] >= 0.316] <- 1
AUCmatrix_pick[,"Myc (216g)"] [ AUCmatrix_pick[,"Myc (216g)"] <0.316] <- 0

AUCmatrix_pick[,"Suz12_extended (1868g)"] [ AUCmatrix_pick[,"Suz12_extended (1868g)"] >= 0.39] <- 1
AUCmatrix_pick[,"Suz12_extended (1868g)"] [ AUCmatrix_pick[,"Suz12_extended (1868g)"] <0.39] <- 0

AUCmatrix_pick[,'Pole3_extended (203g)'] [ AUCmatrix_pick[,'Pole3_extended (203g)'] >= 0.29] <- 1
AUCmatrix_pick[,'Pole3_extended (203g)'] [ AUCmatrix_pick[,'Pole3_extended (203g)'] <0.29] <- 0
}
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
scRNAauc@assays$integrated <- NULL
# saveRDS(scRNAauc,'scRNAauc.rds')
# write.csv(t(AUCmatrix),'AUCmatrix.csv')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds") #getIntName(scenicOptions, "aucell_binary_full")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
BINmatrix<-BINmatrix[,-which(colnames(BINmatrix) %in% TZ)]

BINmatrix<-BINmatrix[order(rownames(BINmatrix)),]
AUCmatrix_pick<-AUCmatrix_pick[order(rownames(AUCmatrix_pick)),]

BINmatrix<-cbind(BINmatrix,AUCmatrix_pick)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
scRNAbin@assays$integrated <- NULL
# saveRDS(scRNAbin, 'scRNAbin.rds')
write.csv(t(BINmatrix),'BINmatrix_S_1207.csv')


# 更换展示的tsne为UMAP
if(F){
uu<-scRNA@reductions$umap
uu2<-uu@cell.embeddings
# uu3<-list(uu2)

NN<-readRDS('int/tSNE_AUC_50pcs_30perpl.Rds')
NN[["Y"]]<-uu2
saveRDS(NN,'int/umap.rds')
}
#使用shiny互动调整阈值
if(F){
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all,
                                tSNE_fileName ='int/umap.rds')
savedSelections <- shiny::runApp(aucellApp)

newThresholds <- aucellApp$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
}

######## PICKED Regulon
TT1<-c('E2f8_extended (298g)','E2f8 (185g)','Tfdp1_extended (442g)',
       'Tfdp1 (329g)','Tfdp2 (75g)','Tfdp2_extended (188g)',
       'Tal1 (118g)','Tal1_extended (158g)',"Myc (216g)",
       "Ctcf_extended (2975g)")#T
TT2<-c('Bcl11a_extended (13g)','Rad21_extended (3121g)','Bclaf1_extended (3562g)',
       'Bclaf1 (2512g)','Atf1_extended (2413g)','Ep300_extended (2527g)',
       'Ep300 (1680g)','Hcfc1 (2050g)','Hcfc1_extended (3586g)',
       'Bdp1_extended (1494g)','Nfe2_extended (1133g)','Sin3a_extended (2781g)',
       'Sin3a (1848g)','Etv6_extended (2694g)','Etv6 (1945g)',
       'Rest_extended (2046g)','Gata2_extended (1448g)',
       'Kdm5b_extended (2857g)','Pml_extended (3256g)',
       'Pml (2029g)','Phf8_extended (2587g)','Phf8 (1568g)',
       'Cux1_extended (325g)','Stat5a_extended (201g)','Rest (43g)',
       'Nr3c1_extended (2164g)','Nr3c1 (1506g)')#T&HSC
TT3<-c('Suz12 (15g)','Pole3 (136g)',"Rcor1 (160g)")#T&Tcm
TT4<-c('E2f1 (426g)','Htatip2 (70g)','Rb1_extended (177g)',
       'Rb1 (70g)','Gata1_extended (76g)','Gata1 (73g)',
       'Klf1_extended (60g)','Klf1 (41g)','Tfdp2_extended (188g)',
       'Tfdp2 (75g)','Nfyb (910g)','Tal1_extended (158g)','Tal1 (118g)',
       'E2f1_extended (2743g)')#T&CMP
TT5<-c('Pole3_extended (203g)','Pole3 (136g)')#T&CMP&Tcm
TT6<-c('Kdm5a (215g)','Mafk_extended (1000g)','Taf1_extended (2756g)',
       'Taf1 (2400g)','Irf3_extended (1334g)','Kdm5a_extended (2467g)',
       'Fli1_extended (1019g)','Elk4_extended (2707g)','Elk4 (1987g)',
       'Smarcc2_extended (1598g)')#T & HSC & Tcm
TT7<-c('Ctcf_extended (2975g)','Suz12_extended (1868g)',
       'Hdac6_extended (1320g)')#T & HSC & Tcm &CMP
PATH<-c('onlyT','T&HSC','T&Tcm','T&CMP','T&CMP&Tcm','T&HSC&Tcm','T&HSC&Tcm&CMP')
# UMAP *3
if(F){
TT1 <- gsub(' \\(','_',TT1);TT1 <- gsub('\\)','',TT1)
TT2 <- gsub(' \\(','_',TT2);TT2 <- gsub('\\)','',TT2)
TT3 <- gsub(' \\(','_',TT3);TT3 <- gsub('\\)','',TT3)
TT4 <- gsub(' \\(','_',TT4);TT4 <- gsub('\\)','',TT4)
TT5 <- gsub(' \\(','_',TT5);TT5 <- gsub('\\)','',TT5)
TT6 <- gsub(' \\(','_',TT6);TT6 <- gsub('\\)','',TT6)
TT7 <- gsub(' \\(','_',TT7);TT7 <- gsub('\\)','',TT7)

TT_list<-list(TT1,TT2,TT3,TT4,TT5,TT6,TT7)

ann_colors = c(LT.HSC ="#f0c7b2",CMP="#97dba5",CLP="#b4c5fa",CD4T="#595ee6",
               CD8T ="#5992e6",T.HSC ='#fcfc5b' ,Tcm.S ='#db7878',Tcm.L='#e33232')
p3 = DimPlot(scRNA, reduction = 'umap', group.by = "orig.ident", label=T,cols = ann_colors)
p3



for (i in PATH){
dir.create(i)
}
for (oo in c(3)) {
  for (ft in TT_list[[oo]]){

#FeaturePlot
p1 = FeaturePlot(scRNAauc, features=ft, label=T, reduction = 'umap')
p2 = FeaturePlot(scRNAbin, features=ft, label=T, reduction = 'umap')
p3 = DimPlot(scRNA, reduction = 'umap', group.by = "orig.ident", label=T,cols = ann_colors)
plotc = p1|p2|p3
plotc
ggplot2::ggsave(file=paste0(PATH[oo],'/',ft,'.pdf') , plotc, width=14 ,height=4)
  }
}
}