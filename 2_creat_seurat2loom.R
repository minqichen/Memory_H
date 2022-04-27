library(loomR)
library(Seurat)

## 构建Seurat对象并导出
  rm(list=ls())
  pwd <- getwd()
  library(pheatmap)
  library(grid)
  
  ### load expression and meta data of homeostasis and transplantation cell
  load('./input/01.Homeostasis.Cells.UMI_TPM_metadata.RData')

  cell.type <- "clust"
  cell.subtype <- c('T.HSC','Tcm.S','Tcm.L','HSC','CMP','CLP','CD4T','CD8T') 
  
  input.meta<-subset(meta.Ho,meta.Ho$phenotype %in% cell.subtype)  
  

  input.mtx <- filter.tpm.Ho[, rownames(input.meta)] #行加基因名可筛选基因 如pc.topGenes
  # input.mtx[is.na(input.mtx)]<-0
  
  library(Seurat)
  library(uwot)
  input.umi <- counts.Ho[,rownames(input.meta)]

  nf=5000
  raw.obj <- CreateSeuratObject(input.umi, min.cells = 3,  project = "clust")
  
  CD4.nm <- c(rownames(raw.obj@meta.data)[grep("CD4",rownames(raw.obj@meta.data),invert=F)]) 
  Idents(raw.obj, cells = CD4.nm) <-'CD4T'
  
  raw.obj <- NormalizeData(object = raw.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  raw.obj <- ScaleData(raw.obj, display.progress = F)
  raw.obj <- FindVariableFeatures(object = raw.obj, mean.function = ExpMean,
                                  dispersion.function = LogVMR, x.low.cutoff = 0.25,
                                  x.high.cutoff = 3, y.cutoff = 0.5,
                                  nfeatures = nf)
  saveRDS(raw.obj,file = "220115raw.obj.RDS")
  
sdata <- readRDS(file = "220115raw.obj.RDS")
# seurat对象转换为loop文件
sdata.loom <- as.loom(x = sdata, filename = "sdata220115.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()

