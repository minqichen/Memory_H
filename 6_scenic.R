library(SCENIC)
dir.create("./int")
###### scenicOptions 
if(T){
scenicOptions <- initializeScenic(org="mgi", dbDir="../cisTarget_databases", nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
}
# scenicOptions <- readRDS("../int/scenicOptions.Rds")

############ 导入细胞/基因表达矩阵

scRNA<-readRDS('DEgene_1118.RDS')
exprMat <- scRNA@assays$RNA@counts
cellInfo <- data.frame(scRNA@meta.data)


library(RcisTarget)
dbFilePath <- getDatabases(scenicOptions)[[1]]
motifRankings <- importRankings(dbFilePath)
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells<-rownames(exprMat)
length(genesLeft_minCells)
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)
genesKept <- genesLeft_minCells_inDatabases
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)




class(exprMat_filtered)
runCorrelation(as.matrix(exprMat_filtered), scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 


library(GENIE3)
runGenie3(as.matrix(exprMat_filtered), scenicOptions)
save(exprMat_filtered,scenicOptions,file = paste0("GENIE3_data.Rdata"))


scRNA<-readRDS('DEgene_1118.RDS')
exprMat <- scRNA@assays$RNA@counts
cellInfo <- data.frame(scRNA@meta.data)

logMat <- log2(exprMat+1)
dim(exprMat)



runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, as.matrix(logMat))
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=logMat)




