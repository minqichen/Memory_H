# BiocManager::install('scater')
library(scater)
pwd<-getwd()
load('01_raw_cells');load('01_raw.counts');load('01_raw.tpm')
metadata.Ho <- raw_cells[raw_cells$Ho_Tx == "Ho",]
cells.idx.Ho <- match(rownames(metadata.Ho), colnames(raw.counts))
cells.idx.Ho<-cells.idx.Ho[!is.na(cells.idx.Ho)]
all.counts.Ho <- as.data.frame(raw.counts[,cells.idx.Ho])
dim(all.counts.Ho)

## create output directory
dir.create(paste0(pwd,"/output"))
dir.create(paste0(pwd,"/input"))
### Quality contol 
# Report total counts and expressed genes in each cell
SCQCstat <- function(data=data, is.expr=1){
  total_counts <- colSums(data, na.rm = T)
  total_genes <- apply(data, 2, function(x){
    sum(x>=is.expr)
  })
  return(list("total_counts"=total_counts, "total_genes"=unlist(total_genes)))
}


QC <- function(x,label){
  if(x=="Ho"){
    counts.input <- all.counts.Ho
    meta.input <- metadata.Ho
  } else {
    counts.input <- all.counts.Tx
    meta.input <- metadata.Tx
  }
  # quality control for cells and genes
  counts.qc <- SCQCstat(counts.input)
  size.drop <- scater::isOutlier(counts.qc$total_counts, nmads=3, type="both", log=F, batch = meta.input$phenotype)
  gene.drop <- scater::isOutlier(counts.qc$total_genes, nmads=3, type="both", log=F, batch = meta.input$phenotype)
  counts.filter <- counts.input[, !(size.drop | gene.drop)]
  counts.filter.qc <- SCQCstat(counts.filter)
  # filtering with gene number and umi counts
  gene.count <- 1000
  umi.count <- 20000
  counts.filter <- counts.filter[, counts.filter.qc$total_counts >= umi.count & counts.filter.qc$total_genes >= gene.count]
  dim(counts.filter) # Ho:23154 1270; Tx:23154  1059
  counts.filter.qc <- SCQCstat(counts.filter)
  # filtering low abundance genes
  ave.cut <- 0.2
  ave.counts <- rowMeans(counts.filter)
  counts.filter <- counts.filter[ log10(ave.counts) >= log10(ave.cut),]
  if(x=="Ho"){
    filter.counts.Ho <- counts.filter
    filter.tpm.Ho <- raw.tpm[row.names(counts.filter), colnames(counts.filter)]
    meta.Ho <- meta.input[colnames(counts.filter),]
    counts.Ho <- counts.input[,colnames(counts.filter)]
    tpm.Ho <- raw.tpm[,colnames(counts.filter)]
    save(filter.tpm.Ho, filter.tpm.Ho, counts.Ho, tpm.Ho, meta.Ho, file = paste0(pwd,"/input/01.Homeostasis.Cells.UMI_TPM_metadata.RData"))
    return(meta.Ho)
  } else {
    filter.counts.Tx <- counts.filter
    filter.tpm.Tx <- raw.tpm[row.names(counts.filter), colnames(counts.filter)]
    meta.Tx <- meta.input[colnames(counts.filter),]	 
    counts.Ho <- counts.input[,colnames(counts.filter)]
    tpm.Ho <- raw.tpm[,colnames(counts.filter)]
    save(filter.counts.Tx, filter.tpm.Tx, counts.Ho, tpm.Ho, meta.Tx,file = paste0(pwd,"/input/01.Transplantation.Cells.UMI_TPM_metadata.RData"))  
    return(meta.Tx)
  }
}
meta.Ho <- QC("Ho")
