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

