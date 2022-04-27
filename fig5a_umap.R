imap_20_fil<-readRDS('imap_20_fil_rn_dim30re0.5.RDS') 
p1 <- DimPlot(imap_20_fil, reduction = "umap", group.by = 'orig.ident', pt.size = 1,label = TRUE,cols =col )
p2 <- DimPlot(imap_20_fil, reduction = "umap", label = TRUE, pt.size = 1)

pdf(paste('./Dimplot/DimPlot_dim_',d,'_res_',r,'.pdf'),width = 4.8,height = 4)
p1;p2
dev.off()