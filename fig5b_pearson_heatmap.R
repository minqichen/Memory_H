
imap_20_fil<-readRDS('imap_20_fil_rn_dim30re0.5.RDS')
va_gene<-imap_20_fil@assays[["RNA"]]@var.features

av <-AverageExpression(imap_20_fil,
                       group.by = "orig.ident",
                       assays = "RNA")
av=av[[1]]
aa<-subset(av,rownames(av) %in% va_gene)

av<-aa
av<-av[,c("T.HSC",'LT.HSC','Tcm.L','CMP','Tcm.S','CLP','CD4T','CD8T')]
cg=names(tail(sort(apply(av, 1, sd)),5000))
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),cluster_rows = T,cluster_cols = T)

pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),
                   cluster_rows = F,cluster_cols = F,display_numbers = TRUE,
                   legend_breaks = c(0,0.5,1),legend_labels = c("0",'0.5',"1")) 

bk <- c(seq(0,1,by=0.02))

pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),
                   cluster_rows = F,cluster_cols = F,display_numbers = TRUE,
         color = c(colorRampPalette(colors = c("#0f00b3",'#6d63db',"white"))(length(bk)/2), #down
                   colorRampPalette(colors = c("white",'#f5a1b2',"#fa5778"))(length(bk)/2)), #up
         legend_breaks=seq(0,1,by=0.2),
         breaks=bk,cellwidth = 30, cellheight = 30)
library(ggplot2)
ggsave(aa,filename = "pearson_S_3.pdf",width = 7,height = 7)


