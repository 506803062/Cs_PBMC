############# Cell communication analysis #############

library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)

mycolour2 <- c('#af2337', '#ecc342', '#2967a0', '#2f3c28', '#96b437',
         '#da93ab','#e58932', '#80598f', '#7e331f', '#3b855a',
         '#c0b286', '#a9c9ed', '#ec977f', '#848482', '#604628',
         '#d26034', '#a64c6b', '#dbd245', '#eba83b', '#5d5092',
         '#222222', '#f2f3f4')
		 
######## All level #######

pbmc <- readRDS("path/result/pbmc_harmony_UMAP_TSNE_celltype.rds")
iTalk_data <- as.data.frame(t(pbmc@assays$RNA$counts+1))
iTalk_data$cell_type <- pbmc@meta.data$celltype
iTalk_data$compare_group <- pbmc@meta.data$level		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/iTALK/highly_exprs_genes.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/result/iTALK/ALL_LRpairs_Overview.csv")

pdf(file = "path/image/iTALK/All_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1.2,edge.label.cex = 0.9,vertex.size=30,arrow.width=3,edge.max.width=10,margin=0.2)
dev.off()
pdf(file = "path/image/iTALK/All_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

######## Control vs Mild ########
data2 <- subset(iTalk_data, subset=iTalk_data$compare_group=="Control"|iTalk_data$compare_group=="Mild")
deg_Treg<-DEG(data2 %>% filter(cell_type=="Treg"),method='DESeq2',contrast=c("Control","Mild"))
deg_CD4_T<-DEG(data2 %>% filter(cell_type=='CD4+ T cell'),method='DESeq2',contrast=c("Control","Mild"))
deg_CD8_T<-DEG(data2 %>% filter(cell_type=='CD8+ T cell'),method='DESeq2',contrast=c("Control","Mild"))
deg_Memory_B<-DEG(data2 %>% filter(cell_type=="Memory B cell"),method='DESeq2',contrast=c("Control","Mild"))
deg_Naive_B<-DEG(data2 %>% filter(cell_type=="Naive B cell"),method='DESeq2',contrast=c("Control","Mild"))
deg_Monocyte<-DEG(data2 %>% filter(cell_type=="Monocyte"),method='DESeq2',contrast=c("Control","Mild"))
deg_NK<-DEG(data2 %>% filter(cell_type=="NK"),method='DESeq2',contrast=c("Control","Mild"))
deg_DC<-DEG(data2 %>% filter(cell_type=="DC"),method='DESeq2',contrast=c("Control","Mild"))
deg_Megakaryocyte<-DEG(data2 %>% filter(cell_type=="Megakaryocyte"),method='DESeq2',contrast=c("Control","Mild"))
deg_Plasma_cell<-DEG(data2 %>% filter(cell_type=="Plasma cell"),method='DESeq2',contrast=c("Control","Mild"))
deg_RBC<-DEG(data2 %>% filter(cell_type=="RBC"),method='DESeq2',contrast=c("Control","Mild"))
deg_HSC<-DEG(data2 %>% filter(cell_type=="HSC"),method='DESeq2',contrast=c("Control","Mild"))
deg_all<-rbind(deg_Treg,deg_CD4_T,deg_CD8_T,deg_Memory_B,deg_Naive_B,deg_Monocyte,deg_NK,deg_DC,deg_Megakaryocyte,deg_Plasma_cell,deg_RBC,deg_HSC)
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(deg_all, datatype='DEG',comm_type=comm_type)
  res<-rbind(res,res_cat)
}
write.csv(res, "path/result/iTALK/Control_vs_Mild_res.csv")
res_2<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:50,]
write.csv(res_2, "path/result/iTALK/Control_vs_Mild_res_top50.csv", quote = FALSE)
dat <- read.csv("path/result/iTALK/Control_vs_Mild_res_top50.csv")
pdf(file = "path/image/iTALK/Control_vs_Mild_iTALK_LRPlot.pdf")
LRPlot(dat,datatype='DEG',cell_col=cell_col,link.arr.lwd=dat$cell_from_logFC,link.arr.width=dat$cell_to_logFC)
dev.off()

######## Control vs Severe ########
data3 <- subset(iTalk_data, subset=iTalk_data$compare_group=="Control"|iTalk_data$compare_group=="Severe")
deg_Treg<-DEG(data3 %>% filter(cell_type=="Treg"),method='DESeq2',contrast=c("Control","Severe"))
deg_CD4_T<-DEG(data3 %>% filter(cell_type=='CD4+ T cell'),method='DESeq2',contrast=c("Control","Severe"))
deg_CD8_T<-DEG(data3 %>% filter(cell_type=='CD8+ T cell'),method='DESeq2',contrast=c("Control","Severe"))
deg_Memory_B<-DEG(data3 %>% filter(cell_type=="Memory B cell"),method='DESeq2',contrast=c("Control","Severe"))
deg_Naive_B<-DEG(data3 %>% filter(cell_type=="Naive B cell"),method='DESeq2',contrast=c("Control","Severe"))
deg_Monocyte<-DEG(data3 %>% filter(cell_type=="Monocyte"),method='DESeq2',contrast=c("Control","Severe"))
deg_NK<-DEG(data3 %>% filter(cell_type=="NK"),method='DESeq2',contrast=c("Control","Severe"))
deg_DC<-DEG(data3 %>% filter(cell_type=="DC"),method='DESeq2',contrast=c("Control","Severe"))
deg_Megakaryocyte<-DEG(data3 %>% filter(cell_type=="Megakaryocyte"),method='DESeq2',contrast=c("Control","Severe"))
deg_Plasma_cell<-DEG(data3 %>% filter(cell_type=="Plasma cell"),method='DESeq2',contrast=c("Control","Severe"))
deg_RBC<-DEG(data2 %>% filter(cell_type=="RBC"),method='DESeq2',contrast=c("Control","Severe"))
deg_HSC<-DEG(data2 %>% filter(cell_type=="HSC"),method='DESeq2',contrast=c("Control","Severe"))
deg_all<-rbind(deg_Treg,deg_CD4_T,deg_CD8_T,deg_Memory_B,deg_Naive_B,deg_Monocyte,deg_NK,deg_DC,deg_Megakaryocyte,deg_Plasma_cell,deg_RBC,deg_HSC)
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(deg_all, datatype='DEG',comm_type=comm_type)
  res<-rbind(res,res_cat)
}
res_2<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:50,]
write.csv(res_2, "path/result/iTALK/Control_vs_Severe_res_top50.csv", quote = FALSE)
dat <- read.csv("path/result/iTALK/Control_vs_Severe_res_top50.csv")
pdf(file = "path/image/iTALK/Control_vs_Severe_iTALK_LRPlot.pdf")
LRPlot(dat,datatype='DEG',link.arr.lwd=dat$cell_from_logFC,link.arr.width=dat$cell_to_logFC)
dev.off()

##################### Control #####################

Control <- subset(pbmc, level %in% c("Control"))
iTalk_data <- as.data.frame(t(Control@assays$RNA$counts+1))
iTalk_data$cell_type <- Control@meta.data$celltype
iTalk_data$compare_group <- Control@meta.data$level		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/result/iTALK/highly_exprs_genes_Control.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/result/iTALK/ALL_LRpairs_Overview_Control.csv")
   
pdf(file = "path/image/iTALK/Control_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/image/iTALK/Control_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

##################### Mild #####################

Mild <- subset(pbmc, level %in% c("Mild"))
iTalk_data <- as.data.frame(t(Mild@assays$RNA$counts+1))
iTalk_data$cell_type <- Mild@meta.data$celltype
iTalk_data$compare_group <- Mild@meta.data$level		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/result/iTALK/highly_exprs_genes_Mild.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/result/iTALK/ALL_LRpairs_Overview_Mild.csv")
 
pdf(file = "path/image/iTALK/Mild_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/image/iTALK/Mild_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

##################### Severe #####################

Severe <- subset(pbmc, level %in% c("Severe"))
iTalk_data <- as.data.frame(t(Severe@assays$RNA$counts+1))
iTalk_data$cell_type <- Severe@meta.data$celltype
iTalk_data$compare_group <- Severe@meta.data$level		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/result/iTALK/highly_exprs_genes_Severe.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/result/iTALK/ALL_LRpairs_Overview_Severe.csv")
   
pdf(file = "path/image/iTALK/Severe_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/image/iTALK/Severe_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

##################### Control and Mild #####################

Control_Mild <- subset(pbmc, level %in% c("Control", "Mild"))
iTalk_data <- as.data.frame(t(ControlMild@assays$RNA$counts+1))
iTalk_data$cell_type <- Control_Mild@meta.data$celltype
iTalk_data$compare_group <- Control_Mild@meta.data$level		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/result/iTALK/highly_exprs_genes_Control_Mild.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/result/iTALK/ALL_LRpairs_Overview_Control_Mild.csv")
   
pdf(file = "path/image/iTALK/Control_Mild_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/image/iTALK/Control_Mild_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()