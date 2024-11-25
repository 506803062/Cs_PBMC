######################### All PBMC #########################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scRNAtoolVis)
library(reshape2)
library(ggpubr)
library(cowplot)

mycolour <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
              '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
              '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175')
mycolour2 <- c('#96b437','#da93ab','#a9c9ed', '#2f3c28', '#96b437',
         '#da93ab','#e58932', '#80598f', '#7e331f', '#3b855a',
         '#c0b286', '#a9c9ed', '#ec977f', '#848482', '#604628',
         '#d26034', '#a64c6b', '#dbd245', '#eba83b', '#5d5092',
         '#222222', '#f2f3f4')

pbmc <- readRDS("path/result/pbmc_harmony_singlet.rds")

## Data standardization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## PCA dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:50)
ElbowPlot(pbmc, ndims = 50)

### UMAP dimensionality reduction
pbmc <- pbmc %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.5) %>% identity()
pdf(file = "path/image/All_harmony_UMAP_cluster.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "umap", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

### TSNE dimensionality reduction
pbmc <- RunTSNE(pbmc, reduction = "harmony",dims = 1:50)
pdf(file = "path/image/All_harmony_TSNE_cluster.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "tsne", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
saveRDS(pbmc, file = "path/result/pbmc_harmony_UMAP_TSNE_r0.5.rds")

## Annotated cell cluster
pbmc <- readRDS("path/result/pbmc_harmony_UMAP_TSNE_r0.5.rds")
current.cluster.ids <- c(0:31)
new.cluster.ids <- c("Monocyte", "NK", "CD8+ T cell", "CD4+ T cell", "CD4+ T cell", "NK", "Naive B cell", "CD8+ T cell", "NK", 
                     "Memory B cell", "CD8+ T cell", "Monocyte", "Megakaryocyte", "CD4+ T cell", "Monocyte", "Treg", "NK", 
                     "CD4+ T cell", "Monocyte", "Monocyte", "Monocyte", "RBC", "Monocyte", "HSC", "Monocyte", "NK", "NK", "Plasma cell", 
                     "Naive B cell", "DC", "CD4+ T cell", "Monocyte")
pbmc@meta.data$celltype = plyr::mapvalues(x = pbmc@meta.data[,"RNA_snn_res.0.5"], from = current.cluster.ids, to = new.cluster.ids)
pbmc$celltype <- factor(pbmc$celltype,level = c ("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC"))

## Cell type show
pdf(file = "path/image/All_harmony_UMAP_celltype.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = "celltype", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_TSNE_celltype.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "tsne", group.by = "celltype", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_UMAP_group.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_TSNE_group.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "tsne", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_UMAP_level.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = "level", cols = c('#96b437','#da93ab','#a9c9ed')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/All_harmony_TSNE_level.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "tsne", group.by = "level", cols = c('#96b437','#da93ab','#a9c9ed')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
saveRDS(pbmc, file = "path/result/pbmc_harmony_UMAP_TSNE_celltype.rds")
##Show marker genes
features <- c("CD3D", "FOXP3", "CD4", "CD8A", "CD8B", "MS4A1", "CD27", "IGHD", "CD14", "S100A8", "FCGR3A", "NKG7", "KLRF1", "CST3", "LYZ", "FCER1A", "PPBP", "IGHG1", "MZB1","HBB", "CPA3")
pdf(file = "path/image/Plot_markergene_celltype.pdf", width = 8.5, height = 6)
jjDotPlot(pbmc, gene = features, id = "celltype", ytree = F)
dev.off()

## Find DEGs among cell types
Idents(pbmc) = "celltype"
diff.wilcox = FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(diff.wilcox, "path/result/diff_genes_wilcox_celltype.csv", quote = FALSE)
pdf(file = "path/image/diffgene_volcano_polar_celltype.pdf", width = 8, height = 8)
jjVolcano(diffData = diff.wilcox, log2FC.cutoff = 0.5, topGeneN = 3, tile.col = mycolour2, size = 2, polar = T)
dev.off()
pdf(file = "path/image/diffgene_top10_heatmap.pdf", width = 5, height = 2)
averageHeatmap(object = pbmc,
               markerGene = top10$gene,
               column_split = 1:12,
               border = T, annoCol = TRUE, myanCol = mycolour2[1:12], showRowNames = F)
dev.off()

## Cell percentage analysis

### Sample proportion in each cell type
df1 <- table(pbmc$celltype,pbmc$orig.ident) %>% melt()
colnames(df1) <- c("Cluster","Sample","Number")
df1$Cluster <- factor(df1$Cluster,level = c ("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC"))
df1$Sample <- factor(df1$Sample,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/celltype_sample_percent.pdf", width = 5, height = 3)
ggplot(data = df1, aes(x =Number, y = Cluster, fill =  Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycolour) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
dev.off()

### Group proportion in each cell type
df2 <- table(pbmc$celltype,pbmc$group) %>% melt()
colnames(df2) <- c("Cluster","Group","Number")
df2$Cluster <- factor(df2$Cluster,level = c ("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC"))
df2$Group <- factor(df2$Group,level = c ("H", "P"))
pdf(file = "path/image/celltype_group_percent.pdf", width = 4.5, height = 3)
ggplot(data = df2, aes(x =Number, y = Cluster, fill =  Group)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c('#2967a0', '#af2337')) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
dev.off()

### Level proportion in each cell type
df3 <- table(pbmc$celltype,pbmc$level) %>% melt()
colnames(df3) <- c("Cluster","Level","Number")
df3$Cluster <- factor(df3$Cluster,level = c ("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC"))
df3$Level<- factor(df3$Level,level = c ("Control", "Mild", "Severe"))
pdf(file = "path/image/celltype_level_percent.pdf", width = 4.5, height = 3)
ggplot(data = df3, aes(x =Number, y = Cluster, fill =  Level)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed')) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90))
dev.off()

### Different percentage of each major cell type in group
table(pbmc$orig.ident)
prop.table(table(Idents(pbmc)))
table(Idents(pbmc), pbmc$orig.ident)
Cellratio <- prop.table(table(Idents(pbmc), pbmc$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
sample <- c("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12")
group <- c(rep("H", 17), rep("P", 12))
samples <- data.frame(sample, group)
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']
cellper$group <- samples[rownames(cellper),'group']
pplist = list()
sce_groups = c("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC")

for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))
  colnames(cellper_) = c('sample','group','percent')
  cellper_$percent = as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("H", "P") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/group_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['Treg']],
          pplist[['CD4+ T cell']],
          pplist[['CD8+ T cell']],
          pplist[['Memory B cell']],
          pplist[['Naive B cell']],
          pplist[['Monocyte']],
          pplist[['NK']],
		  pplist[['DC']],
		  pplist[['Megakaryocyte']],
		  pplist[['Plasma cell']],
		  pplist[['RBC']],
		  pplist[['HSC']]
		  )
dev.off()

### Different percentage of each major cell type in level
table(pbmc$orig.ident)
prop.table(table(Idents(pbmc)))
table(Idents(pbmc), pbmc$orig.ident)
Cellratio <- prop.table(table(Idents(pbmc), pbmc$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
sample <- c("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12")
group <- c(rep("Control", 17), "Mild", "Mild", "Mild", "Severe", "Severe", "Severe", 
               "Severe", "Severe", "Mild", "Mild", "Severe", "Mild")
samples <- data.frame(sample, group)
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']
cellper$group <- samples[rownames(cellper),'group']
pplist = list()
sce_groups = c("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC")

for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))
  colnames(cellper_) = c('sample','group','percent')
  cellper_$percent = as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper = quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent,fill=group)) +
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/level_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['Treg']],
          pplist[['CD4+ T cell']],
          pplist[['CD8+ T cell']],
          pplist[['Memory B cell']],
          pplist[['Naive B cell']],
          pplist[['Monocyte']],
          pplist[['NK']],
		  pplist[['DC']],
		  pplist[['Megakaryocyte']],
		  pplist[['Plasma cell']],
		  pplist[['RBC']],
		  pplist[['HSC']]
		  )
dev.off()

######## calculate immune inhibitory geneset core ########
Genes = readLines("path/ImmuneInhibitory_selected.txt")
Genes = intersect(Genes,rownames(pbmc[["RNA"]]$data))
# P group
## 1.construct dataframe to plot
cells = rownames(pbmc@meta.data)[pbmc$group == "P"]
subData = as.matrix(t(pbmc[["RNA"]]$data[,cells]))
avg_matrix = data.frame(
  row.names = rownames(subData),
  immuneInhibitory = apply(subData[,Genes],1,mean)
)
avg_matrix$celltype = pbmc$celltype[rownames(avg_matrix)]
write.table(avg_matrix,"path/result/allImmunCell_InhibitScoreMatrix_onlyP.csv")
avg_matrix$celltype = factor(avg_matrix$celltype,levels=c("Treg", "CD4+ T cell", "CD8+ T cell", "Memory B cell", "Naive B cell", "Monocyte", "NK", "DC", "Megakaryocyte", "Plasma cell", "RBC", "HSC"))

## 2.plot
pdf(file = "path/image/Immuneinhibitory_allImmuneCell_P.pdf", width = 4, height = 6)
boxplot(immuneInhibitory ~ celltype, data = avg_matrix, 
    col=mycolour2,
        par(las="2"),
        cex.lab = 1.5,cex.axis=1.5,
        ylab="",
        xlab="Immune inhibitory score",
        outline = FALSE, horizontal = TRUE
)
axis(side=1, at=1:10, labels=FALSE)
dev.off()