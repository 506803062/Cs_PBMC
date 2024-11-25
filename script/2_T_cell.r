############ T cell subcluster ############

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scRNAtoolVis)
library(reshape2)
library(ggpubr)
library(cowplot)
library(clusterProfiler)
library(TCellSI)

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

# Extract T cell		 
pbmc <- readRDS("path/result/pbmc_harmony_UMAP_TSNE_celltype.rds")
Tcell <- subset(pbmc, celltype %in% c("Treg", "CD4+ T cell", "CD8+ T cell"))
saveRDS(Tcell, file = "path/result/T_cell/pbmc_Tcell.rds")

Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tcell)
Tcell <- ScaleData(Tcell, features = all.genes)

## PCA dimensionality reduction
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
VizDimLoadings(Tcell, dims = 1:2, reduction = "pca")
Tcell <- JackStraw(Tcell, num.replicate = 100, dims = 50)
Tcell <- ScoreJackStraw(Tcell, dims = 1:50)
p1 <- JackStrawPlot(Tcell, dims = 1:50)
p2 <- ElbowPlot(Tcell, ndims = 50)
pdf(file = "path/image/T_cell/ElbowPlot.pdf", width = 13, height = 3)
p1+p2
dev.off()

### UMAP dimensionality reduction
Tcell <- Tcell %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.9) %>% identity()
pdf(file = "path/image/T_cell/Tcell_UMAP_cluster_r0.9.pdf", width = 5.5, height = 4)
DimPlot(Tcell, reduction = "umap", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(Tcell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

### TSNE dimensionality reduction
Tcell <- RunTSNE(Tcell, reduction = "harmony",dims = 1:50)
pdf(file = "path/image/T_cell/Tcell_TSNE_cluster_r0.9.pdf", width = 5.5, height = 4)
DimPlot(Tcell, reduction = "tsne", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(Tcell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
saveRDS(Tcell, file = "path/image/T_cell/Tcell_harmony_UMAP_TSNE_r0.9.rds")

######### Annotated cell cluster #########
current.cluster.ids <- c(0:16)
new.cluster.ids <- c("CD4+ naive T cell", "CD8+ effector T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD4+ memory T cell", "CD8+ naive T cell", "NKT", "CD4+ naive T cell", "CD4+ naive T cell", 
                     "CD4+ regulatory T cell", "CD4+ effector T cell", "MAIT cell", "CD4+ naive T cell", "CD8+ naive T cell", "CD8+ effector T cell", "CD4+ naive T cell", "CD4+ naive T cell")
Tcell@meta.data$celltype_T = plyr::mapvalues(x = Tcell@meta.data[,"RNA_snn_res.0.9"], from = current.cluster.ids, to = new.cluster.ids)
Tcell$celltype_T <- factor(Tcell$celltype_T,level = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT"))

## Cell type show
pdf(file = "path/image/T_cell/Tcell_harmony_UMAP_celltype.pdf", width = 6, height = 4)
DimPlot(Tcell, reduction = "umap", group.by = "celltype_T", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_harmony_TSNE_celltype.pdf", width = 6, height = 4)
DimPlot(Tcell, reduction = "tsne", group.by = "celltype_T", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/T_cell/Tcell_harmony_UMAP_group.pdf", width = 5, height = 4)
DimPlot(Tcell, reduction = "umap", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_harmony_TSNE_group.pdf", width = 5, height = 4)
DimPlot(Tcell, reduction = "tsne", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/T_cell/Tcell_harmony_UMAP_level.pdf", width = 5, height = 4)
DimPlot(Tcell, reduction = "umap", group.by = "level", cols = c('#96b437','#da93ab','#a9c9ed')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_harmony_TSNE_level.pdf", width = 5, height = 4)
DimPlot(Tcell, reduction = "tsne", group.by = "level", cols = c('#96b437','#da93ab','#a9c9ed')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

Tcell$orig.ident <- factor(Tcell$orig.ident,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/T_cell/Tcell_harmony_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(Tcell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_harmony_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(Tcell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Show marker genes
features <- c("CD4", "CX3CR1", "NKG7", "CCR7", "LEF1", "SELL", "FOXP3", "RTKN2", "IKZF2", 
"IL7R", "LTB", "GPR183", "CD27", "CD44", "CCL5", "GZMA", "GZMK", "CXCR4", "FGFBP2", 
"GZMH", "GZMB", "KLRD1", "LRRN3", "SLC4A10", "RORA", "NCR3", "TRDC", "TRGC1")
pdf(file = "path/image/T_cell/Plot_markergene_celltype_T.pdf", width = 11, height = 8)
jjDotPlot(Tcell, gene = features, id = "celltype_T", ytree = F)
dev.off()

## Find DEGs among cell types
Idents(Tcell) = "celltype_T"
diff.wilcox = FindAllMarkers(Tcell, min.pct = 0.25, logfc.threshold = 0.25)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(diff.wilcox, "path/image/T_cell/diff_genes_wilcox_celltype_T.csv", quote = FALSE)
pdf(file = "path/image/T_cell/diffgene_volcano_celltype.pdf", width = 11, height = 4)
jjVolcano(diffData = diff.wilcox, log2FC.cutoff = 0.5, topGeneN = 3, tile.col = mycolour2, size = 2)
dev.off()
pdf(file = "path/image/T_cell/diffgene_top10_celltype_T_heatmap.pdf", width = 5, height = 3)
averageHeatmap(object = Tcell,
               markerGene = top10$gene,
               column_split = 1:9,
               border = T, annoCol = TRUE, myanCol = mycolour2[1:9], showRowNames = F)
dev.off()

## Cell percentage analysis

### Sample proportion in each cell type
df1 <- table(Tcell$celltype_T,Tcell$orig.ident) %>% melt()
colnames(df1) <- c("Cluster","Sample","Number")
df1$Cluster <- factor(df1$Cluster,level = c ("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT"))
df1$Sample <- factor(df1$Sample,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/T_cell/Tcell_celltype_sample_percent.pdf", width = 5, height = 3)
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
df2 <- table(Tcell$celltype_T,Tcell$group) %>% melt()
colnames(df2) <- c("Cluster","Group","Number")
df2$Cluster <- factor(df2$Cluster,level = c ("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT"))
df2$Group <- factor(df2$Group,level = c ("H", "P"))
pdf(file = "path/image/T_cell/Tcell_celltype_group_percent.pdf", width = 4.5, height = 3)
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
df3 <- table(Tcell$celltype_T,Tcell$level) %>% melt()
colnames(df3) <- c("Cluster","Level","Number")
df3$Cluster <- factor(df3$Cluster,level = c ("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT"))
df3$Level<- factor(df3$Level,level = c ("Control", "Mild", "Severe"))
pdf(file = "path/image/T_cell/Tcell_celltype_level_percent.pdf", width = 4.5, height = 3)
ggplot(data = df3, aes(x =Number, y = Cluster, fill =  Level)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed')) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
dev.off()

### Different ercentage of each T cell subtype in group
table(Tcell$orig.ident)
prop.table(table(Idents(Tcell)))
table(Idents(Tcell), Tcell$orig.ident)
Cellratio <- prop.table(table(Idents(Tcell), Tcell$orig.ident), margin = 2)
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
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

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

pdf(file = "path/image/T_cell/Tcell_group_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

### Different ercentage of each T cell subtype in level
table(Tcell$orig.ident)
prop.table(table(Idents(Tcell)))
table(Idents(Tcell), Tcell$orig.ident)
Cellratio <- prop.table(table(Idents(Tcell), Tcell$orig.ident), margin = 2)
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
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

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
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ###Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/T_cell/Tcell_level_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()
saveRDS(Tcell, file = "path/image/T_cell/Tcell_harmony_UMAP_TSNE_celltype.rds")

# Find DEGs in level
diff_T_cell_CvsM <- FindMarkers(Tcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Mild")
diff_T_cell_CvsS <- FindMarkers(Tcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Severe")
write.csv(diff_T_cell_CvsM, "path/image/T_cell/Tcell_diff_genes_wilcox_Control_vs_Mild.csv", quote = FALSE)
write.csv(diff_T_cell_CvsS, "path/image/T_cell/Tcell_diff_genes_wilcox_Control_vs_Severe.csv", quote = FALSE)

## Draw Control_vs_Mild volcano plot
diff_T_cell_CvsM[which(diff_T_cell_CvsM$avg_log2FC  >= 0.25 & diff_T_cell_CvsM$p_val < 0.05),'sig'] <- 'up'
diff_T_cell_CvsM[which(diff_T_cell_CvsM$avg_log2FC  <= -0.25 & diff_T_cell_CvsM$p_val < 0.05),'sig'] <- 'down'
diff_T_cell_CvsM$symbol <- rownames(diff_T_cell_CvsM)
# Find top 10 up-regulated genes
up_data <- filter(diff_T_cell_CvsM, sig == 'up') %>% 
  distinct(symbol, .keep_all = TRUE) %>%               
  top_n(3, avg_log2FC)                           

# Find top 10 down-regulated genes
down_data <- filter(diff_T_cell_CvsM, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                   
  top_n(-3, avg_log2FC)                              

pdf(file = "path/image/T_cell/Tcell_Control_vs_Mild_volcano.pdf", width = 4, height = 4)
ggplot(diff_T_cell_CvsM, aes(avg_log2FC , -log10(p_val), col = sig)) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_color_manual(values = c("#00314f", "#93292d")) +
  labs(x="log2(FoldChange)",y="-log10 (p_val)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.25, 0.25), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  geom_text_repel(data = up_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) +  
  geom_text_repel(data = down_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) 
dev.off()

## Draw Control_vs_Severe volcano plot
diff_T_cell_CvsS[which(diff_T_cell_CvsS$avg_log2FC  >= 0.25 & diff_T_cell_CvsS$p_val < 0.05),'sig'] <- 'up'
diff_T_cell_CvsS[which(diff_T_cell_CvsS$avg_log2FC  <= -0.25 & diff_T_cell_CvsS$p_val < 0.05),'sig'] <- 'down'
diff_T_cell_CvsS$symbol <- rownames(diff_T_cell_CvsS)
# Find top 10 up-regulated genes
up_data <- filter(diff_T_cell_CvsS, sig == 'up') %>%  
  distinct(symbol, .keep_all = TRUE) %>%           
  top_n(3, avg_log2FC)                         

# Find top 10 down-regulated genes
down_data <- filter(diff_T_cell_CvsS, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                  
  top_n(-3, avg_log2FC)                       

pdf(file = "path/image/T_cell/Tcell_Control_vs_Severe_volcano.pdf", width = 4, height = 4)
ggplot(diff_T_cell_CvsS, aes(avg_log2FC , -log10(p_val), col = sig)) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_color_manual(values = c("#00314f", "#93292d")) +
  labs(x="log2(FoldChange)",y="-log10 (p_val)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.25, 0.25), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  geom_text_repel(data = up_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) +
  geom_text_repel(data = down_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) 
dev.off()

# GO analysis

Idents(Tcell) = "level"
diff_Tcell_level  <- FindAllMarkers(Tcell, only.pos = TRUE,
                       min.pct = 0.25, 
                       logfc.threshold = 0.75)
diff_Tcell_level_sig  <- diff_Tcell_level[diff_Tcell_level$p_val_adj < 0.05, ]
write.csv(diff_Tcell_level, "path/image/T_cell/Tcell_diff_genes_wilcox_Level.csv", quote = FALSE)

group <- data.frame(gene=diff_Tcell_level_sig$gene,
                    group=diff_Tcell_level_sig$cluster)

Gene_ID <- bitr(diff_Tcell_level_sig$gene, fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)

pdf(file = "path/image/T_cell/Tcell_level_GO.pdf", width = 3, height = 6.5)
dotplot(data_GO_sim, color = "p.adjust", size = "Count", showCategory=10, font.size = 5)
dev.off()

############ Score the status of T cell subsets at different levels ############

ResultScores_TCSS <- TCSS_Calculate(data.frame(Tcell@assays$RNA$counts))
write.csv(ResultScores_TCSS, "path/image/T_cell/Tcell_ResultScores_TCSS.csv", quote = FALSE)

################# Level #####################

ResultScores_TCSS <- TCSS_Calculate(data.frame(Tcell@assays$RNA$counts))
ResultScores_TCSS_2 <- data.frame(t(ResultScores_TCSS))

###Add group information

ResultScores_TCSS_2$sample <- rownames(ResultScores_TCSS_2)
ResultScores_TCSS_2$group <- Tcell$level
ResultScores_TCSS_2$celltype_T <- Tcell$celltype_T

###### Quiescence ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Quiescence,fill=group)) +
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Quiescence') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Quiescence)
  compare_means(Quiescence ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Quiescence.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Regulating ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Regulating,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Regulating') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Regulating)
  compare_means(Regulating ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Regulating.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Proliferation ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Proliferation,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Proliferation') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Proliferation)
  compare_means(Proliferation ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Proliferation.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Helper ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Helper,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Helper') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Helper)
  compare_means(Helper ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Helper.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Cytotoxicity ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Cytotoxicity,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Cytotoxicity') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ###Inter-group t-test
  labely = max(cellper_$Cytotoxicity)
  compare_means(Cytotoxicity ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Cytotoxicity.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Progenitor_exhaustion ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Progenitor_exhaustion,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Progenitor_exhaustion') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Progenitor_exhaustion)
  compare_means(Progenitor_exhaustion ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Progenitor_exhaustion.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Terminal_exhaustion ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Terminal_exhaustion,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Terminal_exhaustion') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Terminal_exhaustion)
  compare_means(Terminal_exhaustion ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Terminal_exhaustion.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Senescence ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Senescence,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Senescence') +
    scale_fill_manual(values=c('#96b437','#da93ab','#a9c9ed'))
  
  ### Inter-group t-test
  labely = max(cellper_$Senescence)
  compare_means(Senescence ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_Senescence.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

################# Group #####################

### Add group information

ResultScores_TCSS_2$sample <- rownames(ResultScores_TCSS_2)
ResultScores_TCSS_2$group <- Tcell$group
ResultScores_TCSS_2$celltype_T <- Tcell$celltype_T

###### Quiescence ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Quiescence,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Quiescence') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Quiescence)
  compare_means(Quiescence ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Quiescence.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Regulating ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Regulating,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Regulating') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Regulating)
  compare_means(Regulating ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Regulating.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Proliferation ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Proliferation,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Proliferation') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Proliferation)
  compare_means(Proliferation ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Proliferation.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Helper ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Helper,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Helper') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Helper)
  compare_means(Helper ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Helper.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Cytotoxicity ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Cytotoxicity,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Cytotoxicity') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Cytotoxicity)
  compare_means(Cytotoxicity ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Cytotoxicity.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Progenitor_exhaustion ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Progenitor_exhaustion,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Progenitor_exhaustion') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Progenitor_exhaustion)
  compare_means(Progenitor_exhaustion ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Progenitor_exhaustion.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Terminal_exhaustion ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Terminal_exhaustion,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Terminal_exhaustion') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$Terminal_exhaustion)
  compare_means(Terminal_exhaustion ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Terminal_exhaustion.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()

###### Senescence ######

pplist = list()
sce_groups = c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell", "CD8+ effector memory T cell", "CD8+ effector T cell", 
                           "CD8+ naive T cell", "MAIT cell", "NKT")

for(group_ in sce_groups){
  cellper_  = subset(ResultScores_TCSS_2, celltype_T %in% group_)
  pp1 = ggplot(cellper_,aes(x=group,y=Senescence,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Senescence') +
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ###Inter-group t-test
  labely = max(cellper_$Senescence)
  compare_means(Senescence ~ group,  data = cellper_)
  my_comparisons <- list(c("H", "P"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}
pdf(file = "path/image/T_cell/Tcell_subtype_diff_group_Senescence.pdf", width = 10, height = 9)
plot_grid(pplist[['CD4+ effector T cell']],
          pplist[['CD4+ naive T cell']],
          pplist[['CD4+ regulatory T cell']],
          pplist[['CD4+ memory T cell']],
          pplist[['CD8+ effector memory T cell']],
          pplist[['CD8+ effector T cell']],
          pplist[['CD8+ naive T cell']],
		  pplist[['MAIT cell']],
		  pplist[['NKT']]
		  )
dev.off()