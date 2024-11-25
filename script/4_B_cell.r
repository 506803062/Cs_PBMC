############ B cell subcluster ############

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scRNAtoolVis)
library(reshape2)
library(ggpubr)
library(cowplot)
library(clusterProfiler)
library(SCENIC)

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

# Extract B cell		 
pbmc <- readRDS("path/result/pbmc_harmony_UMAP_TSNE_celltype.rds")
Bcell <- subset(pbmc, celltype %in% c("Memory B cell", "Naive B cell"))
saveRDS(Bcell, file = "path/result/B_cell/pbmc_Bcell.rds")

################## Find DEGs in level ##################
diff_B_cell_CvsM <- FindMarkers(Bcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Mild")
diff_B_cell_CvsS <- FindMarkers(Bcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Severe")
write.csv(diff_B_cell_CvsM, "path/result/B_cell/Bcell_diff_genes_wilcox_Control_vs_Mild.csv", quote = FALSE)
write.csv(diff_B_cell_CvsS, "path/result/B_cell/Bcell_diff_genes_wilcox_Control_vs_Severe.csv", quote = FALSE)

## Draw Control_vs_Mild volcano plot
diff_B_cell_CvsM[which(diff_B_cell_CvsM$avg_log2FC  >= 0.25 & diff_B_cell_CvsM$p_val < 0.05),'sig'] <- 'up'
diff_B_cell_CvsM[which(diff_B_cell_CvsM$avg_log2FC  <= -0.25 & diff_B_cell_CvsM$p_val < 0.05),'sig'] <- 'down'
diff_B_cell_CvsM$symbol <- rownames(diff_B_cell_CvsM)
# Find top 10 up-regulated genes
up_data <- filter(diff_B_cell_CvsM, sig == 'up') %>% 
  distinct(symbol, .keep_all = TRUE) %>%              
  top_n(3, avg_log2FC)                        

# Find top 10 down-regulated genes
down_data <- filter(diff_B_cell_CvsM, sig == 'down') %>% 
  distinct(symbol, .keep_all = TRUE) %>%               
  top_n(-3, avg_log2FC)                              

pdf(file = "path/image/B_cell/Bcell_Control_vs_Mild_volcano.pdf", width = 4, height = 4)
ggplot(diff_B_cell_CvsM, aes(avg_log2FC , -log10(p_val), col = sig)) +
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
diff_B_cell_CvsS[which(diff_B_cell_CvsS$avg_log2FC  >= 0.25 & diff_B_cell_CvsS$p_val < 0.05),'sig'] <- 'up'
diff_B_cell_CvsS[which(diff_B_cell_CvsS$avg_log2FC  <= -0.25 & diff_B_cell_CvsS$p_val < 0.05),'sig'] <- 'down'
diff_B_cell_CvsS$symbol <- rownames(diff_B_cell_CvsS)
# Find top 10 up-regulated genes
up_data <- filter(diff_B_cell_CvsS, sig == 'up') %>% 
  distinct(symbol, .keep_all = TRUE) %>%              
  top_n(3, avg_log2FC)                          

# Find top 10 down-regulated genes
down_data <- filter(diff_B_cell_CvsS, sig == 'down') %>%
  distinct(symbol, .keep_all = TRUE) %>%                 
  top_n(-3, avg_log2FC)                             

pdf(file = "path/image/B_cell/Bcell_Control_vs_Severe_volcano.pdf", width = 4, height = 4)
ggplot(diff_B_cell_CvsS, aes(avg_log2FC , -log10(p_val), col = sig)) +
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

Idents(Bcell) = "level"
diff_Bcell_level  <- FindAllMarkers(Bcell, only.pos = TRUE,
                       min.pct = 0.25, 
                       logfc.threshold = 0.75)
diff_Bcell_level_sig  <- diff_Bcell_level[diff_Bcell_level$p_val_adj < 0.05, ]
write.csv(diff_Bcell_level, "path/result/B_cell/Bcell_diff_genes_wilcox_Level.csv", quote = FALSE)

group <- data.frame(gene=diff_Bcell_level_sig$gene,
                    group=diff_Bcell_level_sig$cluster)

Gene_ID <- bitr(diff_Bcell_level_sig$gene, fromType="SYMBOL", 
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

pdf(file = "path/image/B_cell/Bcell_level_GO.pdf", width = 3, height = 5.5)
dotplot(data_GO_sim, color = "p.adjust", size = "Count", showCategory=10, font.size = 5)
dev.off()

########### B cell subcluster ###########
Bcell$orig.ident <- factor(Bcell$orig.ident,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
Bcell <- FindVariableFeatures(Bcell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Bcell)
Bcell <- ScaleData(Bcell, features = all.genes)

##PCA dimensionality reduction
Bcell <- RunPCA(Bcell, features = VariableFeatures(object = Bcell))
VizDimLoadings(Bcell, dims = 1:2, reduction = "pca")
Bcell <- JackStraw(Bcell, num.replicate = 100, dims = 50)
Bcell <- ScoreJackStraw(Bcell, dims = 1:50)
p1 <- JackStrawPlot(Bcell, dims = 1:50)
p2 <- ElbowPlot(Bcell, ndims = 50)
pdf(file = "path/image/B_cell/ElbowPlot.pdf", width = 13, height = 3)
p1+p2
dev.off()

### UMAP dimensionality reduction
Bcell <- Bcell %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.9) %>% identity()
pdf(file = "path/image/B_cell/Bcell_UMAP_cluster_r0.9.pdf", width = 5.5, height = 4)
DimPlot(Bcell, reduction = "umap", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/B_cell/Bcell_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(Bcell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

### TSNE dimensionality reduction
Bcell <- RunTSNE(Bcell, reduction = "harmony",dims = 1:50)
pdf(file = "path/image/B_cell/Bcell_TSNE_cluster_r0.9.pdf", width = 5.5, height = 4)
DimPlot(Bcell, reduction = "tsne", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/B_cell/Bcell_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(Bcell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

######### Annotated cell cluster #########
current.cluster.ids <- c(0:15)
new.cluster.ids <- c("Naive B cell", "Naive B cell", "Naive B cell", "Naive B cell", "Memory B cell", "Follicular B cell", "Naive B cell", "Naive B cell", "Memory B cell", 
                     "Follicular B cell", "Naive B cell", "Memory B cell", "Memory B cell", "Naive B cell", "Megakaryocyte-like cell", "Naive B cell")
Bcell@meta.data$celltype_B = plyr::mapvalues(x = Bcell@meta.data[,"RNA_snn_res.1.2"], from = current.cluster.ids, to = new.cluster.ids)
Bcell$celltype_B <- factor(Bcell$celltype_B,level = c("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell"))

## Cell type show
pdf(file = "path/image/B_cell/Bcell_harmony_UMAP_celltype.pdf", width = 6, height = 4)
DimPlot(Bcell, reduction = "umap", group.by = "celltype_B", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/B_cell/Bcell_harmony_TSNE_celltype.pdf", width = 6, height = 4)
DimPlot(Bcell, reduction = "tsne", group.by = "celltype_B", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/B_cell/Bcell_harmony_UMAP_group.pdf", width = 5, height = 4)
DimPlot(Bcell, reduction = "umap", group.by = "group", cols = c('#af2337', '#2967a0')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/B_cell/Bcell_harmony_TSNE_group.pdf", width = 5, height = 4)
DimPlot(Bcell, reduction = "tsne", group.by = "group", cols = c('#af2337', '#2967a0')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/B_cell/Bcell_harmony_UMAP_level.pdf", width = 5, height = 4)
DimPlot(Bcell, reduction = "umap", group.by = "level", cols = c('#2967a0', '#af2337', '#ecc342')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/B_cell/Bcell_harmony_TSNE_level.pdf", width = 5, height = 4)
DimPlot(Bcell, reduction = "tsne", group.by = "level", cols = c('#2967a0', '#af2337', '#ecc342')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

Bcell$orig.ident <- factor(Bcell$orig.ident,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/B_cell/Bcell_harmony_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(Bcell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/B_cell/Bcell_harmony_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(Bcell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Show marker genes
features <- c("CD79A", "MS4A1", "CD79B", "CD22", "PPBP", "PF4", "IGHD", "IL4R", "TCL1A", "FCER2", "CD19", "CD27", "AIM2")
pdf(file = "path/image/B_cell/VlnPlot_markergene_celltype_B.pdf", width = 8, height = 2)
VlnPlot(Bcell, features = features, group.by = "celltype_B", stack = T, cols = mycolour)
dev.off()

## Find DEGs among cell types
Idents(Bcell) = "celltype_B"
diff.wilcox = FindAllMarkers(Bcell, min.pct = 0.25, logfc.threshold = 0.25)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(diff.wilcox, "path/result/B_cell/diff_genes_wilcox_celltype_B.csv", quote = FALSE)

pdf(file = "path/image/B_cell/diffgene_volcano_polar_celltype.pdf", width = 8, height = 8)
jjVolcano(diffData = diff.wilcox, log2FC.cutoff = 0.5, topGeneN = 3, tile.col = mycolour2, size = 2, polar = T)
dev.off()

pdf(file = "path/image/B_cell/diffgene_top10_celltype_T_heatmap.pdf", width = 5, height = 3)
averageHeatmap(object = Bcell,
               markerGene = top10$gene,
               column_split = 1:4,
               border = T, annoCol = TRUE, myanCol = mycolour2[1:4], showRowNames = F)
dev.off()

## Cell percentage analysis

### Sample proportion in each cell type
df1 <- table(Bcell$celltype_B,Bcell$orig.ident) %>% melt()
colnames(df1) <- c("Cluster","Sample","Number")
df1$Cluster <- factor(df1$Cluster,level = c ("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell"))
df1$Sample <- factor(df1$Sample,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/B_cell/Bcell_celltype_sample_percent.pdf", width = 5, height = 3)
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
df2 <- table(Bcell$celltype_B,Bcell$group) %>% melt()
colnames(df2) <- c("Cluster","Group","Number")
df2$Cluster <- factor(df2$Cluster,level = c ("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell"))
df2$Group <- factor(df2$Group,level = c ("H", "P"))
pdf(file = "path/image/B_cell/Bcell_celltype_group_percent.pdf", width = 4.5, height = 3)
ggplot(data = df2, aes(x =Number, y = Cluster, fill =  Group)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c('#af2337', '#2967a0')) +
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
df3 <- table(Bcell$celltype_B,Bcell$level) %>% melt()
colnames(df3) <- c("Cluster","Level","Number")
df3$Cluster <- factor(df3$Cluster,level = c ("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell"))
df3$Level<- factor(df3$Level,level = c ("Control", "Mild", "Severe"))
pdf(file = "path/image/B_cell/Bcell_celltype_level_percent.pdf", width = 4.5, height = 3)
ggplot(data = df3, aes(x =Number, y = Cluster, fill =  Level)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c('#af2337', '#ecc342', '#2967a0')) +
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
table(Bcell$orig.ident)
prop.table(table(Idents(Bcell)))
table(Idents(Bcell), Bcell$orig.ident)
Cellratio <- prop.table(table(Idents(Bcell), Bcell$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
library(reshape2)
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
sce_groups = c("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
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
    scale_fill_manual(values=c('#af2337', '#2967a0'))
  
  ### Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("H", "P") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/B_cell/Bcell_group_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['Follicular B cell']],
          pplist[['Megakaryocyte-like cell']],
          pplist[['Naive B cell']],
          pplist[['Memory B cell']]
		  )
dev.off()

### Different ercentage of each T cell subtype in level
table(Bcell$orig.ident)
prop.table(table(Idents(Bcell)))
table(Idents(Bcell), Bcell$orig.ident)
Cellratio <- prop.table(table(Idents(Bcell), Bcell$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
library(reshape2)
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
sce_groups = c("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell")

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
    scale_fill_manual(values=c('#2967a0', '#af2337', '#ecc342'))
  
  ### Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/B_cell/Bcell_level_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['Follicular B cell']],
          pplist[['Megakaryocyte-like cell']],
          pplist[['Naive B cell']],
          pplist[['Memory B cell']]
		  )
dev.off()

######## Find DEGs in Naive B cell 

B_naivecell <- subset(B_naivecell, celltype %in% c("Naive B cell"))
saveRDS(B_naivecell, file = "path/result/B_cell/pbmc_B_naive_cell.rds")

diff_B_naivecell_CvsM <- FindMarkers(B_naivecell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Mild")
diff_B_naivecell_CvsS <- FindMarkers(B_naivecell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Severe")
write.csv(diff_B_naivecell_CvsM, "path/result/B_cell/B_naivecell_diff_genes_wilcox_Control_vs_Mild.csv", quote = F)
write.csv(diff_B_naivecell_CvsS, "path/result/B_cell/B_naivecell_diff_genes_wilcox_Control_vs_Severe.csv", quote = F)

## Draw Control_vs_Mild volcano plot
diff_B_naivecell_CvsM[which(diff_B_naivecell_CvsM$avg_log2FC  >= 0.25 & diff_B_naivecell_CvsM$p_val < 0.05),'sig'] <- 'up'
diff_B_naivecell_CvsM[which(diff_B_naivecell_CvsM$avg_log2FC  <= -0.25 & diff_B_naivecell_CvsM$p_val < 0.05),'sig'] <- 'down'
diff_B_naivecell_CvsM$symbol <- rownames(diff_B_naivecell_CvsM)
# Find top 10 up-regulated genes
up_data <- filter(diff_B_naivecell_CvsM, sig == 'up') %>%  
  distinct(symbol, .keep_all = TRUE) %>%              
  top_n(3, avg_log2FC)                           

# Find top 10 down-regulated genes
down_data <- filter(diff_B_naivecell_CvsM, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                   
  top_n(-3, avg_log2FC)                               

pdf(file = "path/image/B_cell/B_naivecell_Control_vs_Mild_volcano.pdf", width = 4, height = 4)
ggplot(diff_B_naivecell_CvsM, aes(avg_log2FC , -log10(p_val), col = sig)) +
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
diff_B_naivecell_CvsS[which(diff_B_naivecell_CvsS$avg_log2FC  >= 0.25 & diff_B_naivecell_CvsS$p_val < 0.05),'sig'] <- 'up'
diff_B_naivecell_CvsS[which(diff_B_naivecell_CvsS$avg_log2FC  <= -0.25 & diff_B_naivecell_CvsS$p_val < 0.05),'sig'] <- 'down'
diff_B_naivecell_CvsS$symbol <- rownames(diff_B_naivecell_CvsS)
# Find top 10 up-regulated genes
up_data <- filter(diff_B_naivecell_CvsS, sig == 'up') %>%  
  distinct(symbol, .keep_all = TRUE) %>%               
  top_n(3, avg_log2FC)                           

# Find top 10 down-regulated genes
down_data <- filter(diff_B_naivecell_CvsS, sig == 'down') %>% 
  distinct(symbol, .keep_all = TRUE) %>%                  
  top_n(-3, avg_log2FC)                               

pdf(file = "path/image/B_cell/B_naivecell_Control_vs_Severe_volcano.pdf", width = 4, height = 4)
ggplot(diff_B_naivecell_CvsS, aes(avg_log2FC , -log10(p_val), col = sig)) +
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

######## Find DEGs in Memory B cell

B_Memorycell <- subset(Bcell, celltype %in% c("Memory B cell"))
saveRDS(B_Memorycell, file = "path/result/B_cell/pbmc_B_Memory_cell.rds")

diff_B_Memorycell_CvsM <- FindMarkers(B_Memorycell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Mild")
diff_B_Memorycell_CvsS <- FindMarkers(B_Memorycell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Severe")
write.csv(diff_B_Memorycell_CvsM, "path/result/B_cell/B_Memorycell_diff_genes_wilcox_Control_vs_Mild.csv", quote = F)
write.csv(diff_B_Memorycell_CvsS, "path/result/B_cell/B_Memorycell_diff_genes_wilcox_Control_vs_Severe.csv", quote = F)

## Draw Control_vs_Mild volcano plot
diff_B_Memorycell_CvsM[which(diff_B_Memorycell_CvsM$avg_log2FC  >= 0.25 & diff_B_Memorycell_CvsM$p_val < 0.05),'sig'] <- 'up'
diff_B_Memorycell_CvsM[which(diff_B_Memorycell_CvsM$avg_log2FC  <= -0.25 & diff_B_Memorycell_CvsM$p_val < 0.05),'sig'] <- 'down'
diff_B_Memorycell_CvsM$symbol <- rownames(diff_B_Memorycell_CvsM)
# Find top 10 up-regulated genes
up_data <- filter(diff_B_Memorycell_CvsM, sig == 'up') %>%  
  distinct(symbol, .keep_all = TRUE) %>%               
  top_n(3, avg_log2FC)                           

# Find top 10 down-regulated genes
down_data <- filter(diff_B_Memorycell_CvsM, sig == 'down') %>% 
  distinct(symbol, .keep_all = TRUE) %>%                   
  top_n(-3, avg_log2FC)                               
library(ggrepel)
pdf(file = "path/image/B_cell/B_Memorycell_Control_vs_Mild_volcano.pdf", width = 4, height = 4)
ggplot(diff_B_Memorycell_CvsM, aes(avg_log2FC , -log10(p_val), col = sig)) +
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
diff_B_Memorycell_CvsS[which(diff_B_Memorycell_CvsS$avg_log2FC  >= 0.25 & diff_B_Memorycell_CvsS$p_val < 0.05),'sig'] <- 'up'
diff_B_Memorycell_CvsS[which(diff_B_Memorycell_CvsS$avg_log2FC  <= -0.25 & diff_B_Memorycell_CvsS$p_val < 0.05),'sig'] <- 'down'
diff_B_Memorycell_CvsS$symbol <- rownames(diff_B_Memorycell_CvsS)
# Find top 10 down-regulated genes
up_data <- filter(diff_B_Memorycell_CvsS, sig == 'up') %>%  
  distinct(symbol, .keep_all = TRUE) %>%               
  top_n(3, avg_log2FC)                           

# Find top 10 down-regulated genes
down_data <- filter(diff_B_Memorycell_CvsS, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                   
  top_n(-3, avg_log2FC)                               

pdf(file = "path/image/B_cell/B_Memorycell_Control_vs_Severe_volcano.pdf", width = 4, height = 4)
ggplot(diff_B_Memorycell_CvsS, aes(avg_log2FC , -log10(p_val), col = sig)) +
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

##################### Analysis TF in B cell TF by SCENIC #####################

## Cell meta information
Bcell <- readRDS(file = "path/result/B_cell/Bcell_harmony_UMAP_TSNE_celltype.rds")
cellInfo <- data.frame(Bcell@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_B")] <-"celltype_B"
cellInfo <- cellInfo[,c("sample","cluster","celltype_B")]
saveRDS(cellInfo, file="path/result/B_cell/SCENIC/int/cellInfo.Rds")

## Expression matrix
exprMat <- as.matrix(Bcell@assays$RNA$counts)
mydbDIR <- "path/result/B_cell/SCENIC/cisTarget"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp","10kb")
scenicOptions <- initializeScenic(org="hgnc",nCores=50,dbDir=mydbDIR,dbs = mydbs,datasetTitle ="B_cell")
saveRDS(scenicOptions, "path/result/B_cell/SCENIC/int/scenicOptions.rds")

## Inference of transcriptional regulatory networks ##

## Gene filtering
genesKept <- geneFiltering(exprMat, scenicOptions,minCountsPerGene=3 * 0.01 * ncol(exprMat),minSamples=ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
## Computed correlation matrix
runCorrelation(exprMat_filtered, scenicOptions)

## TF-Targets correlation regression analysis
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)

## Extrapolate co-expression modules
runSCENIC_1_coexNetwork2modules(scenicOptions)

# Extrapolate transcriptional regulatory networks
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget"))

## regulon activity scoring and visualization ##
## regulons calculate AUC values and perform downstream analysis
exprMat_all<- as.matrix(Bcell@assays$RNA$counts)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

################

runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

## Import the original regulon AUC matrix
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(Bcell, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(Bcell, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select='celltype_B')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
# Draw a heat map using regulon raw AUC values
celltypecolor <- c('#af2337', '#ecc342', '#2967a0', '#2f3c28')
names(celltypecolor) <- c("Follicular B cell", "Megakaryocyte-like cell", "Naive B cell", "Memory B cell")
ann_colors <- list(celltype_B = celltypecolor)
pdf(file = "path/image/B_cell/Bcell_regulon_AUC.pdf", width = 8, height = 5)
pheatmap(AUCmatrix, show_colnames=F, annotation_col=celltype, annotation_colors = ann_colors)
dev.off()

# Draw heat maps using regulon binary AUC values
pdf(file = "path/image/B_cell/Bcell_regulon_Bin.pdf", width = 8, height = 5)
pheatmap(BINmatrix, show_colnames=F, annotation_col=celltype, annotation_colors = ann_colors,color = colorRampPalette(colors = c("white","black"))(100))
dev.off()