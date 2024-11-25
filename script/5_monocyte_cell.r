############ Monocyte cell subcluster ############

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scRNAtoolVis)
library(reshape2)
library(ggpubr)
library(cowplot)
library(qusage)

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
		 
# Extract monocyte cell		 
pbmc <- readRDS("path/result/pbmc_harmony_UMAP_TSNE_celltype.rds")
Monocytecell <- subset(pbmc, celltype %in% c("Monocyte"))
saveRDS(Monocytecell, file = "path/result/Monocyte_cell/pbmc_Monocyte.rds")

########### Monocyte cell subcluster ###########

Monocytecell <- FindVariableFeatures(Monocytecell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Monocytecell)
Monocytecell <- ScaleData(Monocytecell, features = all.genes)

## PCA dimensionality reduction
Monocytecell <- RunPCA(Monocytecell, features = VariableFeatures(object = Monocytecell))
VizDimLoadings(Monocytecell, dims = 1:2, reduction = "pca")
Monocytecell <- JackStraw(Monocytecell, num.replicate = 100, dims = 50)
Monocytecell <- ScoreJackStraw(Monocytecell, dims = 1:50)
p1 <- JackStrawPlot(Monocytecell, dims = 1:50)
p2 <- ElbowPlot(Monocytecell, ndims = 50)
pdf(file = "path/image/monocyte_cell/ElbowPlot.pdf", width = 13, height = 3)
p1+p2
dev.off()

### UMAP dimensionality reduction
Monocytecell <- Monocytecell %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.9) %>% identity()
pdf(file = "path/image/monocyte_cell/Monocytecell_UMAP_cluster_r0.9.pdf", width = 5.5, height = 4)
DimPlot(Monocytecell, reduction = "umap", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/monocyte_cell/Monocytecell_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(Monocytecell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

### TSNE dimensionality reduction
Monocytecell <- RunTSNE(Monocytecell, reduction = "harmony",dims = 1:50)
pdf(file = "path/image/monocyte_cell/Monocytecell_TSNE_cluster_r0.9.pdf", width = 5.5, height = 4)
DimPlot(Monocytecell, reduction = "tsne", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/monocyte_cell/Monocytecell_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(Monocytecell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

features <- c("CD14", "FCGR3A")
pdf(file = "path/image/monocyte_cell/Plot_markergene_0.9.pdf", width = 5, height = 5)
DotPlot(Monocytecell, feature = features, cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()
saveRDS(Monocytecell, file = "path/result/Monocyte_cell/Monocytecell_harmony_UMAP_TSNE_r0.9.rds")

######### Annotated cell cluster #########
current.cluster.ids <- c(0:20)
new.cluster.ids <- c("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14lowCD16high", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", 
                     "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14highCD16low", "Monocyte_CD14lowCD16high", 
					 "Monocyte_CD14highCD16low", "Monocyte_CD14lowCD16high", "Monocyte_CD14highCD16low", "Monocyte_CD14lowCD16high","Monocyte_CD14highCD16high")
Monocytecell@meta.data$celltype_Monocyte = plyr::mapvalues(x = Monocytecell@meta.data[,"RNA_snn_res.0.9"], from = current.cluster.ids, to = new.cluster.ids)
Monocytecell$celltype_Monocyte <- factor(Monocytecell$celltype_Monocyte,level = c("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14lowCD16high"))

## Cell type show
pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_UMAP_celltype.pdf", width = 6, height = 4)
DimPlot(Monocytecell, reduction = "umap", group.by = "celltype_Monocyte", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_Monocytecell_tSNE_celltype.pdf", width = 6, height = 4)
DimPlot(Monocytecell, reduction = "tsne", group.by = "celltype_Monocyte", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_UMAP_group.pdf", width = 5, height = 4)
DimPlot(Monocytecell, reduction = "umap", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_Monocytecell_tSNE_group.pdf", width = 5, height = 4)
DimPlot(Monocytecell, reduction = "tsne", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_UMAP_level.pdf", width = 5, height = 4)
DimPlot(Monocytecell, reduction = "umap", group.by = "level", cols = c('#2967a0', '#ecc342', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_Monocytecell_tSNE_level.pdf", width = 5, height = 4)
DimPlot(Monocytecell, reduction = "tsne", group.by = "level", cols = c('#2967a0', '#ecc342', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

Monocytecell$orig.ident <- factor(Monocytecell$orig.ident,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(Monocytecell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/monocyte_cell/Monocytecell_harmony_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(Monocytecell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Show marker genes
features <- c()
pdf(file = "path/image/monocyte_cell/VlnPlot_markergene_celltype_Monocyte.pdf", width = 7, height = 2)
VlnPlot(Monocytecell, features = features, group.by = "celltype_Monocyte", stack = T, cols = mycolour)
dev.off()

## Find DEGs among cell types
Idents(Monocytecell) = "celltype_Monocyte"
diff.wilcox = FindAllMarkers(Monocytecell, min.pct = 0.25, logfc.threshold = 0.25)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(diff.wilcox, "path/result/Monocyte_cell/diff_genes_wilcox_celltype_Monocyte.csv")

pdf(file = "path/image/monocyte_cell/diffgene_volcano_polar_celltype.pdf", width = 8, height = 8)
jjVolcano(diffData = diff.wilcox, log2FC.cutoff = 0.5, topGeneN = 3, tile.col = mycolour2, size = 2, polar = T)
dev.off()

## Cell percentage analysis

### Sample proportion in each cell type
df1 <- table(Monocytecell$celltype_Monocyte,Monocytecell$orig.ident) %>% melt()
colnames(df1) <- c("Cluster","Sample","Number")
df1$Cluster <- factor(df1$Cluster,level = c ("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14lowCD16high"))
df1$Sample <- factor(df1$Sample,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/monocyte_cell/Monocytecell_celltype_sample_percent.pdf", width = 5, height = 3)
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
df2 <- table(Monocytecell$celltype_Monocyte,Monocytecell$group) %>% melt()
colnames(df2) <- c("Cluster","Group","Number")
df2$Cluster <- factor(df2$Cluster,level = c ("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14lowCD16high"))
df2$Group <- factor(df2$Group,level = c ("H", "P"))
pdf(file = "path/image/monocyte_cell/Monocytecell_celltype_group_percent.pdf", width = 4.5, height = 3)
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
df3 <- table(Monocytecell$celltype_Monocyte,Monocytecell$level) %>% melt()
colnames(df3) <- c("Cluster","Level","Number")
df3$Cluster <- factor(df3$Cluster,level = c ("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14lowCD16high"))
df3$Level<- factor(df3$Level,level = c ("Control", "Mild", "Severe"))
pdf(file = "path/image/monocyte_cell/Monocytecell_celltype_level_percent.pdf", width = 4.5, height = 3)
ggplot(data = df3, aes(x =Number, y = Cluster, fill =  Level)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c('#2967a0', '#ecc342', '#af2337')) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
dev.off()

### Different ercentage of each Monocyte cell subtype in group

table(Monocytecell$orig.ident)
prop.table(table(Idents(Monocytecell)))
table(Idents(Monocytecell), Monocytecell$orig.ident)
Cellratio <- prop.table(table(Idents(Monocytecell), Monocytecell$orig.ident), margin = 2)
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
sce_groups = c("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14lowCD16high")
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
    scale_fill_manual(values=c('#2967a0', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("H", "P") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/monocyte_cell/Monocytecell_group_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['Monocyte_CD14highCD16low']],
          pplist[['Monocyte_CD14highCD16high']],
          pplist[['Monocyte_CD14lowCD16high']]
		  )
dev.off()

### Different ercentage of each Monocyte cell subtype in level
table(Monocytecell$orig.ident)
prop.table(table(Idents(Monocytecell)))
table(Idents(Monocytecell), Monocytecell$orig.ident)
Cellratio <- prop.table(table(Idents(Monocytecell), Monocytecell$orig.ident), margin = 2)
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
sce_groups = c("Monocyte_CD14highCD16low", "Monocyte_CD14highCD16high", "Monocyte_CD14lowCD16high")
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
                                                      median = median(percent))#
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent,fill=group)) + 
    geom_boxplot() + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    scale_fill_manual(values=c('#2967a0', '#ecc342', '#af2337'))
  
  ### Inter-group t-test
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list(c("Control", "Mild"), c("Mild", "Severe"), c("Control", "Severe"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

pdf(file = "path/image/monocyte_cell/Monocytecell_level_celltype_percent_diff.pdf", width = 10, height = 9)
plot_grid(pplist[['Monocyte_CD14highCD16low']],
          pplist[['Monocyte_CD14highCD16high']],
          pplist[['Monocyte_CD14lowCD16high']]
		  )
dev.off()
saveRDS(Monocytecell, file = "path/result/Monocyte_cell/Monocytecell_harmony_UMAP_Monocyte_tSNE_celltype.rds")

####### compare the difference of one geneset between two groups by QuSAGE ########

Antigen_genelist <- readLines("path/result/Monocyte_cell/Antigen_presentation_genelist_human.txt")
Phagocytosis_genelist <- readLines("path/result/Monocyte_cell/Phagocytosis_genelist_human.txt")
Inflammatory_genelist <- readLines("path/result/Monocyte_cell/inflammatory_response_genelist_human.txt")

####### Control vs Severe
Monocytecell_2 <- subset(Monocytecell, level %in% c("Control", "Severe"))

# Monocyte_CD14highCD16low

## 1.calculation
eset=as.data.frame(Monocytecell_2[["RNA"]]$data[,rownames(Monocytecell_2@meta.data)[which(Monocytecell_2$celltype_Monocyte %in% c("Monocyte_CD14highCD16low"))]])
labels=as.character(paste(Monocytecell_2$celltype_Monocyte[colnames(eset)],Monocytecell_2$level[colnames(eset)],sep="."))
sf_Monocyte_Antigen=qusage(eset,labels,"Monocyte_CD14highCD16low.Severe-Monocyte_CD14highCD16low.Control",Antigen_genelist)
sf_Monocyte_Phagocytosis=qusage(eset,labels,"Monocyte_CD14highCD16low.Severe-Monocyte_CD14highCD16low.Control",Phagocytosis_genelist)
sf_Monocyte_Inflammatory=qusage(eset,labels,"Monocyte_CD14highCD16low.Severe-Monocyte_CD14highCD16low.Control",Inflammatory_genelist)

# p value adjustment
p.vals_Monocyte_Antigen=pdf.pVal(sf_Monocyte_Antigen)
q.vals_Monocyte_Antigen = p.adjust(p.vals_Monocyte_Antigen, method="fdr")
p.vals_Monocyte_Phagocytosis=pdf.pVal(sf_Monocyte_Phagocytosis)
q.vals_Monocyte_Phagocytosis = p.adjust(p.vals_Monocyte_Phagocytosis, method="fdr")
p.vals_Monocyte_Inflammatory=pdf.pVal(sf_Monocyte_Inflammatory)
q.vals_Monocyte_Inflammatory = p.adjust(p.vals_Monocyte_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/monocyte_cell/ControlvsSevere_genesetScore_Monocyte_CD14highCD16low_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_Monocyte_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Monocyte_CD14highCD16low: Severe VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.4,0),ylim=c(0,320))
plotDensityCurves(sf_Monocyte_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_Monocyte_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("Monocyte_Phagocytosis","Monocyte_Inflammatory","Monocyte_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(80), labels =transPvalue(q.vals_Monocyte_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(100), labels = transPvalue(q.vals_Monocyte_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(50), labels =transPvalue(q.vals_Monocyte_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()

# Monocyte_CD14highCD16high

## 1.calculation
eset=as.data.frame(Monocytecell_2[["RNA"]]$data[,rownames(Monocytecell_2@meta.data)[which(Monocytecell_2$celltype_Monocyte %in% c("Monocyte_CD14highCD16high"))]])
labels=as.character(paste(Monocytecell_2$celltype_Monocyte[colnames(eset)],Monocytecell_2$level[colnames(eset)],sep="."))
sf_Monocyte_Antigen=qusage(eset,labels,"Monocyte_CD14highCD16high.Severe-Monocyte_CD14highCD16high.Control",Antigen_genelist)
sf_Monocyte_Phagocytosis=qusage(eset,labels,"Monocyte_CD14highCD16high.Severe-Monocyte_CD14highCD16high.Control",Phagocytosis_genelist)
sf_Monocyte_Inflammatory=qusage(eset,labels,"Monocyte_CD14highCD16high.Severe-Monocyte_CD14highCD16high.Control",Inflammatory_genelist)

# p value adjustment
p.vals_Monocyte_Antigen=pdf.pVal(sf_Monocyte_Antigen)
q.vals_Monocyte_Antigen = p.adjust(p.vals_Monocyte_Antigen, method="fdr")
p.vals_Monocyte_Phagocytosis=pdf.pVal(sf_Monocyte_Phagocytosis)
q.vals_Monocyte_Phagocytosis = p.adjust(p.vals_Monocyte_Phagocytosis, method="fdr")
p.vals_Monocyte_Inflammatory=pdf.pVal(sf_Monocyte_Inflammatory)
q.vals_Monocyte_Inflammatory = p.adjust(p.vals_Monocyte_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/monocyte_cell/ControlvsSevere_genesetScore_Monocyte_CD14highCD16high_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_Monocyte_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Monocyte_CD14highCD16high: Severe VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.4,0.2),ylim=c(0,70))
plotDensityCurves(sf_Monocyte_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_Monocyte_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("Monocyte_Phagocytosis","Monocyte_Inflammatory","Monocyte_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(20), labels =transPvalue(q.vals_Monocyte_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(30), labels = transPvalue(q.vals_Monocyte_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(40), labels =transPvalue(q.vals_Monocyte_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()

# Monocyte_CD14lowCD16high

## 1.calculation
eset=as.data.frame(Monocytecell_2[["RNA"]]$data[,rownames(Monocytecell_2@meta.data)[which(Monocytecell_2$celltype_Monocyte %in% c("Monocyte_CD14lowCD16high"))]])
labels=as.character(paste(Monocytecell_2$celltype_Monocyte[colnames(eset)],Monocytecell_2$level[colnames(eset)],sep="."))
sf_Monocyte_Antigen=qusage(eset,labels,"Monocyte_CD14lowCD16high.Severe-Monocyte_CD14lowCD16high.Control",Antigen_genelist)
sf_Monocyte_Phagocytosis=qusage(eset,labels,"Monocyte_CD14lowCD16high.Severe-Monocyte_CD14lowCD16high.Control",Phagocytosis_genelist)
sf_Monocyte_Inflammatory=qusage(eset,labels,"Monocyte_CD14lowCD16high.Severe-Monocyte_CD14lowCD16high.Control",Inflammatory_genelist)

# p value adjustment
p.vals_Monocyte_Antigen=pdf.pVal(sf_Monocyte_Antigen)
q.vals_Monocyte_Antigen = p.adjust(p.vals_Monocyte_Antigen, method="fdr")
p.vals_Monocyte_Phagocytosis=pdf.pVal(sf_Monocyte_Phagocytosis)
q.vals_Monocyte_Phagocytosis = p.adjust(p.vals_Monocyte_Phagocytosis, method="fdr")
p.vals_Monocyte_Inflammatory=pdf.pVal(sf_Monocyte_Inflammatory)
q.vals_Monocyte_Inflammatory = p.adjust(p.vals_Monocyte_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/monocyte_cell/ControlvsSevere_genesetScore_Monocyte_CD14lowCD16high_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_Monocyte_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Monocyte_CD14lowCD16high: Severe VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.4,0),ylim=c(0,130))
plotDensityCurves(sf_Monocyte_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_Monocyte_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("Monocyte_Phagocytosis","Monocyte_Inflammatory","Monocyte_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(80), labels =transPvalue(q.vals_Monocyte_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(100), labels = transPvalue(q.vals_Monocyte_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(50), labels =transPvalue(q.vals_Monocyte_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()

####### Control vs Mild

Monocytecell_3 <- subset(Monocytecell, level %in% c("Control", "Mild"))

# Monocyte_CD14highCD16low

## 1.calculation
eset=as.data.frame(Monocytecell_3[["RNA"]]$data[,rownames(Monocytecell_3@meta.data)[which(Monocytecell_3$celltype_Monocyte %in% c("Monocyte_CD14highCD16low"))]])
labels=as.character(paste(Monocytecell_3$celltype_Monocyte[colnames(eset)],Monocytecell_3$level[colnames(eset)],sep="."))
sf_Monocyte_Antigen=qusage(eset,labels,"Monocyte_CD14highCD16low.Mild-Monocyte_CD14highCD16low.Control",Antigen_genelist)
sf_Monocyte_Phagocytosis=qusage(eset,labels,"Monocyte_CD14highCD16low.Mild-Monocyte_CD14highCD16low.Control",Phagocytosis_genelist)
sf_Monocyte_Inflammatory=qusage(eset,labels,"Monocyte_CD14highCD16low.Mild-Monocyte_CD14highCD16low.Control",Inflammatory_genelist)

# p value adjustment
p.vals_Monocyte_Antigen=pdf.pVal(sf_Monocyte_Antigen)
q.vals_Monocyte_Antigen = p.adjust(p.vals_Monocyte_Antigen, method="fdr")
p.vals_Monocyte_Phagocytosis=pdf.pVal(sf_Monocyte_Phagocytosis)
q.vals_Monocyte_Phagocytosis = p.adjust(p.vals_Monocyte_Phagocytosis, method="fdr")
p.vals_Monocyte_Inflammatory=pdf.pVal(sf_Monocyte_Inflammatory)
q.vals_Monocyte_Inflammatory = p.adjust(p.vals_Monocyte_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/monocyte_cell/ControlvsMild_genesetScore_Monocyte_CD14highCD16low_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_Monocyte_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Monocyte_CD14highCD16low: Mild VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.3,0.1),ylim=c(0,300))
plotDensityCurves(sf_Monocyte_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_Monocyte_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("Monocyte_Phagocytosis","Monocyte_Inflammatory","Monocyte_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(80), labels =transPvalue(q.vals_Monocyte_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(100), labels = transPvalue(q.vals_Monocyte_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(50), labels =transPvalue(q.vals_Monocyte_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()

# Monocyte_CD14highCD16high

## 1.calculation
eset=as.data.frame(Monocytecell_3[["RNA"]]$data[,rownames(Monocytecell_3@meta.data)[which(Monocytecell_3$celltype_Monocyte %in% c("Monocyte_CD14highCD16high"))]])
labels=as.character(paste(Monocytecell_3$celltype_Monocyte[colnames(eset)],Monocytecell_3$level[colnames(eset)],sep="."))
sf_Monocyte_Antigen=qusage(eset,labels,"Monocyte_CD14highCD16high.Mild-Monocyte_CD14highCD16high.Control",Antigen_genelist)
sf_Monocyte_Phagocytosis=qusage(eset,labels,"Monocyte_CD14highCD16high.Mild-Monocyte_CD14highCD16high.Control",Phagocytosis_genelist)
sf_Monocyte_Inflammatory=qusage(eset,labels,"Monocyte_CD14highCD16high.Mild-Monocyte_CD14highCD16high.Control",Inflammatory_genelist)

# p value adjustment
p.vals_Monocyte_Antigen=pdf.pVal(sf_Monocyte_Antigen)
q.vals_Monocyte_Antigen = p.adjust(p.vals_Monocyte_Antigen, method="fdr")
p.vals_Monocyte_Phagocytosis=pdf.pVal(sf_Monocyte_Phagocytosis)
q.vals_Monocyte_Phagocytosis = p.adjust(p.vals_Monocyte_Phagocytosis, method="fdr")
p.vals_Monocyte_Inflammatory=pdf.pVal(sf_Monocyte_Inflammatory)
q.vals_Monocyte_Inflammatory = p.adjust(p.vals_Monocyte_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/monocyte_cell/ControlvsMild_genesetScore_Monocyte_CD14highCD16high_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_Monocyte_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Monocyte_CD14highCD16high: Mild VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.4,0.2),ylim=c(0,50))
plotDensityCurves(sf_Monocyte_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_Monocyte_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("Monocyte_Phagocytosis","Monocyte_Inflammatory","Monocyte_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(20), labels =transPvalue(q.vals_Monocyte_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(30), labels = transPvalue(q.vals_Monocyte_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(40), labels =transPvalue(q.vals_Monocyte_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()

# Monocyte_CD14lowCD16high

## 1.calculation
eset=as.data.frame(Monocytecell_3[["RNA"]]$data[,rownames(Monocytecell_3@meta.data)[which(Monocytecell_3$celltype_Monocyte %in% c("Monocyte_CD14lowCD16high"))]])
labels=as.character(paste(Monocytecell_3$celltype_Monocyte[colnames(eset)],Monocytecell_3$level[colnames(eset)],sep="."))
sf_Monocyte_Antigen=qusage(eset,labels,"Monocyte_CD14lowCD16high.Mild-Monocyte_CD14lowCD16high.Control",Antigen_genelist)
sf_Monocyte_Phagocytosis=qusage(eset,labels,"Monocyte_CD14lowCD16high.Mild-Monocyte_CD14lowCD16high.Control",Phagocytosis_genelist)
sf_Monocyte_Inflammatory=qusage(eset,labels,"Monocyte_CD14lowCD16high.Mild-Monocyte_CD14lowCD16high.Control",Inflammatory_genelist)

# p value adjustment
p.vals_Monocyte_Antigen=pdf.pVal(sf_Monocyte_Antigen)
q.vals_Monocyte_Antigen = p.adjust(p.vals_Monocyte_Antigen, method="fdr")
p.vals_Monocyte_Phagocytosis=pdf.pVal(sf_Monocyte_Phagocytosis)
q.vals_Monocyte_Phagocytosis = p.adjust(p.vals_Monocyte_Phagocytosis, method="fdr")
p.vals_Monocyte_Inflammatory=pdf.pVal(sf_Monocyte_Inflammatory)
q.vals_Monocyte_Inflammatory = p.adjust(p.vals_Monocyte_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/monocyte_cell/ControlvsMild_genesetScore_Monocyte_CD14lowCD16high_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_Monocyte_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Monocyte_CD14lowCD16high: Mild VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.3,0),ylim=c(0,150))
plotDensityCurves(sf_Monocyte_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_Monocyte_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("Monocyte_Phagocytosis","Monocyte_Inflammatory","Monocyte_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(80), labels =transPvalue(q.vals_Monocyte_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(100), labels = transPvalue(q.vals_Monocyte_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(50), labels =transPvalue(q.vals_Monocyte_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()