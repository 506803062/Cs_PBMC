############ NK cell subcluster ############

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scRNAtoolVis)
library(ggrepel)
library(clusterProfiler)
library(SCENIC)
library(pheatmap)

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
		 
# Extract NK cell		 
pbmc <- readRDS("path/result/pbmc_harmony_UMAP_TSNE_celltype_RNA.rds")
NKcell <- subset(pbmc, celltype %in% c("NK"))
saveRDS(NKcell, file = "path/result/NK_cell/pbmc_NK.rds")

################## Find DEGs in NK cell in level ##################
diff_NK_cell_CvsM <- FindMarkers(NKcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Mild")
diff_NK_cell_CvsS <- FindMarkers(NKcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "level",ident.1 = "Control", ident.2 = "Severe")
write.csv(diff_NK_cell_CvsM, "path/result/NK_cell/NKcell_diff_genes_wilcox_Control_vs_Mild.csv", row.names = T)
write.csv(diff_NK_cell_CvsS, "path/result/NK_cell/NKcell_diff_genes_wilcox_Control_vs_Severe.csv", row.names = T)

## Draw Control_vs_Mild volcano plot
diff_NK_cell_CvsM[which(diff_NK_cell_CvsM$avg_log2FC  >= 0.25 & diff_NK_cell_CvsM$p_val < 0.05),'sig'] <- 'up'
diff_NK_cell_CvsM[which(diff_NK_cell_CvsM$avg_log2FC  <= -0.25 & diff_NK_cell_CvsM$p_val < 0.05),'sig'] <- 'down'
diff_NK_cell_CvsM$symbol <- rownames(diff_NK_cell_CvsM)

# Find top 10 up-regulated genes
up_data <- filter(diff_NK_cell_CvsM, sig == 'up') %>%  
  distinct(symbol, .keep_all = TRUE) %>%              
  top_n(3, avg_log2FC)                          

# Find top 10 down-regulated genes
down_data <- filter(diff_NK_cell_CvsM, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                   
  top_n(-3, avg_log2FC)                               

pdf(file = "path/image/NK_cell/NKcell_Control_vs_Mild_volcano.pdf", width = 4, height = 4)
ggplot(diff_NK_cell_CvsM, aes(avg_log2FC , -log10(p_val), col = sig)) +
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
diff_NK_cell_CvsS[which(diff_NK_cell_CvsS$avg_log2FC  >= 0.25 & diff_NK_cell_CvsS$p_val < 0.05),'sig'] <- 'up'
diff_NK_cell_CvsS[which(diff_NK_cell_CvsS$avg_log2FC  <= -0.25 & diff_NK_cell_CvsS$p_val < 0.05),'sig'] <- 'down'
diff_NK_cell_CvsS$symbol <- rownames(diff_NK_cell_CvsS)

# Find top 10 up-regulated genes
up_data <- filter(diff_NK_cell_CvsS, sig == 'up') %>% 
  distinct(symbol, .keep_all = TRUE) %>%             
  top_n(3, avg_log2FC)                          

# Find top 10 down-regulated genes
down_data <- filter(diff_NK_cell_CvsS, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                  
  top_n(-3, avg_log2FC)                              

pdf(file = "path/image/NK_cell/NKcell_Control_vs_Severe_volcano.pdf", width = 4, height = 4)
ggplot(diff_NK_cell_CvsS, aes(avg_log2FC, -log10(p_val), col = sig)) +
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
Idents(NKcell) = "level"
diff_NKcell_level  <- FindAllMarkers(NKcell, only.pos = TRUE,
                       min.pct = 0.25, 
                       logfc.threshold = 0.75)
diff_NKcell_level_sig  <- diff_NKcell_level[diff_NKcell_level$p_val_adj < 0.05, ]
write.csv(diff_NKcell_level, "path/result/NK_cell/NKcell_diff_genes_wilcox_Level.csv", quote = FALSE)

group <- data.frame(gene=diff_NKcell_level_sig$gene,
                    group=diff_NKcell_level_sig$cluster)

Gene_ID <- bitr(diff_NKcell_level_sig$gene, fromType="SYMBOL", 
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

pdf(file = "path/image/NK_cell/NKcell_level_GO.pdf", width = 3, height = 4.5)
dotplot(data_GO_sim, color = "p.adjust", size = "Count", showCategory=5, font.size = 5)
dev.off()

########### NK cell subcluster ###########
NKcell$orig.ident <- factor(NKcell$orig.ident,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
NKcell <- FindVariableFeatures(NKcell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NKcell)
NKcell <- ScaleData(NKcell, features = all.genes)

## PCA dimensionality reduction
NKcell <- RunPCA(NKcell, features = VariableFeatures(object = NKcell))
VizDimLoadings(NKcell, dims = 1:2, reduction = "pca")
NKcell <- JackStraw(NKcell, num.replicate = 100, dims = 50)
NKcell <- ScoreJackStraw(NKcell, dims = 1:50)
p1 <- JackStrawPlot(NKcell, dims = 1:50)
p2 <- ElbowPlot(NKcell, ndims = 50)
pdf(file = "path/image/NK_cell/ElbowPlot.pdf", width = 13, height = 3)
p1+p2
dev.off()

### UMAP dimensionality reduction
NKcell <- NKcell %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.5) %>% identity()
pdf(file = "path/image/NK_cell/NKcell_UMAP_cluster_r0.5.pdf", width = 5.5, height = 4)
DimPlot(NKcell, reduction = "umap", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/NK_cell/NKcell_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(NKcell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

### TSNE dimensionality reduction
NKcell <- RunTSNE(NKcell, reduction = "harmony",dims = 1:50)
pdf(file = "path/image/NK_cell/NKcell_TSNE_cluster_r0.5.pdf", width = 5.5, height = 4)
DimPlot(NKcell, reduction = "tsne", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/NK_cell/NKcell_TSNE_sample.pdf", width = 5.5, height = 4)
DimPlot(NKcell, reduction = "tsne", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
saveRDS(NKcell, file = "path/result/NK_cell/NKcell_harmony_UMAP_TSNE_r0.5.rds")

######### Annotated cell cluster #########
current.cluster.ids <- c(0:15)
new.cluster.ids <- c("NK_CD16high", "NK_CD16high", "NK_CD16low", "NK_CD16high", "NK_CD16high", "NK_CD16high", "NK_CD16high", "NK_CD16high", "NK_CD16low", 
                     "NK_CD16low", "NK_CD16high", "NK_CD16low", "NK_CD16high", "NK_CD16high", "NK_CD16low", "NK_CD16high")
NKcell@meta.data$celltype_NK = plyr::mapvalues(x = NKcell@meta.data[,"RNA_snn_res.0.5"], from = current.cluster.ids, to = new.cluster.ids)
NKcell$celltype_NK <- factor(NKcell$celltype_NK,level = c("NK_CD16high", "NK_CD16low"))

## Cell type show
pdf(file = "path/image/NK_cell/NKcell_harmony_UMAP_celltype.pdf", width = 6, height = 4)
DimPlot(NKcell, reduction = "umap", group.by = "celltype_NK", cols = c('#D6E7A3', '#57C3F3')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/NK_cell/NKcell_harmony_TSNE_celltype.pdf", width = 6, height = 4)
DimPlot(NKcell, reduction = "tsne", group.by = "celltype_NK", cols = c('#D6E7A3', '#57C3F3')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/NK_cell/NKcell_harmony_UMAP_group.pdf", width = 5, height = 4)
DimPlot(NKcell, reduction = "umap", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/NK_cell/NKcell_harmony_TSNE_group.pdf", width = 5, height = 4)
DimPlot(NKcell, reduction = "tsne", group.by = "group", cols = c('#2967a0', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/image/NK_cell/NKcell_harmony_UMAP_level.pdf", width = 5, height = 4)
DimPlot(NKcell, reduction = "umap", group.by = "level", cols = c('#96b437','#da93ab','#a9c9ed')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/NK_cell/NKcell_harmony_TSNE_level.pdf", width = 5, height = 4)
DimPlot(NKcell, reduction = "tsne", group.by = "level", cols = c('#96b437','#da93ab','#a9c9ed')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Show marker genes
features <- c("FCGR3A", "SPON2", "FGFBP2", "NKG7", "CST7", "GZMB", "CCL4", "CX3CR1", "CD247", "S100A4", "IGFBP4", "CAPG", "GZMK", "GPR183", "IL7R", "COTL1", "SELL", "TCF7")
pdf(file = "path/image/NK_cell/Plot_markergene_celltype_NK.pdf", width = 8, height = 6)
jjDotPlot(NKcell, gene = features, id = "celltype_NK", ytree = F)
dev.off()

## Cell percentage analysis

### Sample proportion in each cell type
df1 <- table(NKcell$celltype_NK,NKcell$orig.ident) %>% melt()
colnames(df1) <- c("Cluster","Sample","Number")
df1$Cluster <- factor(df1$Cluster,level = c ("NK_CD16high", "NK_CD16low"))
df1$Sample <- factor(df1$Sample,level = c ("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12"))
pdf(file = "path/image/NK_cell/NKcell_celltype_sample_percent.pdf", width = 5, height = 3)
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
df2 <- table(NKcell$celltype_NK,NKcell$group) %>% melt()
colnames(df2) <- c("Cluster","Group","Number")
df2$Cluster <- factor(df2$Cluster,level = c ("NK_CD16high", "NK_CD16low"))
df2$Group <- factor(df2$Group,level = c ("H", "P"))
pdf(file = "path/image/NK_cell/NKcell_celltype_group_percent.pdf", width = 4.5, height = 3)
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
df3 <- table(NKcell$celltype_NK,NKcell$level) %>% melt()
colnames(df3) <- c("Cluster","Level","Number")
df3$Cluster <- factor(df3$Cluster,level = c ("NK_CD16high", "NK_CD16low"))
df3$Level<- factor(df3$Level,level = c ("Control", "Mild", "Severe"))
pdf(file = "path/image/NK_cell/NKcell_celltype_level_percent.pdf", width = 4.5, height = 3)
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

##################### Analysis TF in NK cell TF bySCENIC #####################

## Cell meta information
Bcell <- readRDS(file = "path/result/NK_cell/Bcell_harmony_UMAP_TSNE_celltype.rds")
cellInfo <- data.frame(Bcell@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_B")] <-"celltype_B"
cellInfo <- cellInfo[,c("sample","cluster","celltype_B")]
saveRDS(cellInfo, file="path/result/NK_cell/SCENIC/int/cellInfo.Rds")

## Expression matrix
exprMat <- as.matrix(Bcell@assays$RNA$counts)
mydbDIR <- "path/result/NK_cell/SCENIC/cisTarget"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp","10kb")
scenicOptions <- initializeScenic(org="hgnc",nCores=50,dbDir=mydbDIR,dbs = mydbs,datasetTitle ="NK_cell")
saveRDS(scenicOptions, "path/result/NK_cell/SCENIC/int/scenicOptions.rds")

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

cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select='celltype_B')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

# Draw a heat map using regulon raw AUC values
celltypecolor <- c('#af2337', '#ecc342', '#2967a0', '#2f3c28')
names(celltypecolor) <- c("Follicular NK cell", "Megakaryocyte-like cell", "Naive NK cell", "Memory NK cell")
ann_colors <- list(celltype_B = celltypecolor)
pdf(file = "path/image/NK_cell/Bcell_regulon_AUC.pdf", width = 8, height = 5)
pheatmap(AUCmatrix, show_colnames=F, annotation_col=celltype, annotation_colors = ann_colors)
dev.off()

# Draw heat maps using regulon binary AUC values
pdf(file = "path/image/NK_cell/Bcell_regulon_Bin.pdf", width = 8, height = 5)
pheatmap(BINmatrix, show_colnames=F, annotation_col=celltype, annotation_colors = ann_colors,color = colorRampPalette(colors = c("white","black"))(100))
dev.off()

####### compare the difference of one geneset between two groups by QuSAGE ########

Antigen_genelist <- readLines("path/result/Antigen_presentation_genelist_human.txt")
Cytotoxicity_genelist <- readLines("path/result/Cytotoxicity_genelist_human.txt")
Inflammatory_genelist <- readLines("path/result/inflammatory_response_genelist_human.txt")

####### Control vs Severe

## 1.calculation
eset=as.data.frame(NKcell[["RNA"]]$data[,rownames(NKcell@meta.data)[which(NKcell$celltype_NK %in% c("NK"))]])
labels=as.character(paste(NKcell$celltype_NK[colnames(eset)],NKcell$level[colnames(eset)],sep="."))
sf_NK_Antigen=qusage(eset,labels,"NK.Severe-NK.Control",Antigen_genelist)
sf_NK_Phagocytosis=qusage(eset,labels,"NK.Severe-NK.Control",Phagocytosis_genelist)
sf_NK_Inflammatory=qusage(eset,labels,"NK.Severe-NK.Control",Inflammatory_genelist)

# p value adjustment
p.vals_NK_Antigen=pdf.pVal(sf_NK_Antigen)
q.vals_NK_Antigen = p.adjust(p.vals_NK_Antigen, method="fdr")
p.vals_NK_Phagocytosis=pdf.pVal(sf_NK_Phagocytosis)
q.vals_NK_Phagocytosis = p.adjust(p.vals_NK_Phagocytosis, method="fdr")
p.vals_NK_Inflammatory=pdf.pVal(sf_NK_Inflammatory)
q.vals_NK_Inflammatory = p.adjust(p.vals_NK_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/NK_cell/ControlvsSevere_genesetScore_NK_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_NK_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("NK: Severe VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.4,0),ylim=c(0,320))
plotDensityCurves(sf_NK_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_NK_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("NK_Phagocytosis","NK_Inflammatory","NK_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(80), labels =transPvalue(q.vals_NK_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(100), labels = transPvalue(q.vals_NK_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(50), labels =transPvalue(q.vals_NK_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()

####### Control vs Mild

## 1.calculation
eset=as.data.frame(NKcell[["RNA"]]$data[,rownames(NKcell@meta.data)[which(NKcell$celltype_NK %in% c("NK"))]])
labels=as.character(paste(NKcell$celltype_NK[colnames(eset)],NKcell$level[colnames(eset)],sep="."))
sf_NK_Antigen=qusage(eset,labels,"NK.Mild-NK.Control",Antigen_genelist)
sf_NK_Phagocytosis=qusage(eset,labels,"NK.Mild-NK.Control",Phagocytosis_genelist)
sf_NK_Inflammatory=qusage(eset,labels,"NK.Mild-NK.Control",Inflammatory_genelist)

# p value adjustment
p.vals_NK_Antigen=pdf.pVal(sf_NK_Antigen)
q.vals_NK_Antigen = p.adjust(p.vals_NK_Antigen, method="fdr")
p.vals_NK_Phagocytosis=pdf.pVal(sf_NK_Phagocytosis)
q.vals_NK_Phagocytosis = p.adjust(p.vals_NK_Phagocytosis, method="fdr")
p.vals_NK_Inflammatory=pdf.pVal(sf_NK_Inflammatory)
q.vals_NK_Inflammatory = p.adjust(p.vals_NK_Inflammatory, method="fdr")

## 2.plot
pdf(file = "path/image/NK_cell/ControlvsMild_genesetScore_NK_QuSAGE.pdf", width = 5.5, height = 4.5)
plotDensityCurves(sf_NK_Phagocytosis,col=mycolour2[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("NK: Mild VS Control"),
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.3,0.1),ylim=c(0,300))
plotDensityCurves(sf_NK_Inflammatory,col=mycolour2[3],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_NK_Antigen,col=mycolour2[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("NK_Phagocytosis","NK_Inflammatory","NK_Antigen"),lty=1,lwd=2,col=mycolour2[c(2,3,1)],cex=1.2)
text(x = c(-0.1), y = c(80), labels =transPvalue(q.vals_NK_Phagocytosis), col = mycolour2[2],cex = c(0.8))
text(x = c(-0.1), y = c(100), labels = transPvalue(q.vals_NK_Inflammatory),col = mycolour2[3],cex = c(0.8))
text(x = c(-0.1), y = c(50), labels =transPvalue(q.vals_NK_Antigen), col = mycolour2[1],cex = c(0.8))
dev.off()