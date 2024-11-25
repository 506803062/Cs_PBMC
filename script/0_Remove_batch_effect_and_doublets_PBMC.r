######################### Remove batch effect and doublets #########################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)
library(DoubletFinder)

mycolour <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
              '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
              '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175')

#Import data and create Seurat object 
dir = c( "path/GSM6588511",
         "path/GSM6588512",
         "path/GSM6588513",
         "path/GSM6588514",
         "path/GSM6588515",
         "path/GSM6588516",
         "path/GSM6588517",
         "path/GSM6588518",
         "path/GSM6588519",
         "path/GSM6588520",
         "path/GSM6588521",
         "path/GSM6588522",
         "path/GSM6588523",
         "path/GSM6588524",
         "path/GSM6588525",
         "path/GSM6588526",
         "path/GSM6588527",
         "path/P1/outs/raw_feature_bc_matrix",
         "path/P2/outs/raw_feature_bc_matrix",
         "path/P3/outs/raw_feature_bc_matrix",
         "path/P4/outs/raw_feature_bc_matrix",
         "path/P5/outs/raw_feature_bc_matrix",
         "path/P6/outs/raw_feature_bc_matrix",
         "path/P7/outs/raw_feature_bc_matrix",
         "path/P8/outs/raw_feature_bc_matrix",
         "path/P9/outs/raw_feature_bc_matrix",
         "path/P10/outs/raw_feature_bc_matrix",
         "path/P11/outs/raw_feature_bc_matrix",
         "path/P12/outs/raw_feature_bc_matrix")
names(dir) = c("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
               "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
               "P7", "P8","P9", "P10", "P11", "P12")
pbmclist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  pbmclist[[i]] <- CreateSeuratObject(counts, min.cells=1)
}
pbmc <- merge(pbmclist[[1]], y=c(pbmclist[[2]], pbmclist[[3]],pbmclist[[4]], pbmclist[[5]], 
                                  pbmclist[[6]], pbmclist[[7]],pbmclist[[8]], pbmclist[[9]], 
                                  pbmclist[[10]], pbmclist[[11]], pbmclist[[12]], pbmclist[[13]],
                                  pbmclist[[14]], pbmclist[[15]], pbmclist[[16]], pbmclist[[17]],
                                  pbmclist[[18]], pbmclist[[19]], pbmclist[[20]], pbmclist[[21]],
                                  pbmclist[[22]], pbmclist[[23]], pbmclist[[24]], pbmclist[[25]], 
                                  pbmclist[[26]], pbmclist[[27]], pbmclist[[28]], pbmclist[[29]]))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
sample.ids <- c("H1", "H2", "H3","H4", "H5", "H6","H7", "H8", "H9","H10", "H11", "H12", 
                "H13","H14", "H15", "H16","H17", "P1", "P2", "P3", "P4", "P5", "P6", 
                "P7", "P8","P9", "P10", "P11", "P12")
group.ids <- c(rep("H", 17), rep("P", 12))
pbmc@meta.data$group = plyr::mapvalues(x = pbmc@meta.data[,"orig.ident"], from = sample.ids, to = group.ids)
level.ids <- c(rep("Control", 17), "Mild", "Mild", "Mild", "Severe", "Severe", "Severe", 
               "Severe", "Severe", "Mild", "Mild", "Severe", "Mild")
pbmc@meta.data$level = plyr::mapvalues(x = pbmc@meta.data[,"orig.ident"], from = sample.ids, to = level.ids)
pbmc <- AddMetaData(pbmc, metadata = c(rep("H", 17), rep("P", 12)), col.name = "group")
saveRDS(pbmc, file = "path/result/pbmc_original.rds")

##Quality control
pbmc <- readRDS("path/result/pbmc_original.rds")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
pdf(file = "path/image/QC.pdf", width = 6, height = 7)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0, cols = mycolour)
dev.off()

##Data standardization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##PCA dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
pdf(file = "path/image/All_PCA.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "pca", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
saveRDS(pbmc, file = "path/result/pbmc_pca.rds")

##Remove batch effect
pbmc <- readRDS("path/result/pbmc_pca.rds")
pbmc <- pbmc %>% RunHarmony("orig.ident", plot_convergence = TRUE)
pdf(file = "path/image/All_harmony.pdf", width = 5.5, height = 4)
DimPlot(pbmc, reduction = "harmony", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:50)
ElbowPlot(pbmc, ndims = 50)
saveRDS(pbmc, file = "path/result/pbmc_harmony.rds")

##Remove doublets
sweep.res.list <- paramSweep(pbmc, PCs = 1:50, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = 0.075
homotypic.prop <- modelHomotypic(pbmc$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(pbmc))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pbmc <- doubletFinder(pbmc, PCs = 1:50, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
pbmc <- subset(pbmc, DF.classifications_0.25_0.005_18230 %in% "Singlet")
saveRDS(pbmc, file = "path/result/pbmc_harmony_singlet.rds")