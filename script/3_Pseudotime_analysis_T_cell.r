############ T cell Pseudotime analysis ############

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(monocle3)
library(CytoTRACE2)

Tcell <- readRDS("path/result/T_cell/Tcell_harmony_UMAP_TSNE_celltype.rds")

############# CD4+ T cells #############
#Data import and processing
Tcell_CD4 <- subset(Tcell, celltype_T %in% c("CD4+ effector T cell", "CD4+ naive T cell", "CD4+ regulatory T cell", "CD4+ memory T cell"))
saveRDS(Tcell_CD4, file = "path/result/T_cell/Tcell_CD4.rds")

#Recluster
Tcell_CD4_2 <- Tcell_CD4
Tcell_CD4_2 <- FindVariableFeatures(Tcell_CD4_2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tcell_CD4_2)
Tcell_CD4_2 <- ScaleData(Tcell_CD4_2, features = all.genes)

##PCA dimensionality reduction
Tcell_CD4_2 <- RunPCA(Tcell_CD4_2, features = VariableFeatures(object = Tcell_CD4_2))
VizDimLoadings(Tcell_CD4_2, dims = 1:2, reduction = "pca")
Tcell_CD4_2 <- JackStraw(Tcell_CD4_2, num.replicate = 100, dims = 50)
Tcell_CD4_2 <- ScoreJackStraw(Tcell_CD4_2, dims = 1:50)

### UMAP dimensionality reduction
Tcell_CD4_2 <- Tcell_CD4_2 %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.9) %>% identity()
pdf(file = "path/image/T_cell/Tcell_CD4_UMAP_celltype.pdf", width = 6, height = 4)
DimPlot(Tcell_CD4_2, reduction = "umap", cols = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'), label = TRUE, group.by = "celltype_T") + theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Build data
data <- GetAssayData(Tcell_CD4_2, assay = 'RNA', slot = 'counts')
cell_metadata <- Tcell_CD4_2@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)

##Import the consolidated umap coordinates from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Tcell_CD4_2, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

## Monocle3 cluster
cds <- cluster_cells(cds)

## Find trajectory
cds <- learn_graph(cds)
pdf(file = "path/image/T_cell/Tcell_CD4_UMAP_celltype_track.pdf", width = 6, height = 4)
plot_cells(cds, color_cells_by = "celltype_T", label_groups_by_cluster = FALSE, label_leaves = FALSE, cell_size = 0.1,
               label_branch_points = FALSE, label_cell_groups = F) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_CD4_UMAP_group_track.pdf", width = 5, height = 4)
plot_cells(cds, color_cells_by = "group", label_groups_by_cluster = FALSE, label_leaves = FALSE, cell_size = 0.1,
               label_branch_points = FALSE, label_cell_groups = F) + scale_colour_manual(values = c('#af2337', '#2967a0'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_CD4_UMAP_level_track.pdf", width = 5, height = 4)
plot_cells(cds, color_cells_by = "level", label_groups_by_cluster = FALSE, label_leaves = FALSE, cell_size = 0.1,
               label_branch_points = FALSE, label_cell_groups = F) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Order cells by pseudotime

pdf(file = "path/image/T_cell/Tcell_CD4_UMAP_celltype_auxiliary_line.pdf", width = 4, height = 4)
plot_cells(cds, color_cells_by = "celltype_T", label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE) + geom_vline(xintercept = seq(-3,-2,0.25)) + geom_hline(yintercept = seq(4,5,0.25))
dev.off()
embed <- data.frame(Embeddings(Tcell_CD4_2, reduction = "umap"))
embed <- subset(embed, umap_1 > -2.75 & umap_1 < -2.5 & umap_2 > 4.2 & umap_2 < 4.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)

pdf(file = "path/image/T_cell/Tcell_CD4_UMAP_celltype_pseudotime.pdf", width = 5.5, height = 4)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE, cell_size = 0.1, label_branch_points = FALSE)+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Find DEGs about pseudotime
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=50)

Track_genes_sig <- c("CCR7", "GZMA", "FOXP3", "LTB")
# Gene expression trend map
pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_celltype_T.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 4)+ scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_group.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="group", 
                              min_expr=0.5, ncol = 4) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_level.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="level", 
                              min_expr=0.5, ncol = 4) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

saveRDS(Tcell_CD4_2, file = "path/result/T_cell/Tcell_CD4_reUMAP.rds")
saveRDS(cds, file = "path/result/T_cell/Tcell_CD4_pseudotime.rds")

gene_list_Naive <- c("CCR7", "RPS14", "EEF1A1", "RPS28", "HLA-E", "RPL29", "RPL13", "RPS18", "RPS19", "RPS7")
gene_list_Cytotoxicity <- c("GZMA", "ACTB", "B2M", "TMSB4X", "PFN1", "HLA-C", "CFL1", "HLA-A", "NKG7", "HLA-B")
gene_list_Immunoregulation <- c("FOXP3", "TRAC", "SMCHD1", "CTSS", "ELK3", "IL32", "TMEM248", "TAF1", "CD4", "CD68")
gene_list_Memory <- c("LTB", "HLA-B", "RPS28", "RPL23", "RPS19", "RPL7", "HLA-A", "RPS29", "MT-ATP6", "MT-ND5")

# Naive
pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_celltype_T_Naive.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Naive,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 6)+ scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_group_Naive.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Naive,], color_cells_by="group", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_level_Naive.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Naive,], color_cells_by="level", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

#Cytotoxicity
pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_celltype_T_Cytotoxicity.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Cytotoxicity,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 6)+ scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_group_Cytotoxicity.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Cytotoxicity,], color_cells_by="group", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_level_Cytotoxicity.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Cytotoxicity,], color_cells_by="level", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

# Immunoregulation
pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_celltype_T_Immunoregulation.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Immunoregulation,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 6)+ scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_group_Immunoregulation.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Immunoregulation,], color_cells_by="group", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_level_Immunoregulation.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Immunoregulation,], color_cells_by="level", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

# Memory
pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_celltype_T_Memory.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Memory,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 6)+ scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0', '#2f3c28'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_group_Memory.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Memory,], color_cells_by="group", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD4_pseudotime_level_Memory.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Memory,], color_cells_by="level", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

#### Assessing differentiation potential score by CytoTRACE2

####### Import seurat object ###########
Tcell_CD4_2 <- readRDS("path/result/T_cell/Tcell_CD4_reUMAP.rds")
cytotrace2_result_sce <- cytotrace2(Tcell_CD4_2, 
                                is_seurat = TRUE, 
                                slot_type = "counts", 
                                species = 'human',
                                seed = 1234)
pdf(file = "path/image/T_cell/Tcell_CD4_CytoTRACE2.pdf", width = 5, height = 4)
FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 0.1) + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                                      "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + theme_bw() + theme(panel.grid=element_blank())
dev.off()

############# CD8+ T cells #############
# Data import and processing
Tcell_CD8 <- subset(Tcell, celltype_T %in% c("CD8+ effector memory T cell", "CD8+ effector T cell", "CD8+ naive T cell"))
saveRDS(Tcell_CD8, file = "path/result/T_cell/Tcell_CD8.rds")

Tcell_CD8 <- readRDS("path/result/T_cell/Tcell_CD8.rds")
Tcell_CD8_2 <- Tcell_CD8

# Recluster
Tcell_CD8_2 <- FindVariableFeatures(Tcell_CD8_2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tcell_CD8_2)
Tcell_CD8_2 <- ScaleData(Tcell_CD8_2, features = all.genes)

## PCA dimensionality reduction
Tcell_CD8_2 <- RunPCA(Tcell_CD8_2, features = VariableFeatures(object = Tcell_CD8_2))
VizDimLoadings(Tcell_CD8_2, dims = 1:2, reduction = "pca")
Tcell_CD8_2 <- JackStraw(Tcell_CD8_2, num.replicate = 100, dims = 50)
Tcell_CD8_2 <- ScoreJackStraw(Tcell_CD8_2, dims = 1:50)

### UMAP Tcell_CD4_2 <- Tcell_CD4
Tcell_CD8_2 <- Tcell_CD8_2 %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.9) %>% identity()
pdf(file = "path/image/T_cell/Tcell_CD8_UMAP_celltype.pdf", width = 6, height = 4)
DimPlot(Tcell_CD8_2, reduction = "umap", cols = c('#96b437', '#da93ab','#e58932'), label = TRUE, group.by = "celltype_T") + theme_bw() + theme(panel.grid=element_blank())
dev.off()

# Build data
data <- GetAssayData(Tcell_CD8_2, assay = 'RNA', slot = 'counts')
cell_metadata <- Tcell_CD8_2@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)

# UMAP dimensionality reduction
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')

## Import the consolidated umap coordinates from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Tcell_CD8_2, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

## Monocle3 cluster
cds <- cluster_cells(cds)

## Find trajectory
cds <- learn_graph(cds)
pdf(file = "path/image/T_cell/Tcell_CD8_UMAP_celltype_track.pdf", width = 6, height = 4)
plot_cells(cds, color_cells_by = "celltype_T", label_groups_by_cluster = FALSE, label_leaves = FALSE, cell_size = 0.1,
               label_branch_points = FALSE, label_cell_groups = F) + scale_colour_manual(values = c('#96b437', '#da93ab','#e58932'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_CD8_UMAP_group_track.pdf", width = 5, height = 4)
plot_cells(cds, color_cells_by = "group", label_groups_by_cluster = FALSE, label_leaves = FALSE, cell_size = 0.1,
               label_branch_points = FALSE, label_cell_groups = F) + scale_colour_manual(values = c('#af2337', '#2967a0'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/image/T_cell/Tcell_CD8_UMAP_level_track.pdf", width = 5, height = 4)
plot_cells(cds, color_cells_by = "level", label_groups_by_cluster = FALSE, label_leaves = FALSE, cell_size = 0.1,
               label_branch_points = FALSE, label_cell_groups = F) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Order cells by pseudotime
pdf(file = "path/image/T_cell/Tcell_CD8_UMAP_celltype_auxiliary_line.pdf", width = 4, height = 4)
plot_cells(cds, color_cells_by = "celltype_T", label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE) + geom_vline(xintercept = seq(8,9,0.25)) + geom_hline(yintercept = seq(-1,0,0.25))
dev.off()
embed <- data.frame(Embeddings(Tcell_CD8_2, reduction = "umap"))
embed <- subset(embed, umap_1 > 8.9 & umap_1 < 9.1 & umap_2 > -0.5 & umap_2 < -0.25)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)

pdf(file = "path/image/T_cell/Tcell_CD8_UMAP_celltype_pseudotime.pdf", width = 5.5, height = 4)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE, cell_size = 0.1, label_branch_points = FALSE)+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

## Find DEGs about pseudotime

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=100)

Track_genes_sig <- c("CCR7", "GZMA")
# Gene expression trend map
pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_celltype_T.pdf", width = 10, height = 4)
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 4)+ scale_colour_manual(values = c('#96b437', '#da93ab','#e58932'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_group.pdf", width = 10, height = 4)
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="group", 
                              min_expr=0.5, ncol = 4) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_level.pdf", width = 10, height = 4)
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="level", 
                              min_expr=0.5, ncol = 4) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

saveRDS(Tcell_CD8_2, file = "path/result/T_cell/Tcell_CD8_reUMAP.rds")
saveRDS(cds, file = "path/result/T_cell/Tcell_CD8_pseudotime.rds")

gene_list_Naive <- c("CCR7", "RPS14", "EEF1A1", "RPS28", "HLA-E", "RPL29", "RPL13", "RPS18", "RPS19", "RPS7")
gene_list_Cytotoxicity <- c("GZMA", "ACTB", "B2M", "TMSB4X", "PFN1", "HLA-C", "CFL1", "HLA-A", "NKG7", "HLA-B")

# Naive
pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_celltype_T_Naive.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Naive,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 6)+ scale_colour_manual(values = c('#96b437', '#da93ab','#e58932'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_group_Naive.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Naive,], color_cells_by="group", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_level_Naive.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Naive,], color_cells_by="level", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

# Cytotoxicity
pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_celltype_T_Cytotoxicity.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Cytotoxicity,], color_cells_by="celltype_T", 
                              min_expr=0.5, ncol = 6)+ scale_colour_manual(values = c('#96b437', '#da93ab','#e58932'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_group_Cytotoxicity.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Cytotoxicity,], color_cells_by="group", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#2967a0'))
dev.off()

pdf(file = "path/image/T_cell/Tcell_CD8_pseudotime_level_Cytotoxicity.pdf", width = 16, height = 4)
plot_genes_in_pseudotime(cds[gene_list_Cytotoxicity,], color_cells_by="level", 
                              min_expr=0.5, ncol = 6) + scale_colour_manual(values = c('#af2337', '#ecc342', '#2967a0'))
dev.off()

#### Assessing differentiation potential score by CytoTRACE2

#######Import seurat object###########
Tcell_CD8_2 <- readRDS("path/result/T_cell/Tcell_CD8_reUMAP.rds")
cytotrace2_result_sce <- cytotrace2(Tcell_CD8_2, 
                                is_seurat = TRUE, 
                                slot_type = "counts", 
                                species = 'human',
                                seed = 1234)
pdf(file = "path/image/T_cell/Tcell_CD8_CytoTRACE2.pdf", width = 5, height = 4)
FeaturePlot(cytotrace2_result_sce, "CytoTRACE2_Relative",pt.size = 0.1) + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                                      "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + theme_bw() + theme(panel.grid=element_blank())
dev.off()