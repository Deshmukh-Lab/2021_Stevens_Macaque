# "Balance between protective and pathogenic immune responses to pneumonia in the neonatal lung enforced by gut microbiota"
# Integrative analysis of pulmonary immune response to bacterial pathogens # 

# Stevens et al., in review, as of 6/24/2021 #

library(Seurat)
library(dplyr)
library(tidyverse)
library(scater)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

setwd("/Path_to_wd/")
# Read count matrices and add metadata 
FT1 <- Read10X(data.dir = "/Path_to_wd/FT1/filtered_feature_bc_matrix")
FT1 <-CreateSeuratObject(counts = FT1, min.cells = 1, min.features = 200) 
FT1@meta.data[, "id"] <- "FT1"

FT2 <- Read10X(data.dir = "/Path_to_wd/FT2/filtered_feature_bc_matrix")
FT2 <-CreateSeuratObject(counts = FT2, min.cells = 1, min.features = 200) 
FT2@meta.data[, "id"] <- "FT2"

NABX1 <- Read10X(data.dir = "/Path_to_wd/NABX1/filtered_feature_bc_matrix")
NABX1 <-CreateSeuratObject(counts = NABX1, min.cells = 1, min.features = 200) 
NABX1@meta.data[, "id"] <- "NABX1"
 
NABX2 <- Read10X(data.dir = "/Path_to_wd/NABX2/filtered_feature_bc_matrix")
NABX2 <-CreateSeuratObject(counts = NABX2, min.cells = 1, min.features = 200) 
NABX2@meta.data[, "id"] <- "NABX2"

ABX1 <- Read10X(data.dir = "/Path_to_wd/ABX1/filtered_feature_bc_matrix")
ABX1 <-CreateSeuratObject(counts = ABX1, min.cells = 1, min.features = 200) 
ABX1@meta.data[, "id"] <- "ABX1"

ABX2 <- Read10X(data.dir = "/Path_to_wd/ABX2/filtered_feature_bc_matrix")
ABX2 <-CreateSeuratObject(counts = ABX2, min.cells = 1, min.features = 200) 
ABX2@meta.data[, "id"] <- "ABX2"

# Merge all matrices and add metadata #
FT = merge(x = FT1, y = FT2)
FT@meta.data[, "treat"] <- "FT"
ABX = merge(x = ABX1, y = ABX2)
ABX@meta.data[, "treat"] <- "ABX"
NABX = merge(x = NABX1, y = NABX2)
NABX@meta.data[, "treat"] <- "NABX"
combined.A = merge(x = NABX, y = c(FT, ABX))
combined.A

# Regress out MT, Ribosomal (rps/rpl genes) #
combined.A <- PercentageFeatureSet(combined.A, pattern = "^MT-", col.name = "percent.mt")
combined.A <- PercentageFeatureSet(combined.A, pattern = "^RPS", col.name = "percent.rps")
combined.A <- PercentageFeatureSet(combined.A, pattern = "^RPL", col.name = "percent.rpl")
combined.A <- PercentageFeatureSet(combined.A, pattern = "^RNA\\d8S5", col.name = "percent.rrna")

# Use SCT to regress out the values # 

lung.list <- SplitObject(combined.A, split.by = "treat")

for (i in 1:length(lung.list)) {
  lung.list[[i]] <- SCTransform(lung.list[[i]], vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "percent.rrna", "nCount_RNA", "nFeature_RNA"), verbose = FALSE)
}
lung.features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000)
options(future.globals.maxSize= 2500000000000)
lung.list <- PrepSCTIntegration(object.list = lung.list, anchor.features = lung.features, 
                                verbose = FALSE)
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, normalization.method = "SCT", 
                                       anchor.features = lung.features, verbose = FALSE)
lung.integrated <- IntegrateData(anchorset = lung.anchors, normalization.method = "SCT", 
                                 verbose = FALSE)

lung.integrated <- RunPCA(lung.integrated, dims  = 1:20, verbose = FALSE)
lung.integrated <- RunUMAP(lung.integrated, dims = 1:20, verbose = FALSE)
lung.integrated <- FindNeighbors(lung.integrated, dims = 1:20, verbose = FALSE)
lung.integrated <- FindClusters(lung.integrated, resolution = 0.2) 
DimPlot(lung.integrated, label = T, split.by = "treat")

# find markers for every cluster compared to all remaining cells, report only the positive ones #
Immune.markers <- FindAllMarkers(lung.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(lung.integrated[[]])
Idents(lung.integrated) <- "seurat_clusters"
View(Immune.markers)

# Cluster 12 is bad cluster. Are there any genes other than poor quality genes that are upregulated in them?
bad12.markers <- FindMarkers(lung.integrated, ident.1 = "12")
View((bad12.markers[order(-bad.markers$avg_logFC),]))     
          
# Remove bad cluster 12
lung.integrated <- subset(lung.integrated, ident = c("12"), invert = T)
          
# Map to reference using Seurat Web app 

DefaultAssay(lung.integrated) <-"RNA"
wnn.1.obj <- DietSeurat(lung.integrated, assays = "RNA")
saveRDS(wnn.1.obj, file = "/Path_to_save/wnn.1.obj.rds")

# Use Azimuth to map to PBMC reference on HCA and download the results
predictions <- read.delim(file = "Integrated.Object/azimuth_pred.tsv", row.names = 1)
lung.integrated <- AddMetaData(object = lung.integrated,metadata = predictions)
projected.umap <- readRDS(file = "Integrated.Object/azimuth_umap.Rds")
lung.integrated <- lung.integrated[, Cells(projected.umap)]
lung.integrated[['umap.proj']] <- projected.umap
head(lung.integrated[[]])
          
# Change ids of clusters 
lung.integrated <- RenameIdents(lung.integrated, "0" ="Neutrophil", "1" = "T cell", "2" = "Alveolar Macrophage", "3" = "Neutrophil", "4" = "Intersitial Macrophage", "5" = "B cell", "6" = "NK cell", "7" = "Platelet", "8" = "Plasmablast", "9" = "Proliferating NK Cells", "10" = "cDC", "12" = "pDC")
new.cluster.ids <- c("Neutrophil", "T cell", "Alveolar Macrophage", "Intersitial Macrophage", "B cell", "NK cell", "Platelet", "Plasmablast", "Proliferating NK Cells", "cDC1", "cDC2","pDC")
names(new.cluster.ids) <- levels(lung.integrated)                     
lung.integrated <- RenameIdents(lung.integrated, new.cluster.ids)
head(lung.integrated[[]])
saveRDS(lung.integrated, file = "Path_to_save/lung.integrated.Rds")

# For cellchat analysis, will restrict to major cell.types, for example 
# Neutrophil, T cells, Macrophage Subsets, NK cells, B cells, cDC1 and pDC

Idents(lung.integrated) <- "seurat_clusters"
Cell.Chat <- subset(lung.integrated, idents = c("0","1","2","3","4","5","9","11"))
Idents(Cell.Chat) <- "exp"
Cell.Chat.FecalTrans <- subset(Cell.Chat, idents = "FT")
saveRDS(Cell.Chat.FecalTrans, file = "Path_to_save/Cell.Chat.FecalTrans.rds")
Cell.Chat.ABX <- subset(Cell.Chat, idents = "ABX")

Cell.Chat.NABX <- subset(Cell.Chat, idents = "NABX")
saveRDS(Cell.Chat.NABX, file = "Path_to_save/Cell.Chat.NABX.rds")

# PROCEED WITH NEUTROPHIL ANALYSIS # 

Neutrophil <- subset(lung.integrated, idents = "Neutrophil")
DefaultAssay(Neutrophil) = "RNA"

Neutrophil.list <- SplitObject(Neutrophil, split.by = "treat")

for (i in 1:length(Neutrophil.list)) {
  Neutrophil.list[[i]] <- NormalizeData(Neutrophil.list[[i]], verbose = FALSE)
  Neutrophil.list[[i]] <- FindVariableFeatures(Neutrophil.list[[i]], selection.method = "vst", nfeatures = 3000,verbose = FALSE)}

Neutrophil.anchors <- FindIntegrationAnchors(object.list = Neutrophil.list, dims = 1:10, k.filter = 10)
Neutrophil <- IntegrateData(anchorset = Neutrophil.anchors, dims = 1:20)
Neutrophil <- ScaleData(Neutrophil, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"))
Neutrophil <- RunPCA(Neutrophil, npcs = 25, verbose = FALSE)
Neutrophil <- RunUMAP(Neutrophil, reduction = "pca", dims = 1:15)
Neutrophil <- FindNeighbors(Neutrophil, reduction = "pca", dims = 1:15)
Neutrophil <- FindClusters(Neutrophil, resolution = 0.12)
DimPlot(Neutrophil, split.by = "treat")
DimPlot(Neutrophil)
saveRDS(Neutrophil, file = "Path_to_save/Neutrophil_Manuscript.rds")

# Find markers for every cluster compared to all remaining cells, report only the positive ones
Cluster.markers <- FindAllMarkers(Neutrophil, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.30)
View(Cluster.markers)
write.table(Cluster.markers, file = "Path_to_save/Neutrophil/All.Cluster.markers.tsv", sep="\t")
Cluster.top.markers <-Cluster.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) %>% print(n=250)
View(Cluster.top.markers)
write.table(Cluster.top.markers, file = "Path_to_save/Top.Cluster.markers.tsv", sep="\t")
DefaultAssay(Neutrophil) <- "RNA"
Neutrophil.scale <-ScaleData(Neutrophil, features = Cluster.top.markers$gene)
# Remove duplicates for making heatmaps
List <- Cluster.top.markers[!duplicated(Cluster.top.markers$gene), ]
List <- column_to_rownames(List, var = "gene")
List <- rownames(List)
library(viridis)
DoHeatmap(subset(Neutrophil.scale), features = List, size = 4) + scale_fill_viridis_b() + NoLegend()

# To make figure specific for Aged neutrophil Cluster
FeaturePlot(Neutrophil, features = c("CXCR2", "CXCR4", "CD63", "SELL"), label = TRUE, min.cutoff = "q05", max.cutoff = "q95", split.by = "treat")
FeaturePlot(Neutrophil, features = c("G6PD", "RNF13", "PMSA7", "GADD45B"), label = TRUE, min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(Neutrophil, features = c("TRAF1", "TNFAIP6", "MAP4K4", "TNFSF14"), label = TRUE, min.cutoff = "q05", max.cutoff = "q95", pt.size = 1, cols = c("grey", "seagreen"))

# To make figure specific for naive neutrophils 
FeaturePlot(Neutrophil, features = "OLFM4", label = F, min.cutoff = "q05", max.cutoff = "q95", pt.size = 2, cols = c("grey", "slateblue"))
FeaturePlot(Neutrophil, features = "GCA", label = F, min.cutoff = "q50", max.cutoff = "q95", pt.size = 2, cols = c("grey", "slateblue"))
FeaturePlot(Neutrophil, features = "SELL", label = F, min.cutoff = "q60", max.cutoff = "q85", pt.size = 2, cols = c("grey", "slateblue"))

# To make figure specific for activated neutrophils 
FeaturePlot(Neutrophil, features = c("AIF1", "IRF1", "TNFAIP3", "PTAFR"), label = TRUE, min.cutoff = "q05", max.cutoff = "q95", pt.size = 2, cols = c("grey", "salmon"))
FeaturePlot(Neutrophil, features = "AIF1", label = F, min.cutoff = "q20", max.cutoff = "q95", pt.size = 2, cols = c("grey", "salmon"))
FeaturePlot(Neutrophil, features = "IRF1", label = F, min.cutoff = "q5", max.cutoff = "q95", pt.size = 2, cols = c("grey", "salmon"))
FeaturePlot(Neutrophil, features = "CSF1", label = TRUE, min.cutoff = "q5", max.cutoff = "q95", pt.size = 2, cols = c("grey", "salmon"))
FeaturePlot(Neutrophil, features = "NFKB1", label = TRUE, min.cutoff = "q5", max.cutoff = "q95", pt.size = 2, cols = c("grey", "salmon"))
FeaturePlot(Neutrophil, features = "CSF3R", label = F, min.cutoff = "q5", max.cutoff = "q95", pt.size = 2, cols = c("grey", "salmon"))
FeaturePlot(Neutrophil, features = "LGALS3", label = TRUE, min.cutoff = "q5", max.cutoff = "q50", pt.size = 2, cols = c("grey", "salmon"))

# To make figure specific for Aged neutrophil Cluster

FeaturePlot(Neutrophil, features = "CTSB", label = FALSE, min.cutoff = "q10", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "CD63", label = FALSE, min.cutoff = "q05", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "GADD45B", label = FALSE, min.cutoff = "q05", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "PSMA7", label = FALSE, min.cutoff = "q05", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "G6PD", label = FALSE, min.cutoff = "q01", max.cutoff = "q40", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "RNF13", label = FALSE, min.cutoff = "q10", max.cutoff = "q40", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "CXCR2", label = FALSE, min.cutoff = "q01", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "CXCR4", label = FALSE, min.cutoff = "q01", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "CD63", label = FALSE, min.cutoff = "q20", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "SELL", label = FALSE, min.cutoff = "q20", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "IL1B", label = FALSE, min.cutoff = "q60", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "CXCL8", label = FALSE, min.cutoff = "q60", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "NFKB1", label = FALSE, min.cutoff = "q10", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "HK1", label = FALSE, min.cutoff = "q05", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "CD274", label = FALSE, min.cutoff = "q10", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "PADI4", label = FALSE, min.cutoff = "q10", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)
FeaturePlot(Neutrophil, features = "PFKL", label = TRUE, min.cutoff = "q01", max.cutoff = "q95", cols = c("grey", "seagreen"), pt.size = 2)

# Heatmap of NADPH oxidase comples related genes (Not shown in the manuscript) # 
NADPH.genes <-c("CYBB", "CYBA", "RAC2", "RAC1", "NCF2", "NCF1", "NCF4")
DoHeatmap(subset(Neutrophil.all.scale, downsample = 1000), features = NADPH.genes, size = 2)
C2.NFKB <-c("ID2", "CXCL1", "S100A4", "GADD45A", "C5AR1", "SELL", "SERPIB1", "TXNIP", "IRF1", "CCL2", "CD69")
DoHeatmap(subset(Neutrophil.all.scale, downsample = 1000), features = C2.NFKB, size = 2)

# Heatmap of glycolysis related genes
Energetics <-c("GPI1", "ALDOA", "GAPDH", "PKM", "ENO1", "PGAM1", "TPI1", "PGK1", "HK1", "HK2", "HKDC1", "PFKL", "PFKM", "PGAM2", "ENO2", "PKLR")
DoHeatmap(subset(Neutrophil.scale, downsample = 1000), features = Energetics, size = 2)
FeaturePlot(Neutrophil, features = Energetics)

# Barchart of change in frequency of cells #

table(Neutrophil@meta.data$seurat_clusters, Neutrophil@meta.data$treat)
View(Table)
pt <- table(Neutrophil@meta.data$seurat_clusters, Neutrophil@meta.data$treat)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 25) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_color_manual(color) +
  theme(legend.title = element_blank())

# Calculate Module Scores # based on https://www.biorxiv.org/content/10.1101/2020.12.18.423363v1.supplementary-material
Genes <-rownames(Neutrophil)
Module <- read.csv2(file = "path_to_module_score(SupplementalTable12)", sep = ",")
View(Module)

ISG <- Module[,1]
ISG  <- ISG[ISG  %in% Genes]
Phagocytosis <- Module[,2]
Phagocytosis  <- Phagocytosis[Phagocytosis  %in% Genes]
Degranulation <- Module[,3]
Degranulation  <- Degranulation[Degranulation  %in% Genes]
PDL1_Neutrophil <- Module[,4]
PDL1_Neutrophil  <- PDL1_Neutrophil[PDL1_Neutrophil  %in% Genes]
Sepsis <-Module[,5]
Sepsis  <- Sepsis[Sepsis  %in% Genes]

DefaultAssay(Neutrophil) <-"RNA"
Neutrophil <- AddModuleScore(Neutrophil, features = list(ISG), ctrl.size = 5, name = 'ISGgenes')
Neutrophil <- AddModuleScore(Neutrophil, features = list(Phagocytosis), ctrl.size = 5, name ='Phagoc')
Neutrophil <- AddModuleScore(Neutrophil, features = list(Degranulation), ctrl.size = 5, name ='Degranulation')
Neutrophil <- AddModuleScore(Neutrophil, features = list(PDL1_Neutrophil), ctrl.size = 5, name ='PDL1_Neutrophil')
Neutrophil <- AddModuleScore(Neutrophil, features = list(Sepsis), ctrl.size = 5, name ='Sepsis')

# Visualized Neutrophil metadata slots
head(Neutrophil[[]])
Neutrophil[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophil, vars = 'ISGgenes1')))
Neutrophil[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophil, vars = 'Phagoc1')))
Neutrophil[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophil, vars = 'Degranulation1')))
Neutrophil[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophil, vars = 'PDL1_Neutrophil1')))
Neutrophil[['module']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophil, vars = 'Sepsis1')))

# Visualize Module scores
VlnPlot(object = Neutrophil, features = 'ISGgenes1', assay = 'module', slot = 'data', pt.size = 0.25) + NoLegend()
VlnPlot(object = Neutrophil, features = 'Phagoc1', assay = 'module', slot = 'data', pt.size = 0.25) + NoLegend()
VlnPlot(object = Neutrophil, features = 'Degranulation1', assay = 'module', slot = 'data', pt.size = 0.25) + NoLegend()
VlnPlot(object = Neutrophil, features = 'PDL1_Neutrophil1', assay = 'module', slot = 'data', pt.size = 0.25) + NoLegend()
VlnPlot(object = Neutrophil, features = 'Sepsis1', assay = 'module', slot = 'data', pt.size = 0.25) + NoLegend()

# Visualize using Feature Plot
A <- FeaturePlot(object = Neutrophil, features = 'ISGgenes1', min.cutoff = "q10", max.cutoff = "q50", pt.size = 1, cols = c("grey", "#440154FF")) + NoLegend()
B <- FeaturePlot(object = Neutrophil, features = 'Phagoc1', min.cutoff = "q20", max.cutoff = "q95", pt.size = 1, cols = c("grey", "#2D708EFF")) + NoLegend()
C <- FeaturePlot(object = Neutrophil, features = 'Degranulation1', min.cutoff = "q20", max.cutoff = "q95", pt.size = 1, cols = c("grey","#B8627DFF")) + NoLegend()
D <- FeaturePlot(object = Neutrophil, features = 'PDL1_Neutrophil1', min.cutoff = "q20", max.cutoff = "q95", pt.size = 1, cols = c("grey", "#403891FF")) + NoLegend()
E <- FeaturePlot(object = Neutrophil, features = 'Sepsis1', min.cutoff = "q10", max.cutoff = "q95", pt.size = 1, cols = c("grey", "#F68F46FF")) + NoLegend()
plot_grid(A,B,C,D,E)

# Pseudotime analysis
DefaultAssay(Neutrophil) <- "integrated"
cds <- as.cell_data_set(Neutrophil)
cds <- cluster_cells(cds)
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", group_label_size = 3, show_trajectory_graph = TRUE)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = T)

cds <- estimate_size_factors(cds)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "seurat_clusters", label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 2, show_trajectory_graph = T)
rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, cell_size = 1.5)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
Idents(integrated.sub) <- "treat"
ABX <- subset(integrated.sub, ident = "ABX")
FeaturePlot(ABX, "monocle3_pseudotime", pt.size = 1.5)
NABX <- subset(integrated.sub, ident = "NABX")
FeaturePlot(NABX, "monocle3_pseudotime", pt.size = 1.5)
FeaturePlot(integrated.sub, "monocle3_pseudotime", pt.size = 1.5)

# WORKFLOW for identifying genes which change as function of pseudotime #

cds <- estimate_size_factors(cds)
cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
head(cds_graph_test_results)
deg_ids <-rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.1))
deg_ids <-rownames(subset(cds_graph_test_results, q_value < 0.1))
view(deg_ids)
write.table(deg_ids, file = "Path_to_save/Neutrophil/deg_ids.csv", sep =",")
plot_cells(cds, genes = head(deg_ids), show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)

# Collect the trajectory variable genes into modules

gene_modules <- find_gene_modules(cds[deg_ids,], resolution = c(10^seq(-6,-1)))
View(gene_modules)
write.table(gene_modules, file = "Path_to_save/Neutrophil/genes_changing_as_pseudotime.csv", sep=",")

#Visualize the GENES that change as function of pseudotime #

plot_cells(cds, genes=gene_modules,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE, cell_size = 1.5)

module4 <- filter(gene_modules, module ==4)
plot_cells(cds, genes=module4,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE, cell_size = 2)

module3 <- filter(gene_modules, module ==3)
plot_cells(cds, genes=module3,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE, cell_size = 2)

Interesting_genes <- c("CXCR2", "SELL", "GCA", "CXCL2")
Interesting_genes.2 <- c("CXCR2", "CXCR4", "SELL", "GCA", "CXCL8", "BCL2", "OLFM4", "MMP8", "LTF", "CD47", "CXCL2")
Interesting_cds <- cds[rowData(cds)$gene_short_name %in% Interesting_genes.2]

plot_genes_in_pseudotime(Interesting_cds, min_expr= 0.01, cell_size = 0.75)

# SCENIC analysis 

## Load packages if not already done ##
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(ComplexHeatmap)
library(R2HTML)
library(Rtsne)
library(Seurat)
library(SCopeLoomR)
library(grid)
library(data.table)

# set working directory
setwd("Path_to_save/Neutrophil/SCENIC")

# Make SingleCellExperiment object as it is easier to export counts and barcodes
Neutrophil.sce <- as.SingleCellExperiment(Neutrophil)
expMat <- counts(Neutrophil.sce)
cellInfo <- colData(Neutrophil.sce)
cellInfo <- data.frame(seuratCluster=Idents(Neutrophil))

dbDir <-"Path_to_db/SCENIC"

# set organism 
org <- "hgnc"
# set name
myDatasetTitle <- "Neutrophil"
data("defaultDbNames")
dbs <-defaultDbNames[[org]]

# set options
scenicOptions <- initializeScenic(org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle, nCores = 10)
# save to use later
saveRDS(scenicOptions, file = "Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")

# GeneFilter/Selection # 

genesKept <- geneFiltering(expMat, scenicOptions = scenicOptions)
dim(as.matrix(genesKept))
expMat_filtered <- expMat[genesKept, ]
dim(expMat_filtered)

# Build a co-expression network # 
# Calculate correlation to detect positive and negtaive associations
runCorrelation(as.matrix(expMat_filtered), scenicOptions)
expMat_filtered_log <- log2(expMat_filtered+1) 
View(expMat_filtered_log)

# Run GRN to detect regulons (TF and cis target gene modules)
runGenie3(as.matrix(expMat_filtered_log), scenicOptions) 
saveRDS(scenicOptions, file = "Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")

# For next 4 steps, load the scenicOptions.rds file before each operation.
scenicOptions <- readRDS("Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")

# Build and score the GRN # 

expMat <-readRDS(file ="/Path_to_data/expMat")
expMat_log <- log2(expMat+1)
View(as.matrix(expMat_log))

# Step 1 - runSCENIC_1 
scenicOptions <- readRDS("Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# Step 2 - runSCENIC_2
scenicOptions <- readRDS("Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

# Step 3 - Analyze network activity in each individual cells using AU Cell.
scenicOptions <- readRDS("Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expMat) # I used expMat instead of expMat_log, this worked. 

# Export 
scenicOptions <- readRDS("Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo"))

# Regulators of clusters # 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
saveRDS(regulonActivity_byCellType_Scaled, file = "Path_to_save/Neutrophil/SCENIC/int/regulonActivity_byCellType_Scaled")
View(regulonActivity_byCellType_Scaled)
write.table(regulonActivity_byCellType_Scaled, file = "Path_to_save/Neutrophil/SCENIC/int/regulonActivity_by_Cluster_Scaled.tsv", sep = '\t')

# Save RDS 
saveRDS(scenicOptions, file="Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")

# Plot the heatmap
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity") 

# List of top regulators 
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
saveRDS(topRegulators, file = "Path_to_save/Neutrophil/SCENIC/int/topRegulators")
write.csv2(topRegulators, file = "Path_to_save/Neutrophil/SCENIC/int/topRegulators.csv")
# Also save as html Top_Regulators

# Visualizing the regulon activities on embedding/trajectories calculated with other methods.
#Load dependent libraries if not already done#
library(Seurat)
library(viridis)
scenicOptions <- readRDS("Path_to_save/Neutrophil/SCENIC/int/scenicOptions.rds")
dr_coords <- Embeddings(treatment.integrated.E, reduction="umap")

# Plot activity for C3 cluster
C3.tfs <- c("ATF6","IRF9","STAT1","JUNB","ELK3","IRF1","IKZF1","FLI1", "IRF3","ETV6","RARA","CTCF")
# Plot for C2 cluster 
C2.tfs <- c("SPI1","RUNX1","BCLAF1","KDM5A","CPEB1","FLI1","NFIL3")
# Plot for C1 cluster
C1.tfs <- c("SP4","HDAC2","FOSL","IRF2","NFKB2","CREM","REL", "RELB","MAX","EGR2")
par(mar = c(2,2,2,2), mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, C3.tfs), plots = "AUC", cex = 0.5)
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, C2.tfs), plots = "AUC", cex = 0.5)
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, C1.tfs), plots = "AUC", cex = 0.5)


## Alveolar Macrophage analsysis ## 
# Load and visualize data if starting from this point #

lung.integrated <- readRDS(file = "/Path_to_data/lung.integrated.rds")
Idents(lung.integrated) <-"seurat_clusters"

# Subset Lung.intregarted for NABX and ABX
Idents(lung.integrated) <-"exp"
obj <- subset(lung.integrated, idents = c("NABX", "ABX"))
obj # 104604 features across 9220 samples within 6 assays 
Idents(obj) <-"seurat_clusters"

# Subset Alveolar Macrophages for further analysis#  --
A.Macs <- subset(obj, idents = "Alevolar Macrophage")
A.Macs # 1104604 features  795 samples within 6 assays
DimPlot(A.Macs, reduction = "umap", split.by = "exp", label = TRUE)

# Begin processing data #

options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(A.Macs) <-"RNA"
A.Macs.list <- SplitObject(A.Macs, split.by = "exp")
for (i in 1:length(A.Macs.list)) {
  A.Macs.list[[i]] <- SCTransform(A.Macs.list[[i]], verbose = FALSE)
}

A.Macs.features <- SelectIntegrationFeatures(object.list = A.Macs.list, nfeatures = 2000)
A.Macs.list <- PrepSCTIntegration(object.list = A.Macs.list, anchor.features = A.Macs.features, 
                                  verbose = FALSE)
A.Macs.anchors <- FindIntegrationAnchors(object.list = A.Macs.list, normalization.method = "SCT", 
                                         anchor.features = A.Macs.features, verbose = FALSE)
A.Macs.integrated <- IntegrateData(anchorset = A.Macs.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)
A.Macs.integrated <- RunPCA(A.Macs.integrated, verbose = FALSE)
A.Macs.integrated <- RunUMAP(A.Macs.integrated, dims = 1:10)
A.Macs.integrated <- FindNeighbors(A.Macs.integrated, dims = 1:10)
A.Macs.integrated  <- FindClusters(A.Macs.integrated, resolution = 0.2)
plots <- DimPlot(A.Macs.integrated, split.by = "exp")
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))
saveRDS(A.Macs.integrated, file = "/Path_to_save/A.Macs.integrated.rds")

A.Macs.markers <- FindAllMarkers(A.Macs.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(A.Macs.markers) # 965 genes
write.table(A.Macs.markers, file = "/Path_to_save/A.Macs.markers1.csv", sep = ",")

#report top 30 markers for each group, by avglogFC
A.Macs.top.markers <-A.Macs.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) %>% print(n=250)
View(A.Macs.top.markers) # 90 genes
write.table(A.Macs.top.markers, file = "/Path_to_save/A.Macs.top.markers.csv", sep = ",")

# To generate figure 
DimPlot(A.Macs.integrated, reduction = "umap", label = TRUE, repel = TRUE, split.by = "exp")

# Extract raw counts and metadata for Mac Spectrum
counts <- A.Macs.integrated@assays$RNA@counts
write.table(counts, file = "/Path_to_save/A.Mac_spectrum.csv", sep = ",", row.names = T)
metadata <- A.Macs.integrated@meta.data
metadata  <- metadata [, c("seurat_clusters", "exp")]
write.table(metadata, file = "/Path_to_save/A.Mac_spectrum.metadata.csv", sep = ",", row.names = T)

# Complete MacSprectum Analysis online at http://155.37.255.130:3838/macspec/

# Download MPI/AMDI table 

Macs.Spectrum <- read_csv(file = "/Path_to_read/MPI_AMDI_table.csv")
# Set Feature as grouping varibale 
Macs.Spectrum$Feature <- as.factor(Macs.Spectrum$Feature)

B1 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) 
B2 <- B1 + geom_point(aes(color = Feature), alpha = 0.5)
B3 <- B2 + geom_density2d(color = "grey", linetype =1)
B3
# Set Exp as grouping varibale 
Macs.Spectrum$Exp <- as.factor(Macs.Spectrum$Exp)
A1 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) + facet_grid(~Exp)
A2 <- A1 + geom_point(aes(color = Exp), alpha = 0.5, size = 1.5)
A3 <- A2 + geom_density2d(color = "grey", linetype =1)
A3

A4 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) + facet_grid(~Feature)
A5 <- A4 + geom_point(aes(color = Exp), alpha = 0.5, size = 1.5)
A6 <- A5 + geom_density2d(color = "grey", linetype =1)
A6

A4 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) + facet_grid(~Feature)
A5 <- A4 + geom_point(aes(color = Feature), alpha = 0.5, size = 1.5, color = c("salmon", "green", "blue"))
A6 <- A5 + geom_density2d(color = "grey", linetype =1)
A6


# Now by Clusters
# Subset by treatment groups 
C0 <- subset(Macs.Spectrum, Feature=="0")
C0$Exp <- as.factor(C0$Exp)
C1 <- subset(Macs.Spectrum, Feature=="1")
C1$Feature <- as.factor(C1$Feature)
C2 <- subset(Macs.Spectrum, Feature=="2")
C2$Feature <- as.factor(C2$Feature)

C0.1 <- ggplot(C0, aes(x= MPI, y = AMDI)) + xlim(-10,10) + ylim(-30,30) + geom_point(alpha = 0.5, size = 1.5, color = "salmon") + geom_density2d(color = "grey", linetype =1) + theme_classic()
C1.1 <- ggplot(C1, aes(x= MPI, y = AMDI)) + xlim(-10,10) + ylim(-30,30) + geom_point(alpha = 0.5, size = 1.5, color = "limegreen") + geom_density2d(color = "grey", linetype =1) + theme_classic()
C2.1 <- ggplot(C2, aes(x= MPI, y = AMDI)) + xlim(-10,10) + ylim(-30,30) + geom_point(alpha = 0.5, size = 1.5, color = "steelblue1") + geom_density2d(color = "grey", linetype =1) + theme_classic()
C0.1+C1.1+C2.1

# Module score as done for neutrophils#
Genes <-rownames(A.Macs.integrated)
View(Genes)
Module <- read.csv2(file = "/Critical_Scripts/Module_Score/Cytokine_Inf_genes2.csv", sep = ",")

Mild.RDS <- Module[,1]
Mild.RDS  <- Mild.RDS[Mild.RDS  %in% Genes]

Severe.RDS <- Module[,2]
Severe.RDS  <- Severe.RDS[Severe.RDS  %in% Genes]

DefaultAssay(AM) <-"RNA"
A.Macs.integrated <- AddModuleScore(A.Macs.integrated, features = list(Mild.RDS), ctrl.size = 5, name = 'Mild.RDS', nbin = 20)
A.Macs.integrated <- AddModuleScore(A.Macs.integrated, features = list(Severe.RDS), ctrl.size = 5, name ='Severe.RDS', nbin = 20)

A.Macs.integrated[['module']] <- CreateAssayObject(data = t(x = FetchData(object = A.Macs.integrated, vars = 'Mild.RDS1'))) 
A.Macs.integrated[['module']] <- CreateAssayObject(data = t(x = FetchData(object = AA.Macs.integrated, vars = 'Severe.RDS1')))

# Plot heatmaps for for marker genes # 
DefaultAssay(A.Macs.integrated) <- "RNA"
A.Macs.scale <-ScaleData(A.Macs.integrated, features = A.Macs.top.markers$gene)

# Remove duplicates for making heatmaps #
List <- A.Macs.top.markers[!duplicated(A.Macs.top.markers$gene), ]
List <- column_to_rownames(List, var = "gene")
List <- rownames(List)

# Plot the heatmap 
library(viridis)
DoHeatmap(subset(A.Macs.scale, downsample = 100), features = List, size = 3, angle = 90) + scale_fill_viridis_b() + NoLegend()
# Saved as PDF.

# Plot DimPlot 
P1 <- DimPlot(A.Macs.integrated, reduction = "umap", pt.size = 2, label = FALSE, split.by = "exp")+NoLegend()
P1 <- DimPlot(A.Macs.integrated, reduction = "umap", pt.size = 2, label = FALSE) + NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1

A2 <-LabelClusters(A1, id = "ident", color = unique(ggplot_build(p1)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)
A1
A2 # (600 X 600)
saveRDS(A.Macs.integrated, file = "/Path_to_save/A.Macs.integrated.rds")

# Plots for changes with exp # 
P1 <- DimPlot(A.Macs.integrated, reduction = "umap", pt.size = 1.5, split.by = "exp")+NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1 # (600 X 600)
saveRDS(A.Macs.integrated, file = "/Path_to_save/A.Macs.integrated.rds")


# Make plots for cell distribution # 
pt <- table(Idents(A.Macs.integrated), A.Macs.integrated$exp)
pt <- as.data.frame(pt)
pt <- read_delim(file = "/Path_to_read/pt.csv", delim = ",")
pt$Var1 <- as.character(pt$Var1)
color <-c("#B79F00","#F564E3","#F8766D","#00BA38","#619CFF","#00BFC4")

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 25) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_color_manual(color) +
  theme(legend.title = element_blank())

# # SCENIC analysis for Alveolar Macropahges

library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(ComplexHeatmap)
library(R2HTML)
library(Rtsne)
library(Seurat)
library(SCopeLoomR)
library(grid)
library(data.table)
# set working directory
setwd("Path_to_save/A.Macs/SCENIC")

# Make SingleCellExperiment object as it is easier to export counts and barcodes
Neutrophil.sce <- as.SingleCellExperiment(A.Macs.integrated)
expMat <- counts(A.Macs.integrated.sce)
cellInfo <- colData(A.Macs.integrated.sce)
cellInfo <- data.frame(seuratCluster=Idents(A.Macs.integrated))
dbDir <-"Path_to_db/SCENIC"

# set organism 
org <- "hgnc"
# set name
myDatasetTitle <- "A.Macs"
data("defaultDbNames")
dbs <-defaultDbNames[[org]]

# set options
scenicOptions <- initializeScenic(org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle, nCores = 10)
# save to use later
saveRDS(scenicOptions, file = "Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")

# GeneFilter/Selection # 
genesKept <- geneFiltering(expMat, scenicOptions = scenicOptions)
dim(as.matrix(genesKept))
expMat_filtered <- expMat[genesKept, ]
dim(expMat_filtered)

# Build a co-expression network # 
runCorrelation(as.matrix(expMat_filtered), scenicOptions)
expMat_filtered_log <- log2(expMat_filtered+1) 
View(expMat_filtered_log)

# Run GRN to detect regulons (TF and cis target gene modules)
runGenie3(as.matrix(expMat_filtered_log), scenicOptions) 
saveRDS(scenicOptions, file = "Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")

# For next 4 steps, load the scenicOptions.rds file before each operation.
scenicOptions <- readRDS("Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")

# Build and score the GRN # 

expMat <-readRDS(file ="/Path_to_data/expMat")
expMat_log <- log2(expMat+1)
View(as.matrix(expMat_log))

# Step 1 - runSCENIC_1 
scenicOptions <- readRDS("Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# Step 2 - runSCENIC_2
scenicOptions <- readRDS("Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

# Step 3 - Analyze network activity in each individual cells using AU Cell.
scenicOptions <- readRDS("Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expMat)

# Export 
scenicOptions <- readRDS("Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo"))

# Regulators of clusters # 
scenicOptions <- readRDS("v")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
saveRDS(regulonActivity_byCellType_Scaled, file = "Path_to_save/A.Macs/SCENIC/int/regulonActivity_byCellType_Scaled")
View(regulonActivity_byCellType_Scaled)
write.table(regulonActivity_byCellType_Scaled, file = "Path_to_save/A.Macs/SCENIC/int/regulonActivity_by_Cluster_Scaled.tsv", sep = '\t')
# Save RDS 
saveRDS(scenicOptions, file="Path_to_save/NA.Macs/SCENIC/int/scenicOptions.rds")

# Plot the heatmap
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity") 
# List of top regulators 
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
saveRDS(topRegulators, file = "Path_to_save/A.Macs/SCENIC/int/topRegulators")
write.csv2(topRegulators, file = "Path_to_save/A.Macs/SCENIC/int/topRegulators.csv")
# Also save as html Top_Regulators

# Visualizing the regulon activities on embedding/trajectories calculated with other methods.

#Load dependent libraries if not already done#
library(Seurat)
library(viridis)
scenicOptions <- readRDS("Path_to_save/A.Macs/SCENIC/int/scenicOptions.rds")
dr_coords <- Embeddings(treatment.integrated.E, reduction="umap")

# Plot activity for C3 cluster
C3.tfs <- c("ATF6","IRF9","STAT1","JUNB","ELK3","IRF1","IKZF1","FLI1", "IRF3","ETV6","RARA","CTCF")
# Plot for C2 cluster 
C2.tfs <- c("SPI1","RUNX1","BCLAF1","KDM5A","CPEB1","FLI1","NFIL3")
# Plot for C1 cluster
C1.tfs <- c("SP4","HDAC2","FOSL","IRF2","NFKB2","CREM","REL", "RELB","MAX","EGR2")
par(mar = c(2,2,2,2), mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, C3.tfs), plots = "AUC", cex = 0.5)
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, C2.tfs), plots = "AUC", cex = 0.5)
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, C1.tfs), plots = "AUC", cex = 0.5)

## Interstitial Macrophage analysis ##
Idents(lung.integrated) <-"seurat_clusters"
DimPlot(lung.integrated, reduction = "umap", split.by = "treat", label = TRUE)

# Subset Lung.integrated for NABX and ABX
Idents(lung.integrated) <-"treat"
obj <- subset(lung.integrated, idents = c("NABX", "ABX"))
obj # 104604 features across 9220 samples within 6 assays 
Idents(obj) <-"seurat_clusters"
DimPlot(obj, reduction = "umap", label = TRUE)

# Subset Interstitial Macrophages for further analysis#  
I.Macs <- subset(obj, idents = "Intersitial Macrophage")
I.Macs # 1104604 features 984 samples within 6 assays
DimPlot(I.Macs, reduction = "umap", split.by = "treat", label = TRUE)

# Begin data processing #
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(I.Macs) <-"RNA"
I.Macs.list <- SplitObject(I.Macs, split.by = "treat")

for (i in 1:length(I.Macs.list)) {
  I.Macs.list[[i]] <- SCTransform(I.Macs.list[[i]], verbose = FALSE)
}

I.Macs.features <- SelectIntegrationFeatures(object.list = I.Macs.list, nfeatures = 2000)

I.Macs.list <- PrepSCTIntegration(object.list = I.Macs.list, anchor.features = I.Macs.features, 
                                  verbose = FALSE)

I.Macs.anchors <- FindIntegrationAnchors(object.list = I.Macs.list, normalization.method = "SCT", 
                                         anchor.features = I.Macs.features, verbose = FALSE)

I.Macs.integrated <- IntegrateData(anchorset = I.Macs.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)

I.Macs.integrated <- RunPCA(I.Macs.integrated, verbose = FALSE)
I.Macs.integrated <- RunUMAP(I.Macs.integrated, dims = 1:15)
I.Macs.integrated <- FindNeighbors(I.Macs.integrated, dims = 1:15)
I.Macs.integrated  <- FindClusters(I.Macs.integrated, resolution = 0.15)

plots <- DimPlot(I.Macs.integrated, split.by = "treat")
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))

saveRDS(I.Macs.integrated, file = "/Path_to_save/I.Macs.integrated.rds")

#Find markers for every cluster compared to all remaining cells, report only the positive ones
I.Macs.markers <- FindAllMarkers(I.Macs.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(I.Macs.markers) # 896 genes
write.table(I.Macs.markers, file = "/Path_to_save/I.Macs.markers.csv", sep = ",")

#report top 30 markers for each group, by avglogFC
I.Macs.top.markers <-I.Macs.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) %>% print(n=250)
View(I.Macs.top.markers)
write.table(I.Macs.top.markers, file = "/Path_to_save/I.Macs.top.markers.csv", sep = ",")

# To generate figure 
DimPlot(I.Macs.integrated, reduction = "umap", label = TRUE, repel = TRUE, split.by = "exp")

# Extract raw counts and metadata for Mac Spectrum
counts <- I.Macs.integrated@assays$RNA@counts
write.table(counts, file = "/Path_to_save/Mac_spectrum.csv", sep = ",", row.names = T)
metadata <- I.Macs.integrated@meta.data
metadata  <- metadata [, c("seurat_clusters", "exp")]
write.table(metadata, file = "/Path_to_save/Mac_spectrum.metadata.csv", sep = ",", row.names = T)

# Complete MacSprectum Analysis online at http://155.37.255.130:3838/macspec/
# Download MPI/AMDI table 

Macs.Spectrum <- read_csv(file = "/Path_to_read/MPI_AMDI_table.csv")
# Set Feature as grouping varibale 
Macs.Spectrum$Feature <- as.factor(Macs.Spectrum$Feature)

B1 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) 
B2 <- B1 + geom_point(aes(color = Feature), alpha = 0.5)
B3 <- B2 + geom_density2d(color = "grey", linetype =1)
B3

# Set Exp as grouping varibale #

Macs.Spectrum$Exp <- as.factor(Macs.Spectrum$Exp)
A1 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) + facet_grid(~Exp)
A2 <- A1 + geom_point(aes(color = Exp), alpha = 0.5, size = 1.5)
A3 <- A2 + geom_density2d(color = "grey", linetype =1)
A3

A4 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) + facet_grid(~Feature)
A5 <- A4 + geom_point(aes(color = Exp), alpha = 0.5, size = 1.5)
A6 <- A5 + geom_density2d(color = "grey", linetype =1)
A6

A4 <- ggplot(Macs.Spectrum, aes(x= MPI, y = AMDI)) + facet_grid(~Feature)
A5 <- A4 + geom_point(aes(color = Feature), alpha = 0.5, size = 1.5, color = c("salmon", "green", "blue"))
A6 <- A5 + geom_density2d(color = "grey", linetype =1)
A6

# Subset by treatment groups #

NABX <- subset(Macs.Spectrum, Exp=="NABX")
NABX$Feature <- as.factor(NABX$Feature)
ABX <- subset(Macs.Spectrum, Exp=="ABX")
ABX$Feature <- as.factor(ABX$Feature)
C1 <- ggplot(NABX, aes(x= MPI, y = AMDI)) 
C2 <-C1 + geom_point(aes(color = Feature), alpha = 0.5, size = 1.5)
C3 <- C2 + geom_density2d(color = "grey", linetype =1)
C3

D1 <- ggplot(ABX, aes(x= MPI, y = AMDI)) 
D2 <-D1 + geom_point(aes(color = Feature), alpha = 0.5, size = 1.5)
D3 <- D2 + geom_density2d(color = "grey", linetype =1)
D3

C3+D3

# Now by Clusters #

# Subset by treatment groups#
C0 <- subset(Macs.Spectrum, Feature=="0")
C0$Exp <- as.factor(C0$Exp)

C1 <- subset(Macs.Spectrum, Feature=="1")
C1$Feature <- as.factor(C1$Feature)

C2 <- subset(Macs.Spectrum, Feature=="2")
C2$Feature <- as.factor(C2$Feature)

C0.1 <- ggplot(C0, aes(x= MPI, y = AMDI)) + xlim(-10,10) + ylim(-30,30) + geom_point(alpha = 0.5, size = 1.5, color = "salmon") + geom_density2d(color = "grey", linetype =1) 
C1.1 <- ggplot(C1, aes(x= MPI, y = AMDI)) + xlim(-10,10) + ylim(-30,30) + geom_point(alpha = 0.5, size = 1.5, color = "limeggreen") + geom_density2d(color = "grey", linetype =1)
C2.1 <- ggplot(C2, aes(x= MPI, y = AMDI)) + xlim(-10,10) + ylim(-30,30) + geom_point(alpha = 0.5, size = 1.5, color = "steelblue1") + geom_density2d(color = "grey", linetype =1)
C0.1+C1.1+C2.1


CD4TCM <- c("CD3D", "LTB", "TRAC","LDHB")

# Module score #
I <- AddModuleScore(object = I.Macs.integrated, features = Phagocytosis, nbin = 4, ctrl = 2, assay = "RNA", name = "Phagocytosis")

FeaturePlot(I, features = "Phagocytosis1", min.cutoff = "q10", cols = c("gray90", "red"), pt.size = 2)

# Make feature plots # 

FeaturePlot(I.Macs.integrated, features = Macrophage.migration, reduction = "umap", min.cutoff = "q10", cols = c("grey", "red"))
FeaturePlot(I.Macs.integrated, features = Cell.death, reduction = "umap", min.cutoff = "q10", cols = c("gray90", "salmon"))
FeaturePlot(I.Macs.integrated, features = Phagocytosis, reduction = "umap", min.cutoff = "q10", cols = c("gray90", "salmon"))

P1 <- FeaturePlot(I.Macs.integrated, features = c("CSF1R", "MS4A7"), reduction = "umap", min.cutoff = "q1", cols = c("gray90", "blue"), pt.size = 2)
P1 <- FeaturePlot(I.Macs.integrated, features = c("S100A4", "NLRP3"), reduction = "umap", min.cutoff = "q1", cols = c("gray90", "darkgreen"), pt.size = 2)
P1 <- FeaturePlot(I.Macs.integrated, features = c("LCN2", "MSR1"), reduction = "umap", min.cutoff = "q1", cols = c("gray90", "red"), pt.size = 2)
P1 <- DimPlot(I.Macs.integrated, split.by = "treat", pt.size = 2, label = FALSE) + NoLegend()
P1 <- DimPlot(I.Macs.integrated, pt.size = 2, label = FALSE) + NoLegend()
alpha.use <- 3/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1 =
# Plot heatmaps for for marker genes # 

DefaultAssay(I.Macs.integrated) <- "RNA"
I.Macs.scale <-ScaleData(I.Macs.integrated, features = I.Macs.top.markers$gene)
# Remove duplicates for making heatmaps
List <- I.Macs.top.markers[!duplicated(I.Macs.top.markers$gene), ]
List <- column_to_rownames(List, var = "gene")
List <- rownames(List)

# Plot the heatmap
library(viridis)
DoHeatmap(subset(I.Macs.scale, downsample = 100), features = List, size = 3, group.colors = I.Macs.scale$seurat_clusters, angle = 90) + scale_fill_viridis_b() + NoLegend()

# Plot DimPlot
P1 <- DimPlot(I.Macs.integrated, reduction = "umap", pt.size = 1.5, label = T, split.by = "exp")+NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1

# Plots for changes with treatment #
P1 <- DimPlot(I.Macs.integrated, reduction = "umap", pt.size = 1.5, split.by = "exp")+NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1
saveRDS(I.Macs.integrated, file = "/Path_to_save/I.macs.integrated.rds")

# Make plots for cell distribution # 
pt <- table(Idents(I.Macs.integrated), I.Macs.integrated$exp)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

factor(pt$Var1)

# Change the order of plotting,
my.levels <- c("CD4 Naive","CD4 Central Memory","CD8 Central Memory", "CD8 Naive","CD8 Effector Memory", "CD4 Effector Memory")
pt$Var1 <- factor(pt$Var1, levels= my.levels)
color <-c("#B79F00","#F564E3","#F8766D","#00BA38","#619CFF","#00BFC4")

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 25) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_color_manual(color) +
  theme(legend.title = element_blank())

## T cell analysis ##

# Subset T cells for further analysis based on coarse.id # 
T.cells <- subset(lung.integrated, idents = "T cell")

#Begin processing T cell object #
options(future.globals.maxSize = 4000*1024^2)
DefaultAssay(T.cells) <-"RNA"
T.list <- SplitObject(T.cells, split.by = "treat")
T.list <- T.list[c("NABX", "ABX")]

for (i in 1:length(T.list)) {
  T.list[[i]] <- SCTransform(T.list[[i]], verbose = FALSE)
}

T.features <- SelectIntegrationFeatures(object.list = T.list, nfeatures = 2000)
T.list <- PrepSCTIntegration(object.list = T.list, anchor.features = T.features, 
                             verbose = FALSE)
T.anchors <- FindIntegrationAnchors(object.list = T.list, normalization.method = "SCT", 
                                    anchor.features = T.features, verbose = FALSE)
T.integrated <- IntegrateData(anchorset = T.anchors, normalization.method = "SCT", 
                              verbose = FALSE)
T.integrated <- RunPCA(T.integrated, verbose = FALSE)
T.integrated <- RunUMAP(T.integrated, dims = 1:30)
T.integrated <- FindNeighbors(T.integrated , dims = 1:30)
T.integrated  <- FindClusters(T.integrated , resolution = 0.2)

plots <- DimPlot(T.integrated, split.by = "treat")
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))

saveRDS(T.integrated, file = "/Path_to_save/T.integrated.rds")

# Move the file to Seurat4 for mapping # 

# Reference set can be found at https://azimuth.hubmapconsortium.org/; Human - PBMC set was used # 

DefaultAssay(T.integrated) <-"SCT"
anchors <- FindTransferAnchors(
  reference = reference,
  query = T.integrated,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# Load the T.integrated object with predicted id from Seurat V4 reference mapping # 

T.integrated <- readRDS(file = "/Path_to_read/T.integrated.rds")
head(T.integrated[[]])
Idents(T.integrated) <- "predicted.celltype.l2"
p1 <- DimPlot(T.integrated, reduction = "umap", label = TRUE, repel = TRUE)
Idents(T.integrated) <- "seurat_clusters"
p2 <- DimPlot(T.integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2

#find markers for every cluster compared to all remaining cells, report only the positive ones
T.markers <- FindAllMarkers(T.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(T.markers) # 529 genes
write.csv(T.markers, file = "/Path_to_save/T.markers.csv")

#report top 20 markers for each group, by avglogFC
T.top.markers <-T.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% print(n=250)
View(T.top.markers)
write.csv(T.top.markers, file = "/Path_to_save/T.top.markers.csv")

# Rename idents in T cell subset" 
# '0' = "CD4 Naive", '1'= "CD4 Central Memory", '2' ="CD8 Central Memory", '3'="CD8 Naive", '4'="CD8 Effector Memory", '5'="CD4 Effector Memory" 

T.integrated <-RenameIdents(T.integrated, '0' = "CD4 Naive", '1'= "CD4 Central Memory", '2' ="CD8 Central Memory", '3'="CD8 Naive", '4'="CD8 Effector Memory", '5'="CD4 Effector Memory")

# To generate figure 
DimPlot(T.integrated, reduction = "umap", label = TRUE, repel = TRUE, split.by = "exp")

# Specific Markers
NaiveCD4 <- c("CCR7", "SELL", "TCF7", "IL7R", "GIMAP7", "LEF1", "FCMR", "CD3E", "LEF1", "DGKA", "TCF7")
NaiveCD8 <- c("CD8A", "ACKR4", "C4BPA","FCRLA")
Tregs <- c("FOXP3", "TIGIT", "DUSP4", "IKZF2", "IL2RA", "CTLA4", "CD86")
CD8TEM <- c("CCL5", "GZMB", "GZMK","CXCR6", "ITGAX", "CST7", "TIGIT") 
CD8TCM <- c("IFNG", "CCL4", "NCF4","CXCL3")
CD4TCM <- c("CD3D", "LTB", "TRAC","LDHB")
cell <- c("CD3E", "CD4", "CD8A","CD52", "CD69", "CCR7", "IL7R", "EOMES")

# Calculate Module Scores from GO #
DNA.damage <- c("SOX4","RGCC","RBM38","GADD45A","PLK3")
Apoptosis <- c ("TXNIP", "SOX4", "RGCC", "CTLA4", "PLAC8", "DDX3X", "DEDD2", "STK17A", "SQSTM1", "PIK3CD", "CASP3", "PMAIP1", "GADD45A", "TNFRS10A")
Catabolic.process <- c("RPL23A", "RPL19", "ISG15", "PELO", "UBE2S", "PELI1", "DEDD2", "ABHD17B", "HERPUD1", "RPS15", "SQSTM1", "RPSS4X", "CASP3", "MYLIP", "ZC3H12A", "PMAIP1", "CEMIP2")
Cell.Migration <- c("CD36", "THY1", "FN1", "EMP2", "FBLN1", "IGF2", "PODXL", "AGER", "FBLN1", "AFDN", "TNFSF14")
Cell.adhesion <- c("CCN2", "THY1", "COL3A1", "FN1", "EMP2")
TFGB1.signaling <-c("COL1A2", "CLDN5", "ACVRL1", "COL3A1", "ID1")
Tcell.maturation <-c("HSPG2","COLIA1","CCN2","EDNRB","WASF3","CAV2","AFDN","IGF2","PODXL","THY1","NFIB","HOPX","AGER","ECM1","EGFL7","SLFN5","CLDN5","STMN3","RORC","CCR6","JAGN1","SERPINH1","ACVRL1","CSTA","NDRG2","IRX2","COL4A1","FN1","IGFBP5")

# Module score # 
T.integrated <- AddModuleScore(object = T.integrated, features = TFGB1.signaling, nbin = 4, ctrl = 2, assay = "RNA", name = 'Apoptosis')

# Make feature plots # 
FeaturePlot(T.integrated, features = NaiveCD4, reduction = "umap")
FeaturePlot(T.integrated, features = cell, reduction = "umap")
FeaturePlot(T.integrated, features = NaiveCD8, reduction = "umap")
FeaturePlot(T.integrated, features = Tregs, reduction = "umap")
FeaturePlot(T.integrated, features = CD4TCM, reduction = "umap")
FeaturePlot(T.integrated, features = CD8TCM, reduction = "umap")
FeaturePlot(T.integrated, features = Apoptosis, reduction = "umap")

# Plot heatmaps for for marker genes # 

#Change default assay to RNA #
DefaultAssay(T.integrated) <- "RNA"
T.scale <-ScaleData(T.integrated, features = T.top.markers$gene)

# Remove duplicates for making heatmaps #
List <- T.top.markers[!duplicated(T.top.markers$gene), ]
List <- column_to_rownames(List, var = "gene")
List <- rownames(List)

# Plot the heatmap #
library(viridis)
DoHeatmap(subset(T.scale, downsample = 100), features = List, size = 3, group.colors = T.scale$integrated_snn_res.0.19, angle = 90) + scale_fill_viridis_b() + NoLegend()
# Saved as PDF

# Plot DimPlot #

P1 <- DimPlot(T.integrated, reduction = "umap", pt.size = 1.5)+NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1

A2 <-LabelClusters(A1, id = "ident", color = unique(ggplot_build(p1)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)
A1
A2 # (600 X 600)
saveRDS(T.integrated, file = "/Path_to_save/T.integrated.rds")

# Plots for changes with exp # 
p1 <- DimPlot(T.integrated, reduction = "umap", pt.size = 1.5)+NoLegend()
alpha.use <- 2/5
P2$layers[[1]]$mapping$alpha <- alpha.use
A1 <- p1 + scale_alpha_continuous(range = alpha.use, guide = F)

P2 <- DimPlot(T.integrated, reduction = "umap", split.by = "exp", pt.size = 2.5)+NoLegend()
P2$layers[[1]]$mapping$alpha <- alpha.use
A2 <- P2 + scale_alpha_continuous(range = alpha.use, guide = F)
A3 <-LabelClusters(A2, id = "ident", color = unique(ggplot_build(P2)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)

A2 # (600 X 600)

saveRDS(T.integrated, file = "/Path_to_save/T.integrated.rds")

# Make feature plots # 
FeaturePlot(T.integrated, features = NaiveCD4, reduction = "umap", cols = c("salmon", "grey"), repel = TRUE, min.cutoff = "q25", max.cutoff = "q90")
FeaturePlot(T.integrated, features = NaiveCD8, reduction = "umap")
FeaturePlot(T.integrated, features = Tregs, reduction = "umap", cols = c("salmon", "grey"), repel = TRUE, min.cutoff = "q10", max.cutoff = "q50")
FeaturePlot(T.integrated, features = CD4TCM, reduction = "umap")
FeaturePlot(T.integrated, features = CD8TCM, reduction = "umap")
FeaturePlot(T.integrated, features = Apoptosis, reduction = "umap")

# Make plots for cell distribution # 

pt <- table(Idents(T.integrated), T.integrated$exp)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
factor(pt$Var1)

# Change the order of plotting #

my.levels <- c("CD4 Naive","CD4 Central Memory","CD8 Central Memory", "CD8 Naive","CD8 Effector Memory", "CD4 Effector Memory")
pt$Var1 <- factor(pt$Var1, levels= my.levels)
color <-c("#B79F00","#F564E3","#F8766D","#00BA38","#619CFF","#00BFC4")

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 25) +
  geom_col(position = "fill", width = 0.9) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_color_manual(color) +
  theme(legend.title = element_blank())

# Module score for exhausted T cells #
Genes <-rownames(T.integrated)
Module <- read.csv2(file = "path_to_module_score(SupplementalTable12)", sep = ",")
View(Module)
Tcell.exhaustion <-Module[,8]
Tcell.exhaustion  <- Tcell.exhaustion[Tcell.exhaustion  %in% Genes]

DefaultAssay(T.integrated) <-"RNA"
T.integrated <- AddModuleScore(T.integrated, features = list(Tcell.exhaustion), ctrl.size = 10,name = 'Tcell.exhaustion', nbin = 20)
T.integrated[['module']] <- CreateAssayObject(data = t(x = FetchData(object = T.integrated, vars = 'Tcell.exhaustion1')))

# For viridis colors #
VlnPlot(object = T.integrated, features = 'Tcell.exhaustion1', assay = 'module', slot = 'data', pt.size = 0.25) + NoLegend()
VlnPlot(object = T.integrated, features = 'PLK3', assay = 'integrated', slot = 'data', pt.size = 0.25, group.by = "exp") + NoLegend()

## NK Cell Analysis ##
NK.cells <- subset(lung.integrated, idents = "NK cell") 
NK.cells # 1104369 features across 1030 samples within 6 assays 
DimPlot(NK.cells, reduction = "umap", split.by = "treat", label = TRUE)

options(future.globals.maxSize = 4000*1024^2)
DefaultAssay(NK.cells) <-"RNA"
NK.list <- SplitObject(NK.cells, split.by = "treat")
NK.list <- NK.list[c("NABX", "ABX")]

for (i in 1:length(NK.list)) {
  NK.list[[i]] <- SCTransform(NK.list[[i]], verbose = FALSE)
}

NK.features <- SelectIntegrationFeatures(object.list = NK.list, nfeatures = 2000)
NK.list <- PrepSCTIntegration(object.list = NK.list, anchor.features = NK.features, 
                              verbose = FALSE)
NK.anchors <- FindIntegrationAnchors(object.list = NK.list, normalization.method = "SCT", 
                                     anchor.features = NK.features, verbose = FALSE)
NK.integrated <- IntegrateData(anchorset = NK.anchors, normalization.method = "SCT", 
                               verbose = FALSE)
NK.integrated <- RunPCA(NK.integrated, verbose = FALSE)
NK.integrated <- RunUMAP(NK.integrated, dims = 1:30)
NK.integrated <- FindNeighbors(NK.integrated , dims = 1:30)
NK.integrated  <- FindClusters(NK.integrated , resolution = 1)

plots <- DimPlot(NK.integrated, split.by = "treat")
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))
saveRDS(NK.integrated, file = "/Path_to_save/NK.integrated.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones #
NK.markers <- FindAllMarkers(NK.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(NK.markers) # 529 genes
write.csv(NK.markers, file = "/Path_to_save/NK.markers.csv")

#report top 20 markers for each group, by avglogFC #
NK.top.markers <-NK.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% print(n=250)
View(NK.top.markers)
write.csv(NK.top.markers, file = "/Path_to_save/NK.top.markers.csv")

# To generate figure #
P1 <- DimPlot(NK.integrated, reduction = "umap", pt.size = 1.5, label = F, split.by = "exp")+NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1

# Specific Markers #
CD56.NK <- c("CCR7", "SELL", "TCF7", "IL7R", "GIMAP7", "LEF1", "FCMR", "CD3E", "LEF1", "DGKA", "TCF7", "CD56")
Proliferating.NK <- c("PCNA", "MCN2")
Inflammatory.NK <- c("CD38", "CD52", "NK2G7", "NFKB1", "IRF7", "IRF8", "CD69", "PRF1")

# Calculate Module Scores from GO #
DNA.damage <- c("SOX4","RGCC","RBM38","GADD45A","PLK3")
Apoptosis <- c ("TXNIP", "SOX4", "RGCC", "CTLA4", "PLAC8", "DDX3X", "DEDD2", "STK17A", "SQSTM1", "PIK3CD", "CASP3", "PMAIP1", "GADD45A", "TNFRS10A")
Catabolic.process <- c("RPL23A", "RPL19", "ISG15", "PELO", "UBE2S", "PELI1", "DEDD2", "ABHD17B", "HERPUD1", "RPS15", "SQSTM1", "RPSS4X", "CASP3", "MYLIP", "ZC3H12A", "PMAIP1", "CEMIP2")
Cell.Migration <- c("CD36", "THY1", "FN1", "EMP2", "FBLN1", "IGF2", "PODXL", "AGER", "FBLN1", "AFDN", "TNFSF14")
Cell.adhesion <- c("CCN2", "THY1", "COL3A1", "FN1", "EMP2")
TFGB1.signaling <-c("COL1A2", "CLDN5", "ACVRL1", "COL3A1", "ID1")
NKcell.maturation <-c("HSPG2","COLIA1","CCN2","EDNRB","WASF3","CAV2","AFDN","IGF2","PODXL","THY1","NFIB","HOPX","AGER","ECM1","EGFL7","SLFN5","CLDN5","STMN3","RORC","CCR6","JAGN1","SERPINH1","ACVRL1","CSTA","NDRG2","IRX2","COL4A1","FN1","IGFBP5")

# Module score # 
NK.integrated <- AddModuleScore(object = NK.integrated, features = TFGB1.signaling, nbin = 4, ctrl = 2, assay = "RNA", name = 'Apoptosis')

# Make feature plots # 
FeaturePlot(NK.integrated, features = CD56.NK, reduction = "umap")
FeaturePlot(NK.integrated, features = Inflammatory.NK, reduction = "umap")
FeaturePlot(NK.integrated, features = Proliferating.NK, reduction = "umap")
FeaturePlot(NK.integrated, features = DNA.damage, reduction = "umap")
FeaturePlot(NK.integrated, features = Apoptosis, reduction = "umap")

# Plot heatmaps for for marker genes # 
DefaultAssay(NK.integrated) <- "RNA"
NK.scale <-ScaleData(NK.integrated, features = NK.top.markers$gene)
# Remove duplicates for making heatmaps #
List <- NK.top.markers[!duplicated(NK.top.markers$gene), ]
List <- column_to_rownames(List, var = "gene")
List <- rownames(List)

# Plot the heatmap #
library(viridis)
DoHeatmap(subset(NK.scale, downsample = 100), features = List, size = 3, group.colors = NK.scale$integrated_snn_res.1, angle = 90) + scale_fill_viridis_b() + NoLegend()

# Plot UMAP
P1 <- DimPlot(NK.integrated, reduction = "umap", pt.size = 1.5)+NoLegend()
alpha.use <- 2/5
P1$layers[[1]]$mapping$alpha <- alpha.use
A1 <- P1 + scale_alpha_continuous(range = alpha.use, guide = F)
A1

A2 <-LabelClusters(A1, id = "ident", color = unique(ggplot_build(p1)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)
A1
A2 # (425 X 425)

## CellChat analysis ##
# Comparison analysis of NABX and ABX exp conditions using CellChat # 
library(CellChat)
library(patchwork)

# load cell chat objects if necessary #
cellchat.NABX <- readRDS(file = "/Path_to_read/Cell.Chat.NABX.rds")
cellchat.ABX <- readRDS(file = "/Path_to_read/Cell.Chat.ABX.rds")
object.list <- list(NABX = cellchat.NABX, ABX = cellchat.ABX)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c("blue", "red"))
col <- c("#F8766D", "#DE8C00", "#B79F00", "#7CAE00", "#00BA38", "#00C08B", "#C77CFF", "#FF64B0")
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, color.use = col, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.use = col)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.use = col)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Incoming patterns #
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu", color.use = col)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu", color.use = col)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Overall patterns #
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 4, height = 12, color.heatmap = "OrRd",color.use = col)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 4, height = 12, color.heatmap = "OrRd", color.use = col)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm")) # (Fig. 7C)

# Identify the upregulated and down-regulated signaling ligand-receptor pairs (from macrophages to other cells) #
LR <- c("CCL", "CXCL", "MIF", "IL1", "CSF", "NOTCH", "SEMA4", "SELPLG", "GAS", "PROS", "PD-L1", "PDL2", "CD45", "CD6", "COMPLEMENT", "CD86", "CLEC", "THBS", "THY1", "CD80", "VEGF", "IL16", "NECTIN", "TNF")
netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1:2,4:8),  signaling = LR, comparison = c(1, 2), angle.x = 90, color.text = c("blue", "red"), thresh = 0.01)

# Identify the upregulated and down-regulated signaling ligand-receptor pairs (from all to neutrophils) #
netVisual_bubble(cellchat, sources.use = c(2:8), targets.use = c(1),  signaling = LR, comparison = c(1, 2), angle.x = 90, color.text = c("blue", "red"), thresh = 0.01)

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, color.use = col, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("THBS") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("SELPLG") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("CSF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("MIF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("NOTCH") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("CD45") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}


pathways.show <- c("SEMA4") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

pathways.show <- c("PDL2")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 14, vertex.size.max = 30, signaling.name = paste(pathways.show, names(object.list)[i]), vertex.label.cex = 2)
}

# DeSeq2 analysis of pseudobulk data in Fig 6 #

library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

# Bring in Seurat object if not done already #
seurat <- readRDS(file = "/Path_to_data/lung.integrated.rds")
Idents(seurat) <-"treat"
seurat <- subset(seurat, idents = c("NABX", "ABX", "FT"))
Idents(seurat) <-"seurat_clusters"

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA@counts
metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@meta.data$seurat_clusters)
metadata$sample_id <- factor(seurat@meta.data$id)
metadata$group_id <- factor(seurat@meta.data$exp)

# Create single cell experiment object #
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

# Identify groups for aggregation of counts #
groups <- colData(sce)[, c("cluster_id", "id", "group_id")]
head(colData(sce))

# Named vector of cluster names #
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters #
nk <- length(kids)
nk

# Named vector of sample names #
sids <- purrr::set_names(levels(sce$sample_id))
sids

# Total number of samples #
ns <- length(sids)
ns

# Named vector of group names #
tids <- purrr::set_names(levels(sce$group_id))
tids

# Total number of samples #
nt <- length(tids)
nt

# To perform sample-level differential expression analysis,
# we need to generate sample-level metadata. To do this,
# we will reorder samples in the single-cell metadata to match the order of the
# factor levels of the sample ID, then extract only the sample-level information
# from the first cell corresponding to that sample.

## Determine the number of cells per sample ##
table(sce$sample_id)
## Turn named vector into a numeric vector of number of cells per sample ##
n_cells <- as.numeric(table(sce$sample_id))
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector ##
m <- match(sids, sce$sample_id)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample ##
ei <- data.frame(colData(sce)[m, ],
                 n_cells, row.names = NULL) %>%
  select(-"cluster_id")
ei
# QC: remove count outliers and low count genes using functions from the scater package #
dim(sce)
# Calculate quality control (QC) metrics #
sce <- calculateQCMetrics(sce)
# Get cells w/ few/many detected genes #
sce$is_outlier <- isOutlier(
  metric = sce$total,
  nmads = 2, type = "both", log = TRUE)
# Remove outlier cells #
sce <- sce[, !sce$is_outlier]
dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts ##
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Count aggregation to sample level #
# Aggregate the counts per sample_id and cluster_id #
# Subset metadata to only include the cluster and sample IDs to aggregate across #
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Aggregate across cluster-sample groups #
pb <- aggregate.Matrix(t(counts(sce)),
                       groupings = groups, fun = "sum")
# Not every cluster is present in all samples; create a vector that represents how to split samples #
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_",
                                    n = 2),
                 `[`, 1)
# Now we can turn the matrix into a list that is split into count matrices 
# for each cluster, then transform each data frame so that rows are genes and columns are the samples #
# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs #
pb <- split.data.frame(pb,
                       factor(splitf)) %>%
  lapply(function(u)
    set_colnames(t(u),
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
# Print out the table of cells in each cluster-sample group #
options(width = 100)
table(sce$cluster_id, sce$sample_id)

## Differential gene expression with DESeq2 ##

# Get sample names for each of the cell type clusters #
# prep. data.frame for plotting #
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>% unlist()
# Then we can get the cluster IDs corresponding to each of the samples in the vector. #
# Get cluster IDs for each of the samples #
samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
  rep(names(pb)[x],
      each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Finally, let’s create a data frame with the cluster IDs and the
# corresponding sample IDs. We will merge together the condition information. #
# Create a data frame with the sample IDs, cluster IDs and condition #
gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")])
metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id)
metadata
# Subsetting dataset to cluster(s) of interest #
# Now that we have the sample-level metadata, we can run the differential expression analysis with DESeq2. #
# Oftentimes, we would like to perform the analysis on multiple different clusters, so we can set up the workflow to run easily on any of our clusters. #

# Generate vector of cluster IDs #
clusters <- levels(metadata$cluster_id)
clusters
# [1] "0"  "1"  "10" "11" "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9" #

clusters[1] # whatever cell type you want to extract - Neutrophil (0) or T cell (1) or Alveolar Macrophage
# Subset the metadata to only the Neutrophil cells #
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs #
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)
# Subset the counts to only the Neutrophils or T cells #
counts <- pb[[clusters[1]]] # 1 is for Neutrophil
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2 #
all(rownames(cluster_metadata) == colnames(cluster_counts))

# Create DEQ2 object #
dds <- DESeqDataSetFromMatrix(cluster_counts,
                              colData = cluster_metadata,
                              design = ~ group_id)
## Transform counts for data visualization ##
rld <- rlog(dds, blind=TRUE)
# Plot PCA #
DESeq2::plotPCA(rld, intgroup = "group_id")
# Extract the rlog matrix from the object and compute pairwise correlation values #
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
# Plot heatmap #
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
# Run DESeq2 differential expression analysis #
dds <- DESeq(dds)
plotDispEsts(dds)
# Output results of Wald test for contrast for NBAX vs ABX #
contrast <- c("group_id", "NABX", "ABX")
# resultsNames(dds) #
res <- results(dds,
               contrast = contrast,
               alpha = 0.05)
res <- lfcShrink(dds,
                 contrast =  contrast,
                 res=res)
# First let’s generate the results table for all of our results: 
# Turn the results object into a tibble for use with tidyverse functions #
  res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output #
res_tbl
# Write all results to file #
write.csv2(res_tbl, file = "/Path_to_save/Neutrophil_NABXvsABX_all_genes.csv", quote = FALSE, row.names = FALSE)
# Next, we can filter our table for only the significant genes using a p-adjusted threshold of 0.05 #

# Subset the significant results #
pvalue_cutoff <- 0.05 # set p value cutoff.
sig_res <- dplyr::filter(res_tbl, pvalue < pvalue_cutoff) %>%
dplyr::arrange(pvalue)
# Check significant genes output #
sig_res
# Write significant results to file #
write.table(sig_res, file = "/Path_to_save/Neutrophil_NABXvsABX_sig_genes.csv", sep = ",", quote = FALSE, row.names = FALSE)
# Scatter Plor of nromalized expression of top 20 most signficiant genes #
## ggplot of top genes ##

normalized_counts <- counts(dds,
normalized = TRUE)
## Order results by padj values ##
top20_sig_genes <- sig_res %>%
dplyr::arrange(pvalue) %>%
dplyr::pull(gene) %>%
head(n=20)

top20_sig_norm <- data.frame(normalized_counts) %>%
rownames_to_column(var = "gene") %>%
dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))

## plot using ggplot2 ##

ggplot(gathered_top20_sig) +
geom_point(aes(x = gene,
y = normalized_counts,
color = group_id),
position=position_jitter(w=0.1,h=0)) +
scale_y_log10() +
xlab("Genes") +
ylab("log10 Normalized Counts") +
ggtitle("Top 20 Significant DE Genes") +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
theme(plot.title = element_text(hjust = 0.5))

# Heatmap of all significant genes #
# Extract normalized counts for only the significant genes #
sig_norm <- data.frame(normalized_counts) %>%
rownames_to_column(var = "gene") %>%
dplyr::filter(gene %in% sig_res$gene)

# Set a color palette #
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation #
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
color = viridis(100, alpha = 0.5, option = "D"),
cluster_rows = T,
show_rownames = T,
annotation = cluster_metadata[, c("group_id", "cluster_id")],
border_color = NA,
fontsize = 10,
scale = "row",
fontsize_row = 10,
height = 20,
row_split = 3, column_split = 2,
row_gap = unit(1, "mm"))

# Volcano plot of results #
res_table_thres <- res_tbl %>%
mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot ##
ggplot(res_table_thres) +
geom_point(aes(x = log2FoldChange, y = -log10(pvalue), colour = threshold)) +
ggtitle("Volcano plot of stimulated B cells relative to control") +
xlab("log2 fold change") +
ylab("-log10 adjusted p-value") +
scale_y_continuous(limits = c(0,5)) +
theme(legend.position = "none",
plot.title = element_text(size = rel(1.5), hjust = 0.5),
axis.title = element_text(size = rel(1.25)))

# Cytokine and cell frequency data analysis # 
# Load cytokine data 
Cytokine <- read.csv(file ="/path_to_saves/Postnatal_Adaptation_Siganture_Cytokines.csv", sep = ",", header = T)
Metadata <- read.csv(file ="/path_to_save/Postnatal_Adaptation_Siganture_Metadata.csv", sep = ",", header = T)
Cells <- read.csv(file ="/path_to_save/Postnatal_Adaptation_Siganture_Cytokines.csv", sep = ",", header = T)
Cytokine <- Cytokine %>% remove_rownames %>% column_to_rownames(var="ID")
Cells <- Cells %>% remove_rownames %>% column_to_rownames(var="ID")
Metadata <- Metadata %>% remove_rownames %>% column_to_rownames(var="ID")
str(Cytokine)
summary(Cytokine)
summary(Cells)

# Scaling/standardizing the data using the R function scale
Scale.Cytokine <- scale(Cytokine)
Scale.Cells <- scale(Cells)
heat_colors <- brewer.pal(9, "YlOrRd")
Heatmap(Scale.Cytokine)


# Run pheatmap using the metadata data frame for the annotation
pheatmap(Scale.Cytokine, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         kmeans_k = 4,
         annotation = Metadata[,c("Group","Age")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)     

library(pheatmap)

pheatmap(Scale.Cells, 
         col = viridis(100, alpha = 0.5, option = "D"), 
         cluster_rows = T, 
         show_rownames = T,
         kmeans_k = 4,
         annotation = Metadata[,c("Group","Age")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 

# BAL Cytokines
Cytokine <- read.csv(file ="/path_to_save/BAL_Cytokines.csv", sep = ",", header = T)
Metadata <- read.csv(file ="/path_to_save/BAL_Metadata.csv", sep = ",", header = T)
Cytokine <- Cytokine %>% remove_rownames %>% column_to_rownames(var="ID")
Metadata <- Metadata %>% remove_rownames %>% column_to_rownames(var="ID")
str(Cytokine)
summary(Cytokine)

# Scaling/standardizing the data using the R function scale
Scale.Cytokine <- scale(Cytokine)

# Run pheatmap using the metadata data frame for the annotation
pheatmap(Scale.Cytokine, 
         col = viridis(100, alpha = 0.5, option = "D"), 
         cluster_rows = T, 
         show_rownames = T,
         annotation = Metadata[,c("Group")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row",
         row_split = 3, column_split = 2,
         row_gap = unit(1, "mm")
)     

# Correlation between different cytokines and immune cells in BAL Fluid
Use_Data <- read.csv(file ="/path_to_save/BAL_Cytokines_Cells_Correlation.csv", sep = ",", header = T)

# Neutrophil correlation
P1<- ggscatter(Use_Data, x = "IL8", y = "Granulocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL8(pg/ml)", ylab = "Neutrophils (% of CD45+ cells)")
P2 <- ggscatter(Use_Data, x = "CXCL2", y = "Granulocytes", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL2(pg/ml)", ylab = "Neutrophils (% of CD45+ cells)")
P3 <- ggscatter(Use_Data, x = "CXCL8", y = "Granulocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL8(pg/ml)", ylab = "Neutrophils (% of CD45+ cells)")
P4 <- ggscatter(Use_Data, x = "CXCL10", y = "Granulocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL10(pg/ml)", ylab = "Neutrophils (% of CD45+ cells)")

# Macrophage correlation
A1<- ggscatter(Use_Data, x = "CCL3", y = "Macrophages", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CCL3(pg/ml)", ylab = "Macrophages (% of CD45+ cells)")
A2 <- ggscatter(Use_Data, x = "CCL4", y = "Macrophages", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CCL4(pg/ml)", ylab = "Macrophages (% of CD45+ cells)")
A3 <- ggscatter(Use_Data, x = "CCL11", y = "Macrophages", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CCL11(pg/ml)", ylab = "Macrophages (% of CD45+ cells)")
A4 <- ggscatter(Use_Data, x = "GMCSF", y = "Macrophages", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "GMCSF(pg/ml)", ylab = "Macrophages (% of CD45+ cells)")
A5 <- ggscatter(Use_Data, x = "CXCL10", y = "Granulocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL10(pg/ml)", ylab = "Macrophages (% of CD45+ cells)")

# Lymphocyte correlation 
B1<- ggscatter(Use_Data, x = "IL1B", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL1B(pg/ml)", ylab = "Lymphocytes (% of CD45+ cells)")
B2 <- ggscatter(Use_Data, x = "IL6", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL6(pg/ml)", ylab = "Lymphocytes (% of CD45+ cells)")
B3 <- ggscatter(Use_Data, x = "IL8", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL8(pg/ml)", ylab = "Lymphocytes (% of CD45+ cells)")
B4 <- ggscatter(Use_Data, x = "IL10", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL10(pg/ml)", ylab = "Lymphocytes (% of CD45+ cells)")
B5 <- ggscatter(Use_Data, x = "IL12", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL12(pg/ml)", ylab = "Lymphocytes (% of CD45+ cells)")

B1+B2+B3+B5

# Inflammatory cytokine correlation 
L1<- ggscatter(Use_Data, x = "IL1B", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL1B(pg/ml)", ylab = "Peak PEWS")
L2 <- ggscatter(Use_Data, x = "IL6", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL6(pg/ml)", ylab = "Peak PEWS")
L3 <- ggscatter(Use_Data, x = "IL8", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL8(pg/ml)", ylab = "Peak PEWS")
L4 <- ggscatter(Use_Data, x = "TNFA", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "TNFα(pg/ml)", ylab = "Peak PEWS")
L5 <- ggscatter(Use_Data, x = "TGFA", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "TGFβ(pg/ml)", ylab = "Peak PEWS")

L7 <- ggscatter(Use_Data, x = "CXCL2", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL2(pg/ml)", ylab = "Peak PEWS")
L8 <- ggscatter(Use_Data, x = "CXCL10", y = "Lymphocytes", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL10(pg/ml)", ylab = "Peak PEWS")

L9 <- ggscatter(Use_Data, x = "PDGFBB", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "PDGFBB(pg/ml)", ylab = "Peak PEWS")
L10 <- ggscatter(Use_Data, x = "VEGF", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "VEGF(pg/ml)", ylab = "Peak PEWS")
L11 <- ggscatter(Use_Data, x = "IL10", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "red", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL10(pg/ml)", ylab = "Peak PEWS")

L12 <- ggscatter(Use_Data, x = "CXCL8", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL8(pg/ml)", ylab = "Peak PEWS")
L13 <- ggscatter(Use_Data, x = "CXCL10", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CXCL10(pg/ml)", ylab = "Peak PEWS")
L14 <- ggscatter(Use_Data, x = "TGFβ", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "TGFβ(pg/ml)", ylab = "Peak PEWS")


L15 <- ggscatter(Use_Data, x = "IL6", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL6(pg/ml)", ylab = "Peak PEWS")
L16 <- ggscatter(Use_Data, x = "IL8", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "IL8(pg/ml)", ylab = "Peak PEWS")
L17 <- ggscatter(Use_Data, x = "TNFα", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "TNFα(pg/ml)", ylab = "Peak PEWS")


L18 <- ggscatter(Use_Data, x = "M1.AM", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "M1-activated AM (% of MHCII+,CD11C+,CD11B+ cells", ylab = "Peak PEWS")
L19 <- ggscatter(Use_Data, x = "Dysfunctional.T.cells", y = "Peak.PEWS", size = 2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
                 conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "CD279+CD38+ CD4+ Tcells (% of CD4+ T cells)", ylab = "Peak PEWS")


## END ANALYSIS ##





