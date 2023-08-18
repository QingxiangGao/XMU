require(Seurat)
require(dplyr)
require(harmony)
require(cowplot)
##load data 
pbsclp.data <- read.table("pbsclp_pbmc.tsv",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
plasclp.data <- read.table("plasclp_pbmc.tsv",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
pbsclp <- CreateSeuratObject(pbsclp.data,project = "pbsclp",assay = "RNA",min.cells = 3,min.features = 200)
plasclp <- CreateSeuratObject(plasclp.data,project = "plasclp",assay = "RNA",min.cells = 3,min.features = 200)
## merge SeuratObject
CLP <- merge(pbsclp,plasclp,add.cell.ids = c("PBSCLP","PlasCLP"))
VlnPlot(CLP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
plot1 <- FeatureScatter(CLP, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
plot2 <- FeatureScatter(CLP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")
plot1 + plot2
CLP <- subset(CLP, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 5)
VlnPlot(CLP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##Normalization & FindVariableFeature
CLP <- NormalizeData(CLP, normalization.method = "LogNormalize", scale.factor = 10000)
CLP <- FindVariableFeatures(CLP, selection.method = "mvp")
top10 <- head(VariableFeatures(CLP), 10)
top10
plot1 <- VariableFeaturePlot(CLP)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(CLP)
CLP <- ScaleData(CLP, features = all.genes)
CLP <- RunPCA(CLP, features = VariableFeatures(object = CLP),nfeatures.print = 10)
print(CLP[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(CLP, dims = 1:2, reduction = "pca")
DimPlot(CLP, reduction = "pca")
##scale Data with cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(stringr)
s.genes <- tolower(s.genes)
s.genes <- stringr::str_to_title(s.genes)
g2m.genes <- tolower(g2m.genes)
g2m.genes <- str_to_title(g2m.genes)
CLP <- CellCycleScoring(CLP,s.features = s.genes,g2m.features = g2m.genes,set.ident = T)
RidgePlot(CLP, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2,group.by = "orig.ident")
CLP <- ScaleData(CLP, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CLP))
CLP <- RunPCA(CLP, features = VariableFeatures(CLP), nfeatures.print = 10)
CLP <- RunPCA(CLP, features = c(s.genes, g2m.genes))
DimPlot(CLP)
##harmony reduce batch effector
options(repr.plot.height = 2.5, repr.plot.width = 6)
CLP <- CLP %>% 
  RunHarmony(group.by.vars = "orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(CLP, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = CLP, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = CLP, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
##Findcluster 
CLP <- CLP %>% 
  RunUMAP(reduction = "harmony", dims = 1:3) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:3) %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()
##Draw plot
DimPlot(CLP,reduction = "umap",label = T)
DimPlot(CLP,reduction = "umap",group.by = "orig.ident",split.by = "orig.ident")
FeaturePlot(CLP,reduction = "umap",features = c("Cd3g","Cd19",.....),slot = "count")
VlnPlot(CLP,features = c("Cd3g","Cd19",.....),slot = "count",split.by = "orig.ident",group.by = "seurat_clusters")
