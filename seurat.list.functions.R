# https://satijalab.org/seurat/essential_commands.html



library(hdf5r)
library(SingleCellExperiment)
library(biomaRt)
library(Seurat)
library(dplyr)
library(singleCellTK)
library(Seurat)
library(scater)
library(singleCellTK)
library(SCHNAPPs)
library(iSEE)

load("pbmc5k.seurat.RData")

pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- RunPCA(pbmc,npcs = 50)
pbmc <- RunTSNE(pbmc)
pbmc <- RunICA(pbmc)
pbmc <- RunLSI(pbmc)
# RunUMAP(pbmc, dims = 3:50)


#================================================

DimPlot(object = pbmc, reduction = 'pca')
DimPlot(object = pbmc, reduction = 'tsne')
DimPlot(object = pbmc, reduction = 'ica')
DimPlot(object = pbmc, reduction = 'lsi')


# projected cells, i.e. Eigenvector * input vector
Embeddings(pbmc)[1:5,1:5]
pbmc[['pca']][[1:5,1:5]]

# Eigenvalues
Loadings(pbmc)
pbmc[['pca']][1:5,1:5]

dIca <- DimHeatmap(pbmc, dims = 2, nfeatures = 10, reduction = 'ica', fast = F)
dLsi <- DimHeatmap(pbmc, dims = 2, nfeatures = 10, reduction = 'lsi', fast = F)
dPca <- DimHeatmap(pbmc, dims = 2, nfeatures = 10, reduction = 'pca', fast = F)

CombinePlots(plots = list(dIca, dLsi, dPca))

#=============== Plotting =================================
# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "MS4A1")

# Scatter plot across single cells, replaces GenePlot
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "PC_1")
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "CD3D")


FeatureScatter(object = pbmc, feature1 = 'CD9', feature2 = 'CD3E', slot = 'data')
FeatureScatter(object = pbmc, feature1 = 'CD9', feature2 = 'CD3E', slot = 'counts')
FeatureScatter(object = pbmc, feature1 = 'CD9', feature2 = 'CD3E', slot = 'logcounts')

# Scatter plot across individual features, repleaces CellPlot
CellScatter(object = pbmc, cell1 = colnames(pbmc)[1], cell2 =colnames(pbmc)[2])

VariableFeaturePlot(object = pbmc)


# Violin and Ridge plots
VlnPlot(object = pbmc, features = c("LYZ", "CCL5", "IL32"))
RidgePlot(object = pbmc, feature = c("LYZ", "CCL5", "IL32"))

# Heatmaps
DoHeatmap(object = pbmc, features = c("LYZ", "CCL5", "IL32"))
DimHeatmap(object = pbmc, reduction = "pca", cells = 200)


# New things to try!  Note that plotting functions now return ggplot2 objects, so you can add themes, titles, and options
# onto them
FeaturePlot(object = pbmc, features = c("MS4A1", "CD79A"), blend = TRUE)
DimPlot(object = pbmc) + DarkTheme()
DimPlot(object = pbmc) + labs(title = "2,700 PBMCs clustered using Seurat and viewed\non a two-dimensional tSNE")

# HoverLocator replaces the former `do.hover` argument It can also show extra data throught the `information` argument,
# designed to work smoothly with FetchData
plot <- DimPlot(object = pbmc) + NoLegend()
HoverLocator(plot = plot, information = FetchData(object = pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))

# FeatureLocator replaces the former `do.identify`
select.cells <- CellSelector(plot = plot)

# Label points on a ggplot object
LabelPoints(plot = plot, points = TopCells(object = pbmc[["pca"]]), repel = TRUE)

#================================================

# AddModuleScore
cd_features <- list(c(
  'CD79B',
  'CD79A',
  'CD19',
  'CD180',
  'CD200',
  'CD3D',
  'CD2',
  'CD3E',
  'CD7',
  'CD8A',
  'CD14',
  'CD1C',
  'CD68',
  'CD9',
  'CD247'
))
pbmc <- AddModuleScore(
  object = pbmc,
  features = cd_features,
  ctrl = 5,
  name = 'CD_Features'
)
head(x = pbmc[])
# RunALRA 
# Runs ALRA, a method for imputation of dropped out values in scRNA-seq data.
# ALRAChooseKPlot

as.SingleCellExperiment(x, ...)

AverageExpression	Averaged feature expression by identity class
Returns expression for an 'average' single cell in each identity class

#================================================

> cc.genes
$s.genes
[1] "MCM5"     "PCNA"     "TYMS"     "FEN1"     "MCM2"     "MCM4"     "RRM1"     "UNG"      "GINS2"    "MCM6"     "CDCA7"   
[12] "DTL"      "PRIM1"    "UHRF1"    "MLF1IP"   "HELLS"    "RFC2"     "RPA2"     "NASP"     "RAD51AP1" "GMNN"     "WDR76"   
[23] "SLBP"     "CCNE2"    "UBR7"     "POLD3"    "MSH2"     "ATAD2"    "RAD51"    "RRM2"     "CDC45"    "CDC6"     "EXO1"    
[34] "TIPIN"    "DSCC1"    "BLM"      "CASP8AP2" "USP1"     "CLSPN"    "POLA1"    "CHAF1B"   "BRIP1"    "E2F8"    

$g2m.genes
[1] "HMGB2"   "CDK1"    "NUSAP1"  "UBE2C"   "BIRC5"   "TPX2"    "TOP2A"   "NDC80"   "CKS2"    "NUF2"    "CKS1B"   "MKI67"  
[13] "TMPO"    "CENPF"   "TACC3"   "FAM64A"  "SMC4"    "CCNB2"   "CKAP2L"  "CKAP2"   "AURKB"   "BUB1"    "KIF11"   "ANP32E" 
[25] "TUBB4B"  "GTSE1"   "KIF20B"  "HJURP"   "CDCA3"   "HN1"     "CDC20"   "TTK"     "CDC25C"  "KIF2C"   "RANGAP1" "NCAPD2" 
[37] "DLGAP5"  "CDCA2"   "CDCA8"   "ECT2"    "KIF23"   "HMMR"    "AURKA"   "PSRC1"   "ANLN"    "LBR"     "CKAP5"   "CENPE"  
[49] "CTCF"    "NEK2"    "G2E3"    "GAS2L3"  "CBX5"    "CENPA"  
pbmc_small <- CellCycleScoring(
  object = pbmc_small,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)


#================================================
? How to use log counts?
CellScatter(object = pbmc, cell1 = "AGTTAGCTCACTGGGC-1", cell2 = "AGTTAGCTCTTCCCAG-1")

# get cell names
Cells(pbmc)


# cell selection
plot <- DimPlot(object = pbmc_small)
# Follow instructions in the terminal to select points
cells.located <- CellSelector(plot = plot)
cells.located
# Automatically set the identity class of selected cells and return a new Seurat object
pbmc_small <- CellSelector(plot = plot, object = pbmc_small, ident = 'SelectedCells')






