#!/bin/env Rscript

#########################################################################################################
# PART A) prepare input data for the seurat part of the lab
#########################################################################################################
require(anndata)
require(tidyverse)
infile = 'scder.h5ad'
data = read_h5ad(infile)

mat.outfile = 'scCountMat.csv.gz'
meta.outfile = 'scCountMetadata.csv.gz'
dat = as.matrix(data$X) %>% as_tibble()

dat  %>% 
	write_csv(file = mat.outfile)

meta.df = data$obs
colnames(meta.df)[1] = 'cellName'
write_csv(meta.df,
	file = meta.outfile)

#########################################################################################################
# PART B) Do data IO and processing until clustering with seurat
#########################################################################################################

library(Seurat)
library(tidyverse)

matfile  <- 'scCountMat.csv.gz'
metafile <- 'scCountMetadata.csv.gz'

counts   <- read.csv(matfile, header = TRUE, check.names=FALSE)
meta     <- read.csv(metafile)
rownames(meta) 		<- meta$cellName
rownames(counts) 	<- meta$cellName

## create Seurat Object
scde <- CreateSeuratObject(
	t(counts),
	project = "SeuratProject",
	assay = "RNA",
	meta.data = meta
)
# We will not use the provided cell type identity for now, i.e. setting 
# cell identities in the seurat object to an empty string.
Idents(scde) = ''

# Let's check if the meta data are correctly assigned.
# The cluster.id should provide an annotation of cells belonging to the 
# different subtypes TCM, TEM, TEMRA, and TN.
table(scde[['cluster.id']])

#########################################################################################################
# QC
#########################################################################################################
# We will filter out cells with less than 10000 read counts in total.
# First, we calculate the total number of counts and store it as metadata
scde$cTot <- apply(GetAssayData(object = scde, slot = "counts"), 2, sum)
# Now we can apply the filter criterion
scde <- subset(scde, 
	subset = cTot > 10000)

# another common filter criterion is to reomve 
# cells with very low or high number of detected genes, which 
# likely represent empty droplets or e.g. duplets.
# Here, we want to filter for empty droplets by removing cells with less 
# than 200 detected genes.
scde <- subset(scde, subset = nFeature_RNA > 200)


# The next step in quality control is filtering cells by fraction of
# mitochondrial reads. The first step is to calculate the fractions.
scde[["percent.mt"]] <- PercentageFeatureSet(scde, pattern = "^MT-")
# Secondly, the distribution can be visualized as violin plot, here together with the 
# distribution of the number of detected genes per cell, and the number of unique molecules per cell.
VlnPlot(scde, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Let's filter out all cells with more than 5% mitochondrial reads.
scde <- subset(scde, subset = percent.mt < 5)

#########################################################################################################
# Preprocessing
#########################################################################################################
# After filtering out low quality cells from our dataset, we can start with the 
# actual processing. The first step is to normalize the data. This includes 
# dividing the counts by the total number of counts per cell, multiplying by an
# arbitrary scaling factor and then log-transforming the matrix.
# The normalized data are stored in the slot sdce[["RNA"]]@data
scde <- NormalizeData(scde)

# Next, we will define the set of features that exhibit the highest expression 
# variability across the cells. By default the set of 2000 most variable features is 
# selected, but this parameter can be adjusted according to the dataset.
scde <- FindVariableFeatures(scde, selection.method = "vst", nfeatures = 2000)

# A plot showing the variable features adding labels to the
# 10 most highly variable genes
top10 <- head(VariableFeatures(scde), 10)
p <- VariableFeaturePlot(scde)
LabelPoints(plot = p, points = top10, repel = TRUE)

# A standard pre-processing is the scaling of the normalized expression data so that the
# distribution of expression values per gene is centered around 0 across cells and 
# the variance per gene across cells is similar across genes so that genes with large expression
# changes are not over-emphasized.
# The scaled data are stored in the slot scde[["RNA"]]@scale.data
scde <- ScaleData(scde, features = rownames(scde))

# Now we can apply dimensional reduction methods such as the principal component analysis.
# The number of PCA components to be computed can be adjusted to the specific dataset and is set to 50 by default.
# Furthermore, the PCA is computed using only the set of highly variable features by default. This set 
# can be adjusted using the "features" parameter.
scde <- RunPCA(scde, npcs = 50)
# One way of visualizing the PCA result is to plot the loadings of features for the different dimensions.
VizDimLoadings(scde, dims = 1:2, reduction = "pca")
# We can also plot the cell arrangement in PCA space.
DimPlot(scde, reduction = "pca")

# A common heuristic method to examine the underlaying dimensionality of the dataset is to create a scree plot also 
# known as elbow plot which plots the principal components in order 1 to n against the explained variance The idea is 
# to retain only the first principal components with a large explained variance and to discard components 
# after the charateristic drop at the elbow.
ElbowPlot(scde)

# While PCA is a linear dimensional reduction technique, non-linear methods have become standard practice, such as UMAP and tSNE.
scde <- RunUMAP(scde, dims = 1:30)
scde <- RunTSNE(scde, dims = 1:30)
# Again, we can use the DimPlot method to visualize the low dimensional representations.
p1 <- DimPlot(scde, reduction = "umap")
p2 <- DimPlot(scde, reduction = "tsne")
p1 | p2

#########################################################################################################
# Analysis
#########################################################################################################

# We will perform a clustering of the cells using the Louvain algorithm. The necessary prerequisite
# is to identify the set of nearest neighbor cells, i.e. cells with similar transcriptomics profile.
scde <- FindNeighbors(scde, dims = 1:30)
# The Louvain algorithm is then executed, where the resolution parameter can be adjusted to modify the size and 
# thus number of obtained clusters.
scde <- FindClusters(scde, resolution = 1)

# The cluster assignments can now be visualized using the low dimensional reductions as shown above.
p1 <- DimPlot(scde, reduction = "umap", label = TRUE, repel = TRUE)
p1

# We can now compare this de-novo clustering with the provided cell type annotations stored in the meta data.
p2 <- DimPlot(scde, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cell.type')
p3 <- DimPlot(scde, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cluster.id')
p1 | p2 | p3



# Now we are at the point where we can look for marker genes. These are genes that exhibit differential 
# expression between a certain cluster and the remaining cells. Seurat has a convenience function to 
# test all clusters for differentially expressed genes against the remaining cells. Several thresholds 
# for the marker detection can be set, such as a minimal log-fold change (logfc.threshold), a minimum fraction 
# of cells in which a potential marker gene has been detected (min.pct), and which statistical test to use 
# to determine if the expression change is statistically significant (test.use). 
markerTable <- FindAllMarkers(scde, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")

# The resulting table contains the key information.
head(markerTable)

# We can count how many marker genes we find for each of the clusters using a 
# significance threshold of 0.01 and an absolute minimum log-fold change of 1, thereby including up- and 
# down-regulated genes.
markerTable %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.01, abs(avg_log2FC) > 1) %>%
    summarise(n = n())

# We can select the top 10 genes with the largest upregulation in the cluster.
top10.up <- markerTable %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.01) %>%
    top_n(n = 10, wt = avg_log2FC) 
# and the largest downregulation.
top10.down <- markerTable %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.01) %>%
    top_n(n = 10, wt = -avg_log2FC) 

# And plot the expression patterns as heatmap.
p1 <- DoHeatmap(scde, features = top10.up$gene) + NoLegend() + ggtitle('Top 10 upgregulated genes')
p2 <- DoHeatmap(scde, features = top10.down$gene) + NoLegend()+ ggtitle('Top 10 downgregulated genes')
p1 | p2


# There is a range of additional plots implemented in the Seurat package to facilitate data explorations such
# as VlnPlot, FeaturePlot, RidgePlot, CellScatter, DotPlot. 