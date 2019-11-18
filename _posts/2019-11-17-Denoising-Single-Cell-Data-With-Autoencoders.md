---
title: "2019-11-17-Denoising-Single-Cell-Data-With-Autoencoders"
output: html_document
---




```
## -- Attaching packages ---------------------------------- tidyverse 1.2.1 --
```

```
## v ggplot2 3.2.1     v purrr   0.3.3
## v tibble  2.1.3     v dplyr   0.8.3
## v tidyr   1.0.0     v stringr 1.4.0
## v readr   1.3.1     v forcats 0.4.0
```

```
## -- Conflicts ------------------------------------- tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```


```r
# Download the MCA processed data i.e., RDS file directly from the link found in Seurat MCA guided clustering
mca.matrix <- readRDS(file = "C:/Sivome/SingleCellGenomics/MCA/MCA/MCA_merged_mat.rds")
mca.metadata <- read.csv(file = "C:/Sivome/SingleCellGenomics/MCA/MCA/MCA_All-batch-removed-assignments.csv", row.names = 1)
```

Download the MCA processed data i.e., RDS file directly from the link found in Seurat MCA guided clustering and here I stored them in mca.matrix and mca.metadata.

```r
# Quick look at the meta data
head(mca.metadata)
```

```
##                               ClusterID  Tissue     Batch
## Bladder_1.AAAACGAAAACGGGGCGA  Bladder_1 Bladder Bladder_1
## Bladder_1.AAAACGAAGCGGCCGCTA  Bladder_5 Bladder Bladder_1
## Bladder_1.AAAACGAAGTACTAGCAT Bladder_16 Bladder Bladder_1
## Bladder_1.AAAACGACGTTGCTGTGT  Bladder_8 Bladder Bladder_1
## Bladder_1.AAAACGAGCGAGCGAGTA  Bladder_4 Bladder Bladder_1
## Bladder_1.AAAACGAGGGTCAGATGG  Bladder_7 Bladder Bladder_1
##                                    Cell.Barcode
## Bladder_1.AAAACGAAAACGGGGCGA AAAACGAAAACGGGGCGA
## Bladder_1.AAAACGAAGCGGCCGCTA AAAACGAAGCGGCCGCTA
## Bladder_1.AAAACGAAGTACTAGCAT AAAACGAAGTACTAGCAT
## Bladder_1.AAAACGACGTTGCTGTGT AAAACGACGTTGCTGTGT
## Bladder_1.AAAACGAGCGAGCGAGTA AAAACGAGCGAGCGAGTA
## Bladder_1.AAAACGAGGGTCAGATGG AAAACGAGGGTCAGATGG
```


```r
# Seurat analyses (Most of the content in this block is directly copied from the Seurat Vignette)
dim(mca.matrix)
```

```
## [1]  39855 405191
```

```r
# Extract first 10K cells
mca.matrix.10K <- mca.matrix[,1:10000]

mca.10K <- CreateSeuratObject(counts = mca.matrix.10K, meta.data = mca.metadata, project = "MouseCellAtlas")
```

```
## Warning in CreateSeuratObject(counts = mca.matrix.10K, meta.data =
## mca.metadata, : Some cells in meta.data not present in provided counts
## matrix.
```

```
## Warning: Feature names cannot have underscores ('_'), replacing with dashes
## ('-')
```

```r
# Normalize data
mca.10K <- NormalizeData(mca.10K, normalization.method = "LogNormalize", scale.factor = 10000)
# Variable features in the data
mca.10K <- FindVariableFeatures(mca.10K)

# Mitochondrial expression is not used in the variable feature list for classification
# Mitochondrial genes generally start with 'mt-', but sometimes they could be in different case as well
mca.10K[["percent.mt"]] <- PercentageFeatureSet(mca.10K, pattern = "^mt-")
mca.10K <- ScaleData(mca.10K, vars.to.regress = "percent.mt")
```

```
## Regressing out percent.mt
```

```
## Centering and scaling data matrix
```


```r
mca.10K.original <- read.csv("C:\\Sivome\\SingleCellGenomics\\MCA\\MCA\\mca_counts_norm_lc.csv")
rownames(mca.10K.original) <- mca.10K.original$X
mca.10K.original <- mca.10K.original[, -1]

 
mca.10K.autoencoder <- read.csv("C:\\Sivome\\SingleCellGenomics\\MCA\\MCA\\mca_counts_norm_lc_output_pd_colnames.csv")
mca.10K.autoencoder <- mca.10K.autoencoder[, -1]
rownames(mca.10K.autoencoder) <- rownames(mca.10K.original)

mca.10K.original <- data.matrix(mca.10K.original)
mca.10K.autoencoder <- data.matrix(mca.10K.autoencoder)
```


```r
# duplicate mca.10K seurat objects and replace the scale.data with the new original and autoencoder analyses
mca.10K.original.seurat <- mca.10K
mca.10K.autoencoder.seurat <- mca.10K
mca.10K.original.seurat@assays$RNA@scale.data <- mca.10K.original
mca.10K.autoencoder.seurat@assays$RNA@scale.data <- mca.10K.autoencoder
```


```r
mca.10K.original.seurat <- RunPCA(mca.10K.original.seurat, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
```

```
## PC_ 1 
## Positive:  Ighv2-2, Ighv1-22, Ighv10-3, Ighv3-6, Ighv1-81 
## Negative:  Retnlg, S100a6, mt-Cytb, Rps27, Rplp1 
## PC_ 2 
## Positive:  Retnlg, Elane, Mpo, Fcnb, Igkc 
## Negative:  Mgp, Dcn, Gstm1, Gsn, S100a6 
## PC_ 3 
## Positive:  Retnlg, Mgp, Dcn, Gsn, Col1a2 
## Negative:  Gstm1, Elane, Rplp0, Rplp1, mt-Cytb 
## PC_ 4 
## Positive:  Retnlg, S100a6, Gstm1, Gsta4, Ly6d 
## Negative:  Elane, Mpo, Mgp, Dcn, Gsn 
## PC_ 5 
## Positive:  Crip1, Rplp0, Rplp1, Rps27, Rps29 
## Negative:  Elane, Gstm1, S100a6, Mpo, Prtn3
```

```r
mca.10K.original.seurat <- FindNeighbors(mca.10K.original.seurat, reduction = "pca", dims = 1:75, nn.eps = 0.5)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
mca.10K.original.seurat <- FindClusters(mca.10K.original.seurat, resolution = 3, n.start = 10)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 10000
## Number of edges: 486341
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.7731
## Number of communities: 36
## Elapsed time: 1 seconds
```

```r
#t-SNE, a populare classification method
mca.10K.original.seurat <- RunTSNE(mca.10K.original.seurat, dims = 1:75)

#UMAP is relatively new and with some datasets, it is shown to perform better than t-SNE
mca.10K.original.seurat <- RunUMAP(mca.10K.original.seurat, dims = 1:75, min.dist = 0.75)
```

```
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
```

```
## 20:07:02 UMAP embedding parameters a = 0.2734 b = 1.622
```

```
## 20:07:02 Read 10000 rows and found 75 numeric columns
```

```
## 20:07:02 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 20:07:02 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 20:07:05 Writing NN index file to temp file C:\Users\Viswa\AppData\Local\Temp\Rtmpm6apjQ\filee7043871a7
## 20:07:05 Searching Annoy index using 1 thread, search_k = 3000
## 20:07:08 Annoy recall = 100%
## 20:07:08 Commencing smooth kNN distance calibration using 1 thread
## 20:07:09 Initializing from normalized Laplacian + noise
## 20:07:11 Commencing optimization for 500 epochs, with 482572 positive edges
## 20:07:36 Optimization finished
```


```r
> DimPlot(mca.10K.original.seurat, group.by = 'orig.ident') + DarkTheme()
> DimPlot(mca.10K.autoencoder.seurat, group.by = 'orig.ident') + DarkTheme()
```

```
## Error: <text>:1:1: unexpected '>'
## 1: >
##     ^
```




```r
mca.10K.autoencoder.seurat <- RunPCA(mca.10K.autoencoder.seurat, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
```

```
## PC_ 1 
## Positive:  Gm12583, Llph-ps1, Gm5356, Arg1, Gm4925 
## Negative:  Mgp, Gsn, Sparc, Col1a2, Col3a1 
## PC_ 2 
## Positive:  Mgp, Gsn, Col1a2, Sparc, Col3a1 
## Negative:  Plac8, Rplp0, mt-Cytb, Rplp1, Mpo 
## PC_ 3 
## Positive:  Gsta4, Sprr1a, Krt15, Wfdc2, Ly6d 
## Negative:  Dcn, Col1a2, Col3a1, Mgp, Elane 
## PC_ 4 
## Positive:  Gsta4, Sprr1a, Krt15, Mgp, Wfdc2 
## Negative:  Acta2, Myl9, Igfbp7, Tpm2, Hspb1 
## PC_ 5 
## Positive:  Elane, Plac8, Rplp0, Mpo, mt-Cytb 
## Negative:  Hbb-bt, Hba-a1, Car2, Car1, Blvrb
```

```r
mca.10K.autoencoder.seurat <- FindNeighbors(mca.10K.autoencoder.seurat, reduction = "pca", dims = 1:75, nn.eps = 0.5)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
mca.10K.autoencoder.seurat <- FindClusters(mca.10K.autoencoder.seurat, resolution = 3, n.start = 10)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 10000
## Number of edges: 300384
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8062
## Number of communities: 41
## Elapsed time: 0 seconds
```

```r
#t-SNE, a populare classification method
mca.10K.autoencoder.seurat <- RunTSNE(mca.10K.autoencoder.seurat, dims = 1:75)

#UMAP is relatively new and with some datasets, it is shown to perform better than t-SNE
mca.10K.autoencoder.seurat <- RunUMAP(mca.10K.autoencoder.seurat, dims = 1:75, min.dist = 0.75)
```

```
## 20:08:12 UMAP embedding parameters a = 0.2734 b = 1.622
```

```
## 20:08:12 Read 10000 rows and found 75 numeric columns
```

```
## 20:08:12 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 20:08:12 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 20:08:14 Writing NN index file to temp file C:\Users\Viswa\AppData\Local\Temp\Rtmpm6apjQ\filee70577043cf
## 20:08:14 Searching Annoy index using 1 thread, search_k = 3000
## 20:08:16 Annoy recall = 100%
## 20:08:17 Commencing smooth kNN distance calibration using 1 thread
## 20:08:18 Initializing from normalized Laplacian + noise
## 20:08:18 Commencing optimization for 500 epochs, with 395680 positive edges
## 20:08:41 Optimization finished
```
