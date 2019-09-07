---
layout: post
title:  "Single cell genomics -- Work in progress"
date:   2019-09-05
categories: Genomics
---

Half-a-billion reads to sequence 10K cells with a median gene capture of 3K genes at 50K reads/cell is the current trend in single cell genomics. Most of the RNA-seq experiments focus on bulk RNA-seq methods. However, after closely looking at single cell datasets, the information obtained from such experiments can throw light on variety of underlying biological processes.

There are few known protocols for single-cell methods. For a detailed overview of the current experimental techniques, refer to [1]. I will focus on the computational aspect of the single cell data using the available R packages "Seurat" and "SingleR".

More R vignettes on Seurat can be found here: https://satijalab.org/seurat/vignettes.html
SingleR homepage can be found here: http://comphealth.ucsf.edu/sample-apps/SingleR/SingleR_create.html
Publication on usage of "SingleR" along with R script used to generate the figures can be found here: 

Mouse cell atlas (MCA): http://bis.zju.edu.cn/MCA/
Seurat R vignette on analyzing MCA can be found here: https://satijalab.org/seurat/v3.0/mca.html


```{r}
# required packages
library(tidyverse)
library(Seurat)
library(SingleR)

# Download the MCA processed data i.e., RDS file directly from the above link
library(Seurat)
mca.matrix <- readRDS(file = "C:/Sivome/SingleCellGenomics/MCA/MCA/data/MCA_merged_mat.rds")
mca.metadata <- read.csv(file = "C:/Sivome/ingleCellGenomics/MCA/MCA/MCA_All-batch-removed-assignments.csv", row.names = 1)
```
