---
title: "HOW-TO setup Google Colab to run Monocle3 (single cell RNA-seq)"
date: '2020-04-06'
layout: post
categories: Genomics
---

  https://colab.research.google.com/  

```python
# Load local files to google colab
from google.colab import files
uploaded = files.upload()
```



     <input type="file" id="files-7646b3da-897e-47e3-98ed-a3c84b4be686" name="files[]" multiple disabled />
     <output id="result-7646b3da-897e-47e3-98ed-a3c84b4be686">
      Upload widget is only available when the cell has been executed in the
      current browser session. Please rerun this cell to enable.
      </output>
      <script src="/nbextensions/google.colab/files.js"></script>


    Saving trajectory_donors_ipfs_MC_all_res_40.RData to trajectory_donors_ipfs_MC_all_res_40.RData



```python
!apt-get install libgdal-dev
!apt-get install
!apt-get update
!apt-get install libudunits2-dev
!apt-get update
```


```r
%%R
install.packages("units_0.6-6.tar.gz",configure.args='--with-udunits2-include=/usr/include/udunits2')
devtools::install_github("r-spatial/sf")
```


```r
%%R
install.packages("tidyverse")
install.packages("Matrix")
install.packages("cowplot")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
```


```r
%%R
library(monocle3)
load("your_file.RData")
```
