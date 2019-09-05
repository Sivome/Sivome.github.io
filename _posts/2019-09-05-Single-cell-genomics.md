---
layout: post
title:  "Single cell genomics -- Work in progress"
date:   2019-09-05
categories: Genomics
---

Half-a-billion reads to sequence 10K cells with a median gene capture of 3K genes at 50K reads/cell is the current trend in single cell genomics. Most of the RNA-seq experiments focus on bulk RNA-seq methods. However, after closely looking at single cell datasets, the information obtained from such experiments can throw light on variety of underlying biological processes.

There are few known protocols for single-cell methods. For a detailed overview of the current techniques, refer to [1]. In this post, I will use a publicly available dataset, and go through the classification of the cell-types (using both supervised and unsupervised clustering techniques) using 2 R packages, "Seurat" and "SingleR".