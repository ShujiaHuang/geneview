# genomation

A toolkit to summarize, annotate and visualize genomic intervals.

## Welcome to genomation

This is an R package that contains a collection of tools for visualizing and analyzing genome-wide data sets. The package works with a variety of genomic interval file types and enables easy summarization and annotation of high throughput data sets with given genomic annotations.

## Features

Ability to work with multiple flat file formats such as BED, GFF, BAM, bigWig and generic tabular text files
Annotation of genomic intervals. For example, you can see what percentage of your genomic intervals overlaps with exon/intron/promoter annotation.
Summary of genomic scores or read coverages over pre-defined regions. The pre-defined regions can be equi-width such as region around TSS or could be a set of regions with varying widths such as CpG islands. This operation will result in a matrix or set of matrices containing scores for each pre-defined region
Visualisation of summary matrices such as meta-region(meta-gene, meta-promoter, etc.) plots or heatmaps. These functions can employ K-means clustering of rows of summary matrices as well.


# 发表文章可以参考的地方：

## Visualizing Genomic Data Using Gviz and Bioconductor

这是另一篇值得参考的文献：

https://sci-hub.st/https://link.springer.com/protocol/10.1007/978-1-4939-3578-9_16

## Giotto: a toolbox for integrative analysis and visualization of spatial expression data

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2
