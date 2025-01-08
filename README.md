# ClusterID

## What is it, how does it work, how well does it work?

It's R code.

It matches the expression of markers in your clusters to cell type definitions in a spreadsheet.

It's not amazing.

The basic approach is adapted from ScType:
https://github.com/IanevskiAleksandr/sc-type

More details on [Colibri Cytometry](https://www.colibri-cytometry.com/post/what-s-that-cluster-part-ii)

## Ways in which this might be improved (and I'm not going to right now):
1) Expand the celltype_database to cover more cells and markers.
2) Change the scoring system to provide better matching between phenotypes and cell definitions. Something like RMSD might be a good option.
3) Change the marker expression from positive/negative to scaled expression values (e.g., 0 - 1). This will have the disadvantage of mapping local minima and maxima to global minima and maxima in the celltype_database, but may generate better correspondance for many cell types, particularly for well characterized systems such as human blood.
4) Define cell types with boundaries in multi-dimensional space, then test whether a cluster falls within or overlaps with the bounds of a cell type. This requires mapping multi-dimensional expression patterns from real data (probably many datasets). Data from screens of CD marker expression would a good place to start.
5) Change the approach to something more complicated like a CNN.

Example data from OMIP-102 plotted with tSNE and scattermore using EmbedSOM to generate clusters:
![annotated tSNE](https://github.com/DrCytometer/ClusterID/blob/main/tsne_with_cluster_labels.jpg?raw=true)

Ridgeline plots showing expression of markers by cluster:
![annotated histograms](https://github.com/DrCytometer/ClusterID/blob/main/cluster_histograms.jpg?raw=true)
