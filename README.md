# ClusterID

## What is it, how does it work, how well does it work?

It's R code.

It matches the expression of markers in your clusters to cell type definitions in a spreadsheet.

It's not amazing.

The basic approach is adapted from ScType:
https://github.com/IanevskiAleksandr/sc-type

More details to be published on [Colibri Cytometry](https://www.colibri-cytometry.com/post/what-s-that-cluster-part-ii)


Example data from OMIP-102 plotted with tSNE and scattermore using EmbedSOM to generate clusters:
![annotated tSNE](https://github.com/DrCytometer/ClusterID/blob/main/tsne_with_cluster_labels.jpg?raw=true)

Ridgeline plots showing expression of markers by cluster:
![annotated histograms](https://github.com/DrCytometer/ClusterID/blob/main/cluster_histograms.jpg?raw=true)
