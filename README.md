ClusterID

What is it, how does it work, how well does it work?

It's R code.

It matches the expression of markers in your clusters to cell type definitions in a spreadsheet.

It's not amazing.

The basic approach is adapted from ScType:
https://github.com/IanevskiAleksandr/sc-type

More details to be published at https://www.colibri-cytometry.com/blog


Example data from OMIP-102 plotted with tSNE and scattermore using EmbedSOM to generate clusters:
![annotated tSNE](https://github.com/DrCytometer/ClusterID/blob/main/tsne_with_cluster_labels.jpg?raw=true)

Ridgeline plots showing expression of markers by cluster:
![annotated tSNE](https://github.com/DrCytometer/ClusterID/blob/main/cluster_histograms.jpg?raw=true)
