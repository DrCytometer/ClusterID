library(dplyr)
library(Rtsne)
library(EmbedSOM)
library(ggplot2)
library(scattermore)
library(data.table)
library(viridis)
library(ConsensusClusterPlus)


omip102.data <- data.table::fread("OMIP102_Single Cells.csv")
colnames(omip102.data) <- gsub( "-A", "", colnames(omip102.data) )

n.clusters <- 60

date.seed <- 20240520

set.seed(date.seed)

flow.som <- EmbedSOM::SOM(omip102.data, xdim = 24, 
                          ydim = 24, batch = TRUE,
                          parallel = TRUE, threads = 0 )

# get clusters
flow.som.mapping <- flow.som$mapping[ , 1 ]

# get clusters from som mapping
consensus.cluster <- ConsensusClusterPlus( t( flow.som$codes ),
                                           maxK = n.clusters, reps = 100, pItem = 0.9, pFeature = 1,
                                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                           distance = "euclidean", 
                                           seed = set.seed(date.seed) )

flow.som.event.cluster <- consensus.cluster[[ n.clusters ]]$
  consensusClass[ flow.som.mapping ]

# reorder clusters from bigger to smaller
flow.som.cluster.rank <- 1 + n.clusters - 
  rank( table( flow.som.event.cluster ), ties.method = "last" )
flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
names( flow.som.event.cluster ) <- NULL

# set clusters as a factor

omip102.data$Cluster <- factor( flow.som.event.cluster, 
                             levels = 1 : n.clusters )

# set cluster color
color.palette <- viridis(n.clusters)
color.palette <- sample(color.palette, 60)
omip102.data$Clust_color <- color.palette[flow.som.event.cluster]

set.seed( date.seed )

tsne.result <- Rtsne( omip102.data[,1:50], perplexity = 30, 
                      exaggeration_factor = 4,
                      eta = nrow(omip102.data)/4,
                      max_iter = 750, stop_lying_iter = 75,
                      check_duplicates = FALSE, pca = FALSE, 
                      num_threads = 0 )

omip102.data$tsneX <- tsne.result$Y[,1]
omip102.data$tsneY <- tsne.result$Y[,2]

ggplot(omip102.data, aes(x = tsneX, y = tsneY, color = Clust_color)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=omip102.data$Clust_color, 
                        labels=omip102.data$Cluster,
                        guide = "legend") +
  theme_classic()

fwrite(omip102.data, file = "omip102_example_data.csv")
