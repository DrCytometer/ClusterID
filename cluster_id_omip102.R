
# install any required packages

required_packages = c("dplyr", "readxl", "ggplot2", "ggridges", "data.table",
                      "scattermore", "tidyr")

for (req.package in required_packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    install.packages(req.package, repos='http://cran.us.r-project.org')
  }
}

library(dplyr)
library(readxl)
library(ggplot2)
library(ggridges)
library(data.table)
library(scattermore)
library(tidyr)

src.dir <- getwd()


# read in your cluster data-------
## see example data for expect format

cluster.data <- fread("omip102_example_data.csv", check.names = FALSE)

# remove any unwanted columns
input.data <- cluster.data %>% select(-c(Clust_color, tsneX, tsneY))

# reformat marker names


# summarize median expression by cluster
cluster.summary <- input.data %>%
  group_by(Cluster) %>%
  summarize_at(vars(everything()), median)

channel.labels <- colnames(cluster.summary)[-1]

# set the species used-------
# supported options are "Human" and "Mouse"
## to use a different species, e.g., Macaque, create new spreadsheets in the same format

species.used <- "Human"

cell.marker.filename <- paste0(species.used, "_marker_names.csv")
cell.type.filename <- paste0(species.used, "_celltype_database.xlsx")

cell.database <- read_xlsx( file.path(src.dir, cell.type.filename) )


# restrict cell types based on tissues used---------
# To identify tissue-specific cell types, tell the script which tissue sources
# you've used. Immune cells is the default, so include this in addition to any
# extra sources you've used.

tissue.options <- unique(cell.database$Tissue.restricted)
tissue.type <- c("Immune")

if ( length(tissue.type) > 1 ){
  cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)
} else if ( tissue.type == 0 ){
  cell.database <- cell.database
} else {
  cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)
}


# restrict cell ID based on known input-------
## this step is important to get accurate naming if you're only analyzing a subset
## of cells rather than, say, whole PBMC or all splenocytes.
## check the list of "cell.type.options" and set selected.cell.type appropriately

cell.type.options <- cell.database$Parent.cell.type

selected.cell.type <- c("All")


# match markers to database----------
## this section regularizes the naming of your markers, matching them to
## synonyms in the database. If you use a different spelling or convention,
## add that to the database.

marker.synonyms <- read.csv(file.path(src.dir, cell.marker.filename)) 
rownames(marker.synonyms) <- marker.synonyms$Marker.Name
marker.synonyms <- marker.synonyms %>% select(-Marker.Name)
colnames(marker.synonyms) <- NULL
marker.synonyms <- setNames(split(marker.synonyms, 
                                      seq(nrow(marker.synonyms))), 
                                rownames(marker.synonyms))
marker.synonyms <- lapply(marker.synonyms, as.list)

new.marker.names <- sapply( 1:length(channel.labels), function(x){
  pattern <- ifelse(grepl("\\(", channel.labels[x]), 
                gsub("([()])","\\\\\\1", channel.labels[x]), 
                paste0( channel.labels[x], "\\b"))
  names(marker.synonyms)[grep(pattern, marker.synonyms,
                                 ignore.case = TRUE)]
} )


# scale data for cluster matching--------

scaled.cluster.summary <- scale(cluster.summary[,-1], scale = FALSE)[,]

colnames(scaled.cluster.summary) <- new.marker.names


# prepare marker lists and match clusters-----------

source( file.path( src.dir, "prepare_marker_lists.r" ) )

source( file.path( src.dir, "flow_cluster_id_score.r" ) )

marker.list <- prepare_marker_lists( file.path( src.dir, cell.type.filename ), 
                                     tissue.type, selected.cell.type )

cluster.labels <- paste0( "Cluster_", cluster.summary$Cluster)


# match clusters to cell type database

scaled.id.score <- flow_cluster_id_score(t(scaled.cluster.summary),
                                         marker_pos = marker.list$markers_positive, 
                                         marker_neg = marker.list$markers_negative )

for (cluster in 1:ncol(scaled.id.score)){
  cluster.labels[cluster] <- names( which.max(scaled.id.score[,cluster]) )
}

names(cluster.labels) <- cluster.summary$Cluster
cluster.data$Cluster_label <- cluster.labels[cluster.data$Cluster]

# optional: save csv with cluster labels
fwrite(cluster.data, file = "labeled_cluster_data.csv")

# plot heatmaps of cluster ID scores--------

png( file.path(src.dir, "cluster_id_heatmap.png"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, Rowv = NA, Colv = NA, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()

png( file.path(src.dir, "cluster_id_heatmap_dendro.png"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()



# ridgeline plots of the clusters to check how well the labels match the phenotypes
cluster.summary$Cluster_label <- cluster.labels[cluster.summary$Cluster]

cluster.data.long <- cluster.data %>%
  select(-tsneX, -tsneY) %>%
  pivot_longer( -c(Cluster, Cluster_label, Clust_color), names_to = "Marker", values_to = "Expression")

ggplot(cluster.data.long, aes(x = Expression, y = Marker, fill = Clust_color)) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.001) +
  facet_wrap(~ Cluster, labeller = as_labeller(cluster.summary$Cluster_label)) +
  theme_ridges() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8))

ggsave("cluster_histograms.jpg", width = 20, height = 50, limitsize = FALSE)


# if your data contains tSNE or UMAP coordinates, plot that
ggplot(cluster.data, aes(x = tsneX, y = tsneY, color = Clust_color)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=cluster.data$Clust_color, 
                        labels=cluster.data$Cluster_label,
                        guide = "legend") +
  theme_classic()

ggsave("tsne_with_cluster_labels.jpg", width = 18, height = 10, limitsize = FALSE)







# Optional: differentiate similar clusters by variable markers----------------

# find (positions of) clusters with non-unique names
non.unique.cluster.names <- cluster.labels[duplicated(cluster.labels)]
# sort into groups
non.unique.groups <- unique(non.unique.cluster.names)

# append high variance marker labels to all clusters until unique

for (cluster.group in non.unique.groups) {
  
  # collect clusters needing renaming
  data.temp <- cluster.summary[ grepl( cluster.group, cluster.summary$Cluster_label),]
  data.temp <- data.temp %>% select(-Cluster, -Cluster_label)
  # find most variable channels
  variances <- apply(data.temp, 2, var)
  sorted.variances <- sort(variances, decreasing = TRUE)
  
  # find min and max (neg and pos) for each channel
  max.expression <- apply(data.temp, 2, max)
  min.expression <- apply(data.temp, 2, min)
  
  for( cluster in 1:nrow(data.temp)){
    # get position of cluster in full list
    position.in.cluster.list <- cluster.summary$Cluster_label[
      grepl(cluster.group, cluster.summary$Cluster_label)][cluster]
    
    # append markers to name for each positive variable up to the number of clusters in the group or max 4
    
    markers.to.append.n <- ifelse( nrow(data.temp) < 5, nrow(data.temp), 4 )
    
    for( marker in names(sorted.variances)[1:markers.to.append.n]){
      cluster.summary$Cluster_label[as.numeric(names(position.in.cluster.list))] <- 
        if ( abs(cluster.summary[as.numeric(names(position.in.cluster.list)),marker] - max.expression[marker] ) < abs(cluster.summary[as.numeric(names(position.in.cluster.list)),marker] - min.expression[marker] )){
          paste( cluster.summary$Cluster_label[as.numeric(names(position.in.cluster.list))], marker )
        } else{
          cluster.summary$Cluster_label[as.numeric(names(position.in.cluster.list))] <- cluster.summary$Cluster_label[as.numeric(names(position.in.cluster.list))]
        }
    }
  }
}

# add a cluster name column to the original data
cluster.data$Cluster_label <- cluster.summary$Cluster_label[cluster.data$Cluster]

# optional: save csv with cluster labels
fwrite(cluster.data, file = "long_labeled_cluster_data.csv")


# ridgeline plots of the clusters to check how well the labels match the phenotypes

cluster.data.long <- cluster.data %>%
  select(-tsneX, -tsneY) %>%
  pivot_longer( -c(Cluster, Cluster_label, Clust_color), names_to = "Marker", values_to = "Expression")

ggplot(cluster.data.long, aes(x = Expression, y = Marker, fill = Clust_color)) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.001) +
  facet_wrap(~ Cluster, labeller = as_labeller(cluster.summary$Cluster_label)) +
  theme_ridges() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8))

ggsave("cluster_histograms_long_labels.jpg", width = 20, height = 50, limitsize = FALSE)


# if your data contains tSNE or UMAP coordinates, plot that
ggplot(cluster.data, aes(x = tsneX, y = tsneY, color = Clust_color)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=cluster.data$Clust_color, 
                        labels=cluster.data$Cluster_label,
                        guide = "legend") +
  theme_classic()

ggsave("tsne_with_cluster_long_labels.jpg", width = 18, height = 10, limitsize = FALSE)
