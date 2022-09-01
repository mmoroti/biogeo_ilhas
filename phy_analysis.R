# Wed Aug 17 11:37:22 2022 ------------------------------
# Matheus de Toledo Moroti

# This new script is for analysis of anuran biodiversity in islands. In this script, I made the analysis of phyloregions using the phyloregion package.

# Loading packages
library(phyloregion)
library(phytools)
library(sf)
library(tidyverse)

# Loading data data
# Composition matrix
dir()
composition <- as_tibble(read.csv2("anura_islands_without0.csv"))
# put names in rows
head(composition)
composition.select <- composition %>% select(-X, -id)
rownames(composition.select) <- as.character(composition$id) 
View(composition.select)

# Traits
#traits_anura <- read.csv2("anura_islands_traits.csv", sep=";")

# Phylogeny
phy_anura_islands <- ape::read.tree("anura_islands.tre")

# Analysis of phylogenetic endemism and phyloregions 
#Sparse community matrix - requirement from analysis of phyloregion package
matrix_sparse_2 <- dense2sparse(composition.select)

# Takes a community data table and a (rooted) phylogenetic tree (with branch lengths) and calculates either strict or weighted endemism in Phylogenetic Diversity (PD). Strict endemism equates to the total amount of branch length found only in the sample/s and is described by Faith et al. (2004) as PD-endemism. Weighted endemism calculates the "spatial uniqueness" of each branch in the tree by taking the reciprocal of its range, multiplying by branch length and summing for all branch lengths present at a sample/site. Range is calculated simply as the total number of samples/sites at which the branch is present. This latter approach is described by Rosauer et al. (2009) as Phylogenetic endemism.
# phylo_endemism
pe <- phylo_endemism(matrix_sparse_2, phy_anura_islands)
head(pe)
# saving phylo_endemism values per islands
write.table(pe, "phylogenetic_endemism.txt", sep=",")

# Generating phylogenetic beta diversity
beta_diversity <- phylobeta(matrix_sparse_2, phy_anura_islands,  index.family = "jaccard")

# Cluster algorithm selection and validation
y <- select_linkage(beta_diversity[[1]])
barplot(y, horiz = TRUE, las = 1)

# Determining the optimal number of clusters
(d <- optimal_phyloregion(beta_diversity[[1]], method = "average", k = 20))
plot(d$df$k, d$df$ev, ylab = "Explained variances",
     xlab = "Number of clusters")
lines(d$df$k[order(d$df$k)], d$df$ev[order(d$df$k)], pch = 1)
points(d$optimal$k, d$optimal$ev, pch = 21, bg = "red", cex = 2)
points(d$optimal$k, d$optimal$ev, pch = 21, bg = "red", type = "h")

#?optimal_phyloregion
#?plot.phyloregion

###---
# Mon Aug 15 14:50:12 2022 ------------------------------
# Loading shapefile world islands 
dir()
shape <- sf::st_read("Anura_big islands_5000.shp")
shape <- sf:::as_Spatial(shape) 

# handling with geodata
# The column in shapefile needs to rename to "grids" for match with beta_diversity[[1]] object match
names_islands <- as_tibble(cbind(v1=shape$OBJECTID_1))
grids <- names_islands %>%  
  mutate(grids = paste0(v1)) 
shape <- cbind(shape,grids)
#dir()

# determing phyloregions with optimal number of clusters
phyloregions <- phyloregion(beta_diversity[[1]], k = 18, method = "average")
?phyloregion
phyloregions$evol_distinct
phyloregions$membership

# with geospatial data
phyloregions.shape <- phyloregion(beta_diversity[[1]], k = 18, method = "average", shp= shape) #coluna OBJECTID_1 tem os ids das 

plot_NMDS(phyloregions.shape, cex=6)
text_NMDS(phyloregions.shape, cex=2)

#plot(phyloregions, shp=philly_sf_merged, cex=1, palette="NMDS")

# Join in spatial dataset attributes if the function phyloregion with shape doesn't work
names(shape)
geo_data <- as_tibble(merge(shape, phyloregions$membership, by.x = "OBJECTID_1", by.y="grids"))

geo_data$cluster <- as.character(geo_data$cluster)
phylo_evoldist <- as_tibble(phyloregions$evol_distinct)
class(geo_data)

# join full data
geo_data_complete <- left_join(geo_data, phylo_evoldist, by = "cluster")
#Save
write.table(geo_data_complete, "evo_regions.txt", sep=",")
View(geo_data_complete)
