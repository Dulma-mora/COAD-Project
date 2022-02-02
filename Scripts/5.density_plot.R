#### density plot
library("readr")
library("dplyr")
library("igraph")

# paths

path_mir <- "COAD_mi.tsv"


# cargando matrices 

coad_mi <- as.data.frame(readr::read_tsv(path_mir)) # wow!

# coad_mi <- read.table(file = 'COAD_mi.tsv', sep = '\t', header = TRUE) IGNORAR


#### nueva matriz de adyacencia 

coad_mi <- mi_matrix

density(coad_mi)

quantile(coad_mi, 0.01)

y <- ifelse(coad_mi < quantile(coad_mi, 0.1), 0, 1)



# descubriendo que hace esto lol

g <- igraph::graph_from_incidence_matrix(y)

igraph::bipartite.mapping(g)



