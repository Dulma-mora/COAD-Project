#### density plot
library("readr")
library("dplyr")
library("igraph")

# paths

path_mir <- "COAD_mi.tsv"


# cargando matrices 

#coad_mi <- as.data.frame(readr::read_tsv(path_mir)) # opción que no carga

coad_mi <- vroom::vroom(path_mir)


#### nueva matriz de adyacencia 

#coad_mi <- mi_matrix

density(coad_mi)

quantile(coad_mi, 0.01)

y <- ifelse(coad_mi < quantile(coad_mi, 0.1), 0, 1)



# descubriendo que hace esto lol

g <- igraph::graph_from_incidence_matrix(y)

igraph::bipartite.mapping(g)



