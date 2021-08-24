####################
#discretization and 
#mutual information 
#calculation
#functions
####################

#packrat unbundle
#packrat::on()
#packrat::unbundle(bundle = "packrat/bundles/cdre-miR-BrCanSub-2019-04-10.tar.gz", where = ".")
#packrat::restore()

#libraries
library("infotheo")
library("dplyr")
library("readr")

#functions 

#discretization

discretizer <- function(expmatrix, diskz = "equalfreq"){
     #names 
     featureNames <- dplyr::pull(expmatrix, 1)
     sampleNames  <- colnames(expmatrix[,-1])
     
     my_index <- 1:nrow(expmatrix)
     #names(my_index) <- dplyr::pull(expmatrix, 1)
     names(my_index) <- featureNames
     
     my_discrete <- lapply(X = my_index, FUN = function(i){
          my_row <- expmatrix[i,-1] %>% as.numeric()
          names(my_row) <- sampleNames
          my_row <- infotheo::discretize(X = my_row, disc = diskz)
          my_row <- as.data.frame(t(my_row))
          colnames(my_row) <- sampleNames
          return(my_row)
          
     })
     return(my_discrete)
}

par_discretizer <- function(expmatrix, korez = 2, diskz = "equalfreq"){
     #names 
     featureNames <- dplyr::pull(expmatrix, 1)
     sampleNames  <- colnames(expmatrix[,-1])
     
     my_index <- 1:nrow(expmatrix)
     #names(my_index) <- dplyr::pull(expmatrix, 1)
     names(my_index) <- featureNames
     
     my_discrete <- parallel::mclapply(X = my_index, mc.cores = korez, FUN = function(i){
          my_row <- expmatrix[i,-1] %>% as.numeric()
          names(my_row) <- sampleNames
          my_row <- infotheo::discretize(X = my_row, disc = diskz)
          my_row <- as.data.frame(t(my_row))
          colnames(my_row) <- sampleNames
          return(my_row)
          
     })
     return(my_discrete)
}

#MI calculation 

mi_calc <- function(sources, targets, methodz="emp"){
     lapply(X = sources, FUN = function(i){
          sapply(X = targets, FUN = function(j){
               ni <- as.numeric(i)
               nj <- as.numeric(j)
               infotheo::mutinformation(X = ni, 
                                        Y = nj,
                                        method = methodz)
          })
     })
}

par_mi_calc <- function(sources, targets, korez = 2, methodz="emp"){
     parallel::mclapply(X = sources, mc.cores = korez, FUN = function(i){
          sapply(X = targets, FUN = function(j){
               ni <- as.numeric(i)
               nj <- as.numeric(j)
               infotheo::mutinformation(X = ni, 
                                        Y = nj,
                                        method = methodz)
          })
     })
}