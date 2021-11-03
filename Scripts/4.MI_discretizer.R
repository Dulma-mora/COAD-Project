#### 4. TESTEO | CALCULO DE M.I. Y DISCRETIZACION

library("readr")
library("infotheo")
library("dplyr")

# cargando funciones
# source("libs/functions_mi.R")


# paths

path_mir <- "data_mirnas.csv"
path_rna <- "data_mrna_transcriptome.csv"


# cargando matrices 

mir <- as.data.frame(readr::read_csv(path_mir)) # wow!
rna <- as.data.frame(readr::read_csv(path_rna))
# pregunta, por qué es mejor hacerlo así y no usar read.csv() ?


### DISCRETIZACION 

# mir
tempus <- proc.time()  #funcion nativa de R

d.mir <- par_discretizer(mir, korez = 10) #por qué 10?
tempus <- proc.time() - tempus
print(tempus)


# rna
tempus <- proc.time()

d.rna <- par_discretizer(rna, korez = 10)
tempus <- proc.time() - tempus
print(tempus)


### MI CALCULATING

tempus <- proc.time()
mirXrna <- par_mi_calc(sources = d.mir, 
                       targets = d.rna, 
                       korez = 10)
tempus <- proc.time() - tempus
print(tempus)

#creando matriz
mi_matrix <- bind_rows(!!!mirXrna, #explicit splicing
                       .id = "mirna/gen")

# write_tsv(mi_matrix, "COAD_mi.tsv")




