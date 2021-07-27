# En este script se hizo la consulta de los datos clinicos de pacientes con cancer
# de colon (COAD) y se elaboró la tabla de descripción de los datos para poder 
# empezar con los análisis.

library(TCGAbiolinks)
library(dplyr)
library(DT)

#### Query datos clínicos de COAD ####


## Obteniendo datos clínicos | 1

clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")

clinical %>%
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)


# Asegurando que se descargaron los datos correctos

#   En https://portal.gdc.cancer.gov/ se buscan aleatoriamente los id o alícuotas
#   de dos o tres pacientes para saber si los datos del proyecto concuerdan con la
#   base de datos del gdc. 
#   Posteriormente, en el portal se compara el número de muestras del proyecto y se 
#   compara con los datos clínicos obtenidos en el query. De esta manera se corrobora
#   que todo esté en orden.           # LISTO #


# Filtrando datos clínicos de importancia

#   Nos interesan especialmente los datos como prior malignancy, stage, status, entre
#   otros. Se espera quedar como con 10 datos clínicos.

# submitter_id 1, ajcc_pathologic_stage 3, tumor_stage 4, prior_malignancy 13,
# prior_treatment 15, state 16, gender 41, vital_status 42, age_at_index 43, 
# treatments_pharmaceutical_treatment_type 50, 

datos_clinicos <- clinical %>%
                  select(1, 3, 4, 13, 15, 16, 41, 42, 43, 50, 62)
                  
write.csv(datos_clinicos, "datos_clinicos.csv")

####

datos_clinicos <- read.csv("datos_clinicos.csv")

prueba <- datos_clinicos$submitter_id %>% unique()

####### Tablas de descripcion de los datos | RNA-Seq y miRNA-Seq #####

# Ahora trabajaremos con tablas de descripción de los datos, necesitamos obtener:

# grupos de muestras, cuántas/cuáles tienen RNA-seq,cuántas/cuáles tienen miRNASeq
# pareado con las de RNASeq. Éstas últimas serán con las que estemos trabajando.
    # Una vez que tengamos estas dos tablas, generamos las matrices de diseño
    # y empezamos ya con los análisis."


#--- Query RNA-Seq 

query2 <- GDCquery(project = "TCGA-COAD",
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling", #se necesita especificar
                      #data.type = "Gene Expression Quantification",
                      #workflow.type = "HTSeq - Counts"
                      experimental.strategy = "RNA-Seq")


#--- Query miRNA-Seq

query_mirnaseq <- GDCquery(project = "TCGA-COAD",
                   legacy = FALSE,
                   data.category = "Transcriptome Profiling",
                   experimental.strategy = "miRNA-Seq")    
    
####--- Diseño de las tablas, haremos DOS ---####

### Cuantas muestras tienen miRNA-Seq y/o RNA-Seq DENTRO del proyecto COAD

library(stringr)

#--- miRNAseq | numero de MUESTRAS con datos MIRNASEQ y sus ID | 1

mirnaseq <- query_mirnaseq[[1]][[1]] %>% 
            select(experimental_strategy, cases.submitter_id) %>% #930 MUESTRAS con miRNAseq REPETIDOS
            rename(miRNA_Seq = experimental_strategy) %>%
            mutate(miRNA_Seq = str_replace(miRNA_Seq, "miRNA-Seq","1")) %>%
            mutate(miRNA_Seq = as.numeric(miRNA_Seq)) %>%
            distinct(cases.submitter_id, .keep_all = TRUE) #444 muestras UNICAS


#--- RNAseq | no. de MUESTRAS con datos RNASEQ y sus ID | 2

rnaseq <- query2[[1]][[1]] %>%
          select(experimental_strategy, cases.submitter_id) %>% #1563 MUESTRAS con rnaseq REPETIDOS
          rename(RNA_Seq = experimental_strategy) %>%
          mutate(RNA_Seq = str_replace(RNA_Seq, "RNA-Seq","1")) %>%
          mutate(RNA_Seq = as.numeric(RNA_Seq)) %>% #preparandolo para el full_join
          distinct(cases.submitter_id, .keep_all = TRUE) #456 muestras UNICAS



#--- Diseño TABLA 1 | cuantas tienen RNAseq y/o mirnaseq

length(rownames(merge(rnaseq, mirnaseq, by = "cases.submitter_id"))) #441


x <- c("Muestras con RNA-Seq", "Muestras con miRNA-Seq", "Muestras pareadas")
y <- c(456, 444, 441)

tabla1 <- matrix(y, nrow  = 1) 
colnames(tabla1) <- x
rownames(tabla1) <- "Numero"



#--- Diseño TABLA 2 | CUALES muestras tienen RNAseq y/o mirnaseq

a <- c( "submitter_id", "RNA-Seq","miRNA-Seq", "Pareada" )


tabla2 <- rnaseq %>%
          full_join( mirnaseq,
               by = ("cases.submitter_id") ) %>% #se debe hacer un full_join
          mutate_if( is.numeric,coalesce,0 ) %>% 
          group_by( cases.submitter_id ) %>%
          mutate( Paired = if_else (( RNA_Seq == 1 & miRNA_Seq == 1), 1, 0))


tabla2 <- tabla2[,c(2,1,3,4)]  #3405 muestras en total REPETIDAS
                               #459 muestras UNICAS#

colnames(tabla2) <- a


#write.csv(tabla1, "tabla_descripcion_numerado")
write.csv(tabla2, "tabla_de_descripcion_pareados.csv")

# corroborando que hagan match con los pacientes de los datos clínicos
    
    


    
    
    
    
    








