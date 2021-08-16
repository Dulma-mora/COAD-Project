####### Tablas de descripcion de los datos | RNA-Seq y miRNA-Seq #####

# Ahora trabajaremos con tablas de descripción de los datos, necesitamos obtener:

# grupos de muestras, cuántas/cuáles tienen RNA-seq,cuántas/cuales tienen miRNASeq
# pareado con las de RNASeq. Éstas últimas serán con las que estemos trabajando.
# Una vez que tengamos estas dos tablas, generamos las matrices de diseño
# y empezamos ya con los análisis."


#--- TUMORES con RNA-Seq 

#query general
query2 <- GDCquery(project = "TCGA-COAD",
                   legacy = FALSE,
                   data.category = "Transcriptome Profiling", #se necesita especificar
                   #data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts", #esto se ocupa durante la expresion diferencial
                   experimental.strategy = "RNA-Seq")
                   #sample.type = "Primary Tumor")


rnaseq <- (query2$results)
rnaseq <- as.data.frame(rnaseq[[1]])

tumores_rnaseq <- filter(rnaseq, sample_type == "Primary Tumor")

### unicos
length(unique(tumores_rnaseq$cases.submitter_id))
length(tumores_rnaseq$cases.submitter_id)



#--- TUMORES con miRNA-Seq

#query general
query_mirnaseq <- GDCquery(project = "TCGA-COAD",
                           legacy = FALSE,
                           data.category = "Transcriptome Profiling",
                           experimental.strategy = "miRNA-Seq",
                           data.type = "Isoform Expression Quantification") # usamos isoformas y no mature mirna


mirnaseq <- (query_mirnaseq$results)
mirnaseq <- as.data.frame(mirnaseq[[1]])

tumores_mirnaseq <- filter(mirnaseq, sample_type == "Primary Tumor")
length(unique(tumores_mirnaseq$cases.submitter_id)) #444 tumores con mirnaseq


### cuantas muestras SANAS hay con mirnaseq

sanos_mirnaseq <- filter(mirnaseq, sample_type == "Solid Tissue Normal")
length(unique(sanos_mirnaseq$cases.submitter_id)) # hay 8


### muestras SANAS con rnaseq

sanos_rnaseq <- filter(rnaseq, sample_type == "Solid Tissue Normal")
length(unique(sanos_rnaseq$cases.submitter_id)) #41 sanos con rnaseq


#--- TUMORES con rnaseq Y mirnaseq

length(intersect(tumores_rnaseq$cases.submitter_id,
                 tumores_mirnaseq$cases.submitter_id)) # 441 tumores con ambos

#--- SANOS con rnaseq Y mirnaseq

length(intersect(sanos_rnaseq$cases.submitter_id,
                 sanos_mirnaseq$cases.submitter_id)) # hay 8


###----- Diseñando tablas ----- ###

library(stringr)

# TABLA DATOS DE TUMOR

datos_mirnaseq_tumor <- tumores_mirnaseq %>% 
     select(cases.submitter_id, experimental_strategy) %>% 
     rename(miRNA_Seq = experimental_strategy) %>%
     #rename(Primary_Tumor = sample_type) %>%
     mutate(miRNA_Seq = str_replace(miRNA_Seq, "miRNA-Seq","1")) %>%
     #mutate(Primary_Tumor = str_replace(Primary_Tumor, "Primary Tumor","1"))
     mutate(miRNA_Seq = as.numeric(miRNA_Seq)) %>%
     #mutate(Primary_Tumor = as.numeric(Primary_Tumor)) %>%
     distinct(cases.submitter_id, .keep_all = TRUE)


datos_rnaseq_tumor <- tumores_rnaseq %>%
     select(cases.submitter_id, experimental_strategy) %>% 
     rename(RNA_Seq = experimental_strategy) %>%
     #rename(Primary_Tumor = sample_type) %>%
     mutate(RNA_Seq = str_replace(RNA_Seq, "RNA-Seq","1")) %>%
     #mutate(Primary_Tumor = str_replace(Primary_Tumor, "Primary Tumor","1"))
     mutate(RNA_Seq = as.numeric(RNA_Seq)) %>% 
     #mutate(Primary_Tumor = as.numeric(Primary_Tumor)) %>%
     distinct(cases.submitter_id, .keep_all = TRUE)

# haciendo full_join

tabla_tumores <- datos_rnaseq_tumor %>%
           full_join( datos_mirnaseq_tumor,
                by = ("cases.submitter_id") ) %>% #se debe hacer un full_join
           mutate_if( is.numeric,coalesce,0 ) %>% 
           group_by( cases.submitter_id ) %>%
           mutate( Paired = if_else (( RNA_Seq == 1 & miRNA_Seq == 1), 1, 0))

write.csv(tabla_tumores, "tabla_tumores.csv")


### ---- TABLA DATOS SANOS

datos_mirnaseq_sanos <- sanos_mirnaseq %>%
                        select(cases.submitter_id, experimental_strategy) %>%
                        rename(miRNA_Seq = experimental_strategy) %>%
                        #rename(Solid_Tissue_Normal = sample_type) %>%
                        mutate(miRNA_Seq = str_replace(miRNA_Seq, "miRNA-Seq","1")) %>%
                        #mutate(Solid_Tissue_Normal = str_replace(Solid_Tissue_Normal, "Solid Tissue Normal","1"))
                        mutate(miRNA_Seq = as.numeric(miRNA_Seq)) %>%
                        #mutate(Solid_Tissue_Normal = as.numeric(Solid_Tissue_Normal)) %>%
                        distinct(cases.submitter_id, .keep_all = TRUE)

                        
datos_rnaseq_sanos <- sanos_rnaseq %>%
     select(cases.submitter_id, experimental_strategy) %>% 
     rename(RNA_Seq = experimental_strategy) %>%
     mutate(RNA_Seq = str_replace(RNA_Seq, "RNA-Seq","1")) %>%
     mutate(RNA_Seq = as.numeric(RNA_Seq)) %>% 
     distinct(cases.submitter_id, .keep_all = TRUE)


# haciendo full join

tabla_sanos <- datos_rnaseq_sanos %>%
     full_join( datos_mirnaseq_sanos,
                by = ("cases.submitter_id") ) %>% #se debe hacer un full_join
     mutate_if( is.numeric,coalesce,0 ) %>% 
     group_by( cases.submitter_id ) %>%
     mutate( Paired = if_else (( RNA_Seq == 1 & miRNA_Seq == 1), 1, 0)) %>% as.data.frame()

write.csv(tabla_sanos, "tabla_sanos.csv")


### ------------------------- Tablas resultantes

datos_clinicos <- read.csv("datos_clinicos.csv")
datos_clinicos <- datos_clinicos[,-1]

tabla_tumores <- read.csv("tabla_tumores.csv")
tabla_tumores <- tabla_tumores[,-1]

tabla_sanos <- read.csv("tabla_sanos.csv")
tabla_sanos <- tabla_sanos[,-1]


# cambiando los colnames de las tablas para que coincidan con los datos clinicos

colnames(tabla_tumores)[1] <- "submitter_id"
colnames(tabla_sanos)[1] <- "submitter_id"



#### corroborando que los valores de las tablas estén en la lista de datos clinicos


#--- viendo que los datos tumorales estén en los datos clínicos

(unique(datos_clinicos$submitter_id)) %in%(unique(tabla_tumores$cases.submitter_id))
# hay dos de datos clinicos que NO están en la tabla de tumorales 
# correcto!! tenemos 461 pacientes en total y 459 pacientes con tumor
# concuerda con que hay pacientes con muestras tumorales y tejido sano
# es decir, solo hay DOS pacientes en el proyecto son UNICAMENTE datos de tejido sano


#--- viendo que los datos sanos estén en los datos clinicos

(unique(tabla_sanos$cases.submitter_id)) %in% (unique(datos_clinicos$submitter_id))
 #siii

(unique(tabla_tumores$cases.submitter_id)) %in% (unique(datos_clinicos$submitter_id))
# siii x2


#### añadiendo primary diagnosis a las tablas resultantes

#tumores

tabla_tumores <- merge( x = tabla_tumores, y = datos_clinicos[c(1,4)],
                        by = ("submitter_id"))

#sano

tabla_sanos <- merge( x = tabla_sanos, y = datos_clinicos[c(1,4)],
                      by = ("submitter_id"))


write.csv(tabla_tumores, "tabla_tumores.csv")
write.csv(tabla_sanos, "tabla_sanos.csv")

