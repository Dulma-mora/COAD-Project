#### 3. MATRICES DE EXPRESIÓN MIRNAS Y RNA 

library(TCGAbiolinks)
library(dplyr)


######## consulta de datos de expresión de miRNAs 

query_mirna <- GDCquery(project = "TCGA-COAD",
                   legacy = FALSE,
                   data.category = "Transcriptome Profiling", 
                   data.type = "Isoform Expression Quantification",
                   #workflow.type = "HTSeq - Counts", 
                   experimental.strategy = "miRNA-Seq")

GDCdownload(query_mirna, method = "api")
data <- GDCprepare(query_mirna)



############# CONSULTA DATOS RNA -> no quiere jalar

query_rna <- GDCquery(project = "TCGA-COAD",
               legacy = FALSE,
               data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification",
               workflow.type = "HTSeq - Counts",
               experimental.strategy = "RNA-Seq")

GDCdownload(query_rna, method = "api")
data_RNA <- GDCprepare(query_rna)













