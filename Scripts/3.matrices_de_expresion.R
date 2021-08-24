#### 3. MATRICES DE EXPRESIÓN MIRNAS Y RNA 

library(TCGAbiolinks)
library(dplyr)


# consulta de datos de expresión de miRNAs

query_mirna <- GDCquery(project = "TCGA-COAD",
               legacy = FALSE,
               data.category = "Transcriptome Profiling",
               data.type = "Gene Expression Quantification",
               workflow.type = "HTSeq - Counts")

GDCdownload(query_mirna, method = "api")
data <- GDCprepare(query_mirna)



