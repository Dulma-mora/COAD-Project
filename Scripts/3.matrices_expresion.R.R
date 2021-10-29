
library("TCGAbiolinks")
library("SummarizedExperiment")
query_mirna <- GDCquery(project = "TCGA-COAD",
                        legacy = FALSE,
                        data.category = "Transcriptome Profiling", 
                        data.type = "Isoform Expression Quantification",
                        #workflow.type = "HTSeq - Counts", 
                        experimental.strategy = "miRNA-Seq")

GDCdownload(query_mirna, method = "api")
data <- GDCprepare(query_mirna)

write.csv(data, "data_mirnas.csv")

############# CONSULTA DATOS RNA -> no quiere jalar

query_rna <- GDCquery(project = "TCGA-COAD",
                      legacy = FALSE,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts",
                      experimental.strategy = "RNA-Seq")

GDCdownload(query_rna, method = "api")
data_RNA <- GDCprepare(query_rna)
data_rna_exp <- assay(data_RNA)
write.csv(data_rna_exp, "data_mrna_transcriptome.csv")
