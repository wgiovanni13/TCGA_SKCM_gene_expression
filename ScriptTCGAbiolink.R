library(TCGAbiolinks)
library(glmnet)
library(factoextra)
library(FactoMineR)
library(caret)
library(SummarizedExperiment)
library(ggplot2)
library(RColorBrewer)
library(gProfileR)
library(genefilter)
library(GenomicDataCommons)
library(GenomeInfoDbData)
library(keras)
library(tensorflow)
library(dplyr)
library(biomaRt)
library(annotables)
library(org.Hs.eg.db)
library(tidyverse)
library(plyr)
library(DBI)
library(IOBR)

GDCprojects <- getGDCprojects()
TCGAbiolinks:::getProjectSummary("TCGA-SKCM")

query_TCGA <- GDCquery(project = "TCGA-SKCM", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq", barcode = c("TCGA-*"))

query <- GDCquery(project = "TCGA-SKCM", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq", barcode = c("TCGA-*"))

GDCdownload(query = query_TCGA, method = "api", files.per.chunk = 100)
TCGA_data <- GDCprepare(query_TCGA)
TCGAMatrix <- assay(TCGA_data,"unstranded") 

TCGA.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(TCGA_data)

SKCM.Rnaseq <- GDCprepare(query)

datanorm <- TCGAanalyze_Normalization(
  tabDF = TCGA.RNAseq_CorOutliers, 
  geneInfo =  geneInfoHT
)

dataFilt <- TCGAanalyze_Filtering(
  tabDF = datanorm,
  method = "quantile", 
  qnt.cut =  0.25
)

samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

ensembl <- useEnsembl(biomart = "genes")
dataset <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

ensemblname <- as.vector(ensemblname)
genes <- data.frame(rownames(dataFilt))

GeneIDsSkin <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                 filters = "ensembl_gene_id", 
                 values = genes, 
                 mart = ensembl.con)

dataFilt <- data.frame(Ensemblname = row.names(dataFilt), dataFilt)

dataFilt2<- merge(dataFilt,GeneIDsSkin,by.x="Ensemblname",by.y="ensembl_gene_id")

dataFilt2<-dataFilt2[,c(475, 1:474)]

dataFilt2<-dataFilt2 %>%
  # recode empty strings "" by NAs
  na_if("") %>%
  # remove NAs
  na.omit

dataFilt2<-dataFilt2 %>% distinct()

SKCMMatrixcibersort <- dataFilt2[,c(1, 3:475)]
SKCMMatrixcibersort  <- SKCMMatrixcibersort[!duplicated(SKCMMatrixcibersort[,c('external_gene_name')]),]
rownames(SKCMMatrixcibersort) <- SKCMMatrixcibersort$external_gene_name
SKCMMatrixcibersort <- SKCMMatrixcibersort[,-c(1)]
colnames(SKCMMatrixcibersort) <- gsub(x = colnames(SKCMMatrixcibersort), pattern = "\\.", replacement = "-")  
n <- 3
colnames(SKCMMatrixcibersort) <- vapply(strsplit(colnames(SKCMMatrixcibersort), '-'), function(x) paste(x[seq.int(n)], collapse='-'), character(1L))

#Clinnical data

clinicSKCM <- GDCquery_clinic("TCGA-SKCM", "clinical")
clinicSKCMsurvive <- clinicSKCM[,-c(1,3:37, 39:69)]
clinicSKCMgender <- clinicSKCM[,-c(1,3:35, 37:69)]
clinicSKCMsurvive$vital_status <- revalue(clinicSKCMsurvive$vital_status, c("Dead" = "0", "Alive" = "1"))
clinicSKCMgender$gender <- revalue(clinicSKCMgender$gender, c("male" = "0", "female" = "1"))


SKCMMatrixcibersorttranspuesta <- t(SKCMMatrixcibersort)
SKCMMatrixcibersorttranspuesta <- sqlRownamesToColumn(SKCMMatrixcibersorttranspuesta)
SKCMMatrixcibersorttranspuesta <- data.frame(t(SKCMMatrixcibersort))
SKCMMatrixcibersorttranspuesta$IDCase <- rownames(SKCMMatrixcibersorttranspuesta)
SKCMMatrixcibersorttranspuesta$IDCase <- gsub(x = SKCMMatrixcibersorttranspuesta$IDCase, pattern = "\\.", replacement = "-")

SKCMMatrixcibersorttranspuesta <- SKCMMatrixcibersorttranspuesta[,c(which(colnames(SKCMMatrixcibersorttranspuesta)=="IDCase"),which(colnames(SKCMMatrixcibersorttranspuesta)!="IDCase"))]
SKCMMatrixcibersorttranspuesta$IDCase <- gsub(x = SKCMMatrixcibersorttranspuesta$IDCase, pattern = "\\.", replacement = "-")
SKCMMatrixcibersorttranspuesta <- merge(clinicSKCMsurvive, SKCMMatrixcibersorttranspuesta, by.x = "submitter_id", by.y = "IDCase")
SKCMMatrixcibersorttranspuesta <- SKCMMatrixcibersorttranspuesta[,-c(2:3)]
rownames(SKCMMatrixcibersorttranspuesta) <- SKCMMatrixcibersorttranspuesta$submitter_id
SKCMMatrixcibersorttranspuesta <- SKCMMatrixcibersorttranspuesta[-c(1),]
SKCMMatrixcibersorttranspuesta <- t(SKCMMatrixcibersorttranspuesta)
cibersort469casos<-deconvo_tme(eset = SKCMMatrixcibersorttranspuesta, method = "cibersort", arrays = FALSE, perm = 200)

clinicSKCMsurvive <- clinicSKCMsurvive[-c(368),]
inputRF2 <- merge(cibersort469casos,clinicSKCMsurvive, by.x = "ID", by.y = "submitter_id")
inputRF2 <- inputRF2[, -c(24:26)]

write.csv(inputRF2, "\\InputRFcibersort.csv", row.names=FALSE, quote = FALSE)

xcell469casos<-deconvo_tme(eset = SKCMMatrixcibersorttranspuesta, method = "xcell", arrays = FALSE)
xcell469casos <- xcell469casos[, -c(65:68)]
clinicSKCM <- clinicSKCM[-c(325),]
inputRFxcellcibersort <- merge(xcell469casos,inputRF2, by.x = "ID", by.y = "ID")
inputRFxcellcibersort <- inputRFxcellcibersort[, -c(1)]

write.csv(inputRFxcellcibersort, "\\inputRFxcellcibersort.csv", row.names=FALSE, quote = FALSE)

inputRFxcellcibersortgender <-  merge(cibersort469casos,clinicSKCMgender, by.x = "ID", by.y = "submitter_id")
inputRFxcellcibersortgender2 <- merge(xcell469casos,inputRFxcellcibersortgender, by.x = "ID", by.y = "ID")
inputRFxcellcibersortgender2 <- inputRFxcellcibersortgender2[, -c(87:89)]

write.csv(inputRFxcellcibersortgender2, "\\inputRFxcellcibersortgender.csv", row.names=FALSE, quote = FALSE)

inputRFxcellcibersortgendersurvive <-  merge(inputRFxcellcibersort,clinicSKCMgender, by.x = "ID", by.y = "submitter_id")

write.csv(inputRFxcellcibersortgendersurvive, "\\inputRFxcellcibersortgendersurvive.csv", row.names=FALSE, quote = FALSE)

