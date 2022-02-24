#TCGA-FPKM data download

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("limma")
install.packages('DT')
BiocManager::install("GEOquery")

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(limma)
library(ggplot2)
library(GEOquery)


#Search and download data 

query <- GDCquery(project = "TCGA-READ",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM",
                  sample.type = "Primary Tumor")


GDCdownload(query)
READdata <- GDCprepare(query)

query2 <- GDCquery(project = "TCGA-COAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "HTSeq - FPKM",
                   sample.type = "Primary Tumor")

GDCdownload(query2)
COADdata <- GDCprepare(query2)

datatable(as.data.frame(colData(READdata)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
READ_clinical<-colData(READdata)
READ_clinical<-as.data.frame(READ_clinical)

datatable(as.data.frame(colData(COADdata)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

COAD_clinical<-colData(COADdata)
COAD_clinical<-as.data.frame(COAD_clinical)

COADREAD_clinical<- rbind(COAD_clinical,READ_clinical)

rownames(READdata)

#Gene expression data

READ_GE<-datatable(assay(READdata), 
                   options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
                   rownames = TRUE)


COAD_GE<-datatable(assay(COADdata), 
                   options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
                   rownames = TRUE)


#make a expression dataframe
READ_GE<-as.data.frame(assay(READdata))
COAD_GE<-as.data.frame(assay(COADdata))


##ENSEMBL ID to Gene Symbols 
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
install.packages("maftools")
library(org.Hs.eg.db)
library(clusterProfiler)
keytypes(org.Hs.eg.db)


eg = bitr(rownames(READ_GE), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db") #39.39%
eg2 = bitr(rownames(COAD_GE), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")#39.39%

#expression dataframe rownames symbol로 변경 


READ_GE$ENSEMBL<-rownames(READ_GE)
READ_GE_merge <- merge(READ_GE,eg, by='ENSEMBL')
READ_GE_merge<- READ_GE_merge[,-1]
READ_GE_merge<-aggregate(. ~ SYMBOL, data = READ_GE_merge, mean)
rownames(READ_GE_merge) <- READ_GE_merge$SYMBOL
READ_GE_merge<- READ_GE_merge[,-1]

COAD_GE$ENSEMBL<-rownames(COAD_GE)
COAD_GE_merge <- merge(COAD_GE,eg, by='ENSEMBL')
COAD_GE_merge<- COAD_GE_merge[,-1]
COAD_GE_merge<-aggregate(. ~ SYMBOL, data = COAD_GE_merge, mean)
rownames(COAD_GE_merge) <- COAD_GE_merge$SYMBOL
COAD_GE_merge<- COAD_GE_merge[,-1]

View(COAD_GE_merge)

expression_merge_backup<-expression_merge
expression_merge<-cbind(COAD_GE_merge,READ_GE_merge)

colnames(expression_merge)<-substr(colnames(expression_merge),1,12)


w(colnames(expression_merge))

expression_merge<-expression_merge[,-which(duplicated(colnames(expression_merge)))]



# BMI high , low 데이터 구분

clinical_High <- clinical[clinical$Group == 'High',]
clinical_Low <- clinical[clinical$Group == 'Low',]


q<-which(colnames(expression_merge) %in% clinical_High$patient)
w<-which(colnames(expression_merge) %in% clinical_Low$patient)

final_High_expression <- expression_merge[,q]
final_Low_expression <- expression_merge[,w]

colnames(final_High_expression)
str(final_High_expression)

which(duplicated(colnames(final_Low_expression)))

# 발현값 데이터 High=88,Low=90 명

save.image(file = "backup.RData")
# quantiseqr 
# final_left_FPKM 133명 , final_right_FPKM 134명

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("quantiseqr")

library("quantiseqr")
install.packages("tidyr")

library("dplyr")
library("ggplot2")
library("tidyr")
library("tibble")
library("GEOquery")
library("reshape2")
library("SummarizedExperiment")


#1.quantiplot



ti_racle_High <- quantiseqr::run_quantiseq(
  expression_data = limmaVoom_High,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = FALSE
)


quantiplot(ti_racle_High)




ti_racle_Low <- quantiseqr::run_quantiseq(
  expression_data = limmaVoom_Low,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = FALSE
)



quantiplot(ti_racle_Low)

ti_racle_High<-as.data.frame(ti_racle_High)
ti_racle_Low<-as.data.frame(ti_racle_Low)

write.csv(final_High_expression,file = "final_High_expression.csv")
write.csv(final_Low_expression,file = "final_Low_expression.csv")


##################################
colnames(ti_racle_High)
#"B.cells"        "Macrophages.M1"  "Macrophages.M2"  "Monocytes"       "Neutrophils"     "NK.cells"        "T.cells.CD4"     "T.cells.CD8"    
#"Tregs"           "Dendritic.cells" 




############################################### "two.sided", "less", "greater"
#윌콕스 양측검정
wilcox.test(ti_racle_High$Dendritic.cells,ti_racle_Low$Dendritic.cells,alternative="two.sided")

#윌콕스 단측검정. 대립가설 왼쪽 < 오늘쪽 (평균값)
wilcox.test(ti_racle_High$T.cells.CD4,ti_racle_Low$T.cells.CD4,alternative="less")

#윌콕스 단측검정. 대립가설 왼쪽 > 오늘쪽 (평균값)
wilcox.test(ti_racle_High$B.cells,ti_racle_Low$B.cells,alternative="greater")
##############################################




#Boxplot 

library(ggplot2)



p <- ggplot(data = boxplot, aes(x=Location, y=Fraction)) + 
  geom_boxplot(aes(fill=Location),width=0.7)
p + facet_wrap( ~ Cells, scales="free")+theme(axis.text.x = element_text(size=13,
                                                                         face='bold'
),axis.text.y = element_text(size =11,face = 'plain'),legend.text=element_text(size=13,face='bold'),
legend.title=element_text(size=13,face='bold'),strip.text.x = element_text(size = 13),axis.title = element_text(size = 20) )+stat_compare_means(method = "wilcox.test")



