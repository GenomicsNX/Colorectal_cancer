if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("maftools")
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(maftools)


colnames(metaData)[1] <- "Tumor_Sample_Barcode"

?read.maf
?GDCquery_Maf
maf <- GDCquery_Maf("COAD", pipelines = "mutect2")  
maf2 <- GDCquery_Maf("READ", pipelines = "mutect2")


datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


which(maf@data@Hugo_Symbol %in% "EGFR")

maf <- read.maf(maf)
maf2 <- read.maf(maf2)
All = maftools:::merge_mafs(maf=c(maf,maf2), verbose = TRUE)

#####################################################################
test<-All


test@data$Tumor_Sample_Barcode<-substr(test@data$Tumor_Sample_Barcode,1,16)
test@clinical.data$Tumor_Sample_Barcode<-substr(test@clinical.data$Tumor_Sample_Barcode,1,16)

t<-BMImetaData[,c(2,4)]

colnames(t)[1]<-"Tumor_Sample_Barcode"

test_merge<-merge(test@data,t,by="Tumor_Sample_Barcode")

test@data<-test_merge

test_merge<-merge(test@clinical.data,t,by="Tumor_Sample_Barcode")

test@clinical.data<-test_merge

All_merge<-test

d<-cbind(All@clinical.data$Tumor_Sample_Barcode,All@clinical.data$Tumor_Sample_Barcode)
d<-as.data.frame(d)

d$V2<-substr(d$V2,1,16)

colnames(d)[2]<-"Tumor_Sample_Barcode"

d_merge<-merge(d,t,by="Tumor_Sample_Barcode")

d<-d_merge[-1]
clinical.data<-d[1]

colnames(d)[1]<- "Tumor_Sample_Barcode"

test_merge<-merge(test@data,d,by="Tumor_Sample_Barcode")

test@data<-test_merge

Tumor_Sample_Barcode<-d[1]

ddd<-merge(clinical.data,test@clinical.data,by="Tumor_Sample_Barcode")


test_merge<-merge(test@data,d,by="Tumor_Sample_Barcode")
test@data<-test_merge

t<-as.data.frame(clinical.data)

f<-merge(test@clinical.data,Tumor_Sample_Barcode,by="Tumor_Sample_Barcode")

test@clinical.data<-f
All_merge<-test

###################################################3

oncoplot(maf = All_merge, top = 20, showTumorSampleBarcodes = F, removeNonMutated = TRUE, fontSize = 0.5, draw_titv = FALSE,pathways = 'auto',
         clinicalFeatures ='Group')

titv = titv(maf = All_merge, plot = FALSE, useSyn = TRUE)

#plot titv summary
plotTiTv(res = titv)


lollipopPlot(
  maf = maf,
  gene = 'EGFR',
  AACol = 'HGVSp',
  showMutationRate = TRUE,
  labelPos = 858,
  legendTxtSize = 1,
  domainLabelSize = 0.8
  
)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

View(maf)





