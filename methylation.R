library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)

#### GDC 쿼리 설정 

#project = TCGA 프로젝트명 EX) TCGA-BRCA, TCGA-COAD, TCGA-READ 등등 
#data.category = simple nucleotide variation, copy number variation, transcriptome profiling 등등 
#data.type = Gene Expression Quantification, Isoform Expression Quantification, miRNA Expression Quantification 등등 
#그냥 모를때 ?GDCquery 치면 help탭에 설명 나와있음  

query.met <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450")
)

query.exp <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "HTSeq - Counts"
)

GDCdownload(query.exp)
data <- GDCprepare(query.exp)

# Methylation데이터랑 발현값 데이터를 모두 가지고 있는 환자 추출하기 ㅎ
common.patients <- intersect(
  substr(getResults(query.met, cols = "cases"), 1, 12),
  substr(getResults(query.exp, cols = "cases"), 1, 12)
)

####다운받은 TCGA데이터에서 발현값 뽑아오기 
DF_colorectal <- assay(data)
####앙상블 및 유전자 심볼 찾기  
data@rowRanges[,c(1,2)]
DF_coldata<-colData(data)
DF_coldata <- as.data.frame(DF_coldata)
View(DF_coldata)
elementMetadata(data)
metadata(data)

View(DF_coldata)
df<-assay(data)


####################임상정보 
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
clinical.BCRtab.all$clinical_patient_coad ##TNM이나 기초정보 

clinical.BCRtab.all$clinical_patient_coad


View(clinical.BCRtab.all$clinical_patient_coad)


COAD_BMI<-clinical.BCRtab.all$clinical_patient_coad[,c(2,47,48)]

q<-which(COAD_BMI$weight_kg_at_diagnosis %in% "[Not Available]")
w<-which(COAD_BMI$height_cm_at_diagnosis %in% "[Not Available]")

COAD_BMI<-COAD_BMI[-c(q,w),]

COAD_BMI<-COAD_BMI[-c(1,2),]

exp_bmi<-exp_bmi[,which(!duplicated(colnames(exp_bmi)))]

length(df_bmi$bcr_patient_barcode)

COAD_BMI$weight_kg_at_diagnosis <- as.numeric(COAD_BMI$weight_kg_at_diagnosis)
COAD_BMI$height_cm_at_diagnosis<-as.numeric(COAD_BMI$height_cm_at_diagnosis)

COAD_BMI$height_cm_at_diagnosis<-round(COAD_BMI$height_cm_at_diagnosis*0.01,2)

COAD_BMI$weight_kg_at_diagnosis[1]/(COAD_BMI$height_cm_at_diagnosis[1]^2)
COAD_BMI$BMI <-COAD_BMI$weight_kg_at_diagnosis[i]/(COAD_BMI$height_cm_at_diagnosis[i]^2)

for(i in 1:233){
  COAD_BMI$BMI[i] <-COAD_BMI$weight_kg_at_diagnosis[i]/(COAD_BMI$height_cm_at_diagnosis[i]^2)
  
}

z<-which(COAD_BMI$BMI <= 25)
x<-which(COAD_BMI$BMI > 30)
low<-COAD_BMI[z,]
high<-COAD_BMI[x,]

high$group <- "high"
low$group <- "low"

length(high$group)
length(low$group)

df_bmi<-rbind(low,high)

a<-which(colnames(df) %in% high$bcr_patient_barcode)
b<-which(colnames(df) %in% low$bcr_patient_barcode)

unique(colnames((df[,a])))
ncol(df[,b])

ncol(df[,b])
ncol(df[,a])

exp_bmi<-cbind(df[,b],df[,a])


ncol(cbind(df[,b],df[,a]))

colnames(exp_bmi)[1:10]
df_bmi$bcr_patient_barcode[1:10]

df_bmi<-df_bmi[order(df_bmi$bcr_patient_barcode),]

colnames(exp_bmi)[1:10]
df_bmi$bcr_patient_barcode[1:10]


exp_bmi<-exp_bmi[,order(colnames(exp_bmi))]

exp_bmi ### expression_data
df_bmi ### meta_data

#################XML 데이터 
############다운 받은 clinical data에서 원하는 정보만 분리 by clinical.info
############clinical.info에 포함되어있는거(drug, admin, follow_up, radiation, patient, stage_event, new_tumor_event, 
############sample, bio_patient, analyte, aliquot, protocol, portion, slide, msi)
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)

clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")

######### RNA-seq analysis를 위한 Dataframe 형성 
######### Lymph-node status에 따른 데이터프레임 형성 


colnames(df)<-substr(colnames(df),1,12)

group <- data@colData$ajcc_pathologic_n
group <- ifelse(group=="N0","N0","N1")

#####lymph node NA 값 찾기 
which(data@colData$ajcc_pathologic_n %in% NA)
#####lymph node NA 값 제거 
group <- group[-c(131,381)]
df <- df[,-c(131,381)]


exp_bmi ### expression_data
df_bmi ### meta_data

x<-DGEList(exp_bmi)

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

table(rowSums(x$counts==0)==9)

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors


x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

design <- model.matrix(~0+group)
colnames(design) <- c("high","low")
design

contr.matrix <- makeContrasts(high-low,levels = colnames(design))
contr.matrix


vfit1 <- lmFit(x$counts, design)
vfit1 <- contrasts.fit(vfit1, contrasts=contr.matrix)
efit1 <- eBayes(vfit1)
topTable(efit1)

##### data pre-processing 
v <- voom(exp_bmi, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

top.table <- topTable(efit, sort.by = "P", n = Inf)

View(top.table)

hist(top.table$logFC)


plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))

tfit <- treat(vfit, lfc=0)
dt <- decideTests(tfit)
summary(dt)



#######Volcanoplot 

install.packages('ggalt')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
EnhancedVolcano(res_43292,
                lab = rownames(res_43292),
                x = 'logFC',
                y = 'P.Value',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 2.0,
                labSize = 5.0,
                #colAlpha = 1,
                #legend=c('NS','Log (base 2) fold-change','P value',
                #         'P value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)

#####bmi calculation 












