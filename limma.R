# 1. limma packages
#### load package
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(limma)
library(ggplot2)
library(GEOquery)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)

query <- GDCquery(project = c("TCGA-READ","TCGA-COAD"),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Primary Tumor",
                  barcode = c( "TCGA-QG-A5Z2-01A-11R-A28H-07", "TCGA-G4-6298-01A-11R-1723-07", "TCGA-CM-6167-01A-11R-1653-07", "TCGA-G4-6627-01A-11R-1774-07", "TCGA-DM-A1D9-01A-11R-A155-07", "TCGA-CK-6746-01A-11R-1839-07",
                               "TCGA-A6-5661-01A-01R-1653-07", "TCGA-D5-5540-01A-01R-1653-07", "TCGA-F4-6703-01A-11R-1839-07", "TCGA-D5-6930-01A-11R-1928-07", "TCGA-AY-6196-01A-11R-1723-07", "TCGA-CM-5861-01A-01R-1653-07",
                               "TCGA-AY-5543-01A-01R-1653-07", "TCGA-AZ-4616-01A-21R-1839-07", "TCGA-G4-6294-01A-11R-1774-07", "TCGA-AD-6965-01A-11R-1928-07", "TCGA-CM-6675-01A-11R-1839-07", "TCGA-DM-A28G-01A-11R-A16W-07",
                               "TCGA-AA-3492-01A-01R-1410-07", "TCGA-A6-6782-01A-11R-1839-07", "TCGA-NH-A50U-01A-33R-A37K-07", "TCGA-A6-6649-01A-11R-1774-07", "TCGA-5M-AAT6-01A-11R-A41B-07", "TCGA-G4-6311-01A-11R-1723-07",
                               "TCGA-G4-6306-01A-11R-1774-07", "TCGA-D5-6538-01A-11R-1723-07", "TCGA-CK-6751-01A-11R-1839-07", "TCGA-CK-6747-01A-11R-1839-07", "TCGA-CM-5348-01A-21R-1723-07", "TCGA-CM-6162-01A-11R-1653-07",
                               "TCGA-CM-5860-01A-01R-1653-07", "TCGA-AZ-5407-01A-01R-1723-07", "TCGA-CM-6674-01A-11R-1839-07", "TCGA-CM-6168-01A-11R-1653-07", "TCGA-AZ-6599-01A-11R-1774-07", "TCGA-CM-5862-01A-01R-1653-07",
                               "TCGA-NH-A8F8-01A-72R-A41B-07", "TCGA-A6-4105-01A-02R-1774-07", "TCGA-G4-6320-01A-11R-1723-07", "TCGA-DM-A1DA-01A-11R-A155-07", "TCGA-D5-6531-01A-11R-1723-07", "TCGA-CM-6171-01A-11R-1653-07",
                               "TCGA-AD-6888-01A-11R-1928-07", "TCGA-DM-A1D8-01A-11R-A155-07", "TCGA-AD-A5EK-01A-11R-A28H-07", "TCGA-AA-3663-01A-01R-1723-07", "TCGA-AD-6901-01A-11R-1928-07", "TCGA-G4-6302-01A-11R-1723-07",
                               "TCGA-5M-AATE-01A-11R-A41B-07", "TCGA-AA-A01P-01A-21R-A083-07", "TCGA-AZ-6605-01A-11R-1839-07", "TCGA-G4-6586-01A-11R-1774-07", "TCGA-CA-5796-01A-01R-1653-07", "TCGA-DM-A1D4-01A-21R-A155-07",
                               "TCGA-D5-7000-01A-11R-A32Z-07", "TCGA-AA-3495-01A-01R-1410-07", "TCGA-QL-A97D-01A-12R-A41B-07", "TCGA-DM-A28K-01A-21R-A32Y-07", "TCGA-5M-AAT4-01A-11R-A41B-07", "TCGA-DM-A280-01A-12R-A16W-07",
                               "TCGA-AA-3675-01A-02R-0905-07", "TCGA-WS-AB45-01A-11R-A41B-07", "TCGA-A6-5665-01A-01R-1653-07", "TCGA-AZ-6598-01A-11R-1774-07", "TCGA-AD-A5EJ-01A-11R-A28H-07", "TCGA-AD-6895-01A-11R-1928-07",
                               "TCGA-3L-AA1B-01A-11R-A37K-07", "TCGA-A6-6138-01A-11R-1774-07", "TCGA-DM-A285-01A-11R-A16W-07", "TCGA-AY-A54L-01A-11R-A28H-07", "TCGA-AD-6964-01A-11R-1928-07", "TCGA-AD-5900-01A-11R-1653-07",
                               "TCGA-DM-A1HA-01A-11R-A155-07", "TCGA-AZ-4315-01A-01R-1410-07", "TCGA-CM-6680-01A-11R-1839-07", "TCGA-CM-6677-01A-11R-1839-07", "TCGA-DM-A0XD-01A-12R-A155-07", "TCGA-D5-6928-01A-11R-1928-07",
                               "TCGA-AU-6004-01A-11R-1723-07", "TCGA-AA-A02Y-01A-43R-A32Y-07", "TCGA-AD-6889-01A-11R-1928-07", "TCGA-CK-4950-01A-01R-1723-07", "TCGA-D5-5538-01A-01R-1653-07", "TCGA-G4-6626-01A-11R-1774-07",
                               "TCGA-NH-A50V-01A-11R-A28H-07", "TCGA-AY-6197-01A-11R-1723-07", "TCGA-AY-A71X-01A-12R-A37K-07", "TCGA-AZ-4614-01A-01R-1410-07", "TCGA-CA-6717-01A-11R-1839-07", "TCGA-G4-6295-01A-11R-1723-07",
                               "TCGA-CM-4751-01A-02R-1839-07", "TCGA-NH-A6GA-01A-11R-A37K-07", "TCGA-DM-A282-01A-12R-A16W-07", "TCGA-CM-4744-01A-01R-A32Z-07", "TCGA-D5-6540-01A-11R-1723-07", "TCGA-AZ-6600-01A-11R-1774-07",
                               "TCGA-4N-A93T-01A-11R-A37K-07", "TCGA-DM-A288-01A-11R-A16W-07", "TCGA-AA-A02K-01A-03R-A32Y-07", "TCGA-AA-A01Z-01A-11R-A083-07", "TCGA-AZ-4323-01A-21R-1839-07", "TCGA-CM-5349-01A-21R-1723-07",
                               "TCGA-CM-4743-01A-01R-1723-07", "TCGA-G4-6321-01A-11R-1723-07", "TCGA-QG-A5YW-01A-11R-A28H-07", "TCGA-CK-4951-01A-01R-1410-07", "TCGA-G4-6314-01A-11R-1723-07", "TCGA-AD-6890-01A-11R-1928-07",
                               "TCGA-DM-A28A-01A-21R-A32Y-07", "TCGA-G4-6588-01A-11R-1774-07", "TCGA-DM-A28H-01A-11R-A16W-07", "TCGA-CM-6166-01A-11R-1653-07", "TCGA-G4-6628-01A-11R-1839-07", "TCGA-G4-6297-01A-11R-1723-07",
                               "TCGA-CK-5912-01A-11R-1653-07", "TCGA-CA-6716-01A-11R-1839-07", "TCGA-DM-A0X9-01A-11R-A155-07", "TCGA-G4-6323-01A-11R-1723-07", "TCGA-AD-6963-01A-11R-1928-07", "TCGA-CK-4952-01A-01R-1723-07",
                               "TCGA-RU-A8FL-01A-11R-A37K-07", "TCGA-CM-5864-01A-01R-1653-07", "TCGA-CM-6169-01A-11R-1653-07", "TCGA-D5-5539-01A-01R-1653-07", "TCGA-AA-3496-01A-21R-1839-07", "TCGA-CM-5863-01A-21R-1839-07",
                               "TCGA-CK-5916-01A-11R-1653-07", "TCGA-A6-6137-01A-11R-1774-07", "TCGA-F4-6856-01A-11R-1928-07", "TCGA-CK-5913-01A-11R-1653-07", "TCGA-AZ-4615-01A-01R-1410-07", "TCGA-D5-6530-01A-11R-1723-07",
                               "TCGA-CM-4747-01A-01R-1410-07", "TCGA-AZ-6606-01A-11R-1839-07"))


GDCdownload(query)
COADREADdata <- GDCprepare(query)



datatable(as.data.frame(colData(COADREADdata)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
COADREAD_clinical<-colData(COADREADdata)
COADREAD_clinical<-as.data.frame(COADREAD_clinical)

#COADREAD_clinical<-COADREAD_clinical[-which(duplicated(COADREAD_clinical$patient)),]

COADREAD_clinical$barcode

COADREAD_clinical<- merge(right_dataset[,c(1310,1320)],COADREAD_clinical, by="barcode")


#Gene expression data

COADREAD_GE<-datatable(assay(COADREADdata), 
                       options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
                       rownames = TRUE)



#make a expression dataframe
COADREAD_GE<-as.data.frame(assay(COADREADdata))


##ENSEMBL ID to Gene Symbols 
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(clusterProfiler)
keytypes(org.Hs.eg.db)


eg = bitr(rownames(COADREAD_GE), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db") #39.39%

#expression dataframe rownames symbol로 변경 


COADREAD_GE$ENSEMBL<-rownames(COADREAD_GE)
COADREAD_GE_merge <- merge(COADREAD_GE,eg, by='ENSEMBL')
COADREAD_GE_merge<- COADREAD_GE_merge[,-1]
COADREAD_GE_merge<-aggregate(. ~ SYMBOL, data = COADREAD_GE_merge, mean) #56602->34429
rownames(COADREAD_GE_merge) <- COADREAD_GE_merge$SYMBOL
COADREAD_GE_merge<- COADREAD_GE_merge[,-1]




#metaData<-COADREAD_clinical
eData<-COADREAD_GE_merge

# 발현값 데이터에서 clinical 중복 환자만 가져오기

q<-BMImetaData[-1]
w<-merge(metaData,q, by='sample')

metaData<-w

# 발현값데이터와 메타데이터의 바코드 순서를 일치시켜줘야한다.

?order
sort(x, decreasing = F) 
eData<-eData[, order(colnames(eData))]
COADREAD_clinical<-COADREAD_clinical[order(COADREAD_clinical$barcode),]

colnames(eData)
COADREAD_clinical$barcode


#eData, metaData
#limma



BiocManager::install("Glimma")
BiocManager::install("edgeR")
BiocManager::install("affy")
library(limma)
library(Glimma)
library(edgeR)
library(affy)


metaData<-COADREAD_clinical
#design matrix 
?model.matrix
design <- model.matrix(~0+metaData$right_Macrophages.M1)
head(design)

colnames(design) <- c('High','Low')
head(design)

##############
#Normalization
##############
x<-eData

dim(x)

keep.exprs <- filterByExpr(x,design,min.count = 10, min.total.count = 15, min.prop = 0.7)

x <- x[keep.exprs,]

dim(x)
################# low count 제거 다른 방법
#nrow(x)
#keep <- rowSums(x)>50
#x <- x[keep,]
#dim(x)
############################## 
#Normalization 중에 1가지 선택

v <- voom(x, design, plot=F)

limmaVoom<-(v$E)


###################################
lcpm <- cpm(x, log=TRUE)

TMM <- cpm(x, normalized.lib.sizes=TRUE,log=T)


nrow(v$E) #17164
###################################

# 정규분포 확인
hist(as.matrix((v)))
hist(as.matrix((lcpm)))
hist(as.matrix((TMM)))

#boxplot

boxplot(lcpm,col=rgb(0.8,0.8,0.3,0.5),las=2) 


#############################################################################
#########	MDS plot 		###############################################
#############################################################################

plotMDS(lcpm, labels=metaData$right_Macrophages.M1, top = 10000)

######################################################
vfit <- lmFit(v, design)
cont <- makeContrasts(diff=High-Low,levels=design)
vfit.cont <- contrasts.fit(vfit,cont)
vfit.cont <- eBayes(vfit.cont)
res <- topTable(vfit.cont,number=Inf,adjust="BH")
?topTable
plotSA(vfit.cont, main="Final model: Mean-variance trend")

#####################################################
View(v$E)
res

up_regulated<-res[res$adj.P.Val < 0.05 & res$logFC > 1,]
down_regulated<-res[res$adj.P.Val < 0.05 & res$logFC < -1,]

hist(res$logFC)

rownames(up_regulated)
rownames(down_regulated)


nrow(up_regulated) #434
nrow(down_regulated) #33
rownames(up_regulated)
rownames(down_regulated)

write.csv(up_regulated,file = "/home/dongmin/colorectal cancer microbiome/adddata/limma_up_regulated.csv")
write.csv(down_regulated,file = "/home/dongmin/colorectal cancer microbiome/adddata/limma_down_regulated.csv")


#GO KEGG

#Input genes; convert to ENTREZID 
eg = bitr(rownames(up_regulated), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg1 = bitr(rownames(down_regulated), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#GO-Enrichment analysis 

go <- enrichGO(eg$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") #keyType ='SYMBOL',pAdjustMethod = "BH",
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") # ont= "BP","CC","MF"

dotplot(go, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=7) + facet_grid(ONTOLOGY~., scale="free")
dotplot(go1, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=6) + facet_grid(ONTOLOGY~., scale="free")
##############################
go <- enrichGO(eg$SYMBOL,keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") #keyType ='SYMBOL',pAdjustMethod = "BH",
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1) # ont= "BP","CC","MF"

dotplot(go, split="ONTOLOGY",color="pvalue", font.size =10,showCategory=7) + facet_grid(ONTOLOGY~., scale="free")

View(go@result)

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = "",
                x = 'logFC',
                y = 'P.Value',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 5.0,
                #colAlpha = 1,
                #legend=c('NS','Log (base 2) fold-change','P value',
                #         'P value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)





kg <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)

kg1 <- enrichKEGG(gene         = eg1$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 1)

dotplot(kg,color='pvalue',showCategory=15)
dotplot(kg1,color='pvalue',showCategory=10)

View(kg@result)

ggsave(filename ='/home/dongmin/colorectal cancer microbiome/adddata/GO_up_regulated.png',
       device = 'png',
       dpi = 1200
)


dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
plotMA(res)


gseaplot(kg,geneSetID = 1)
