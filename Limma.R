
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
                  sample.type = "Primary Tumor")


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
################# 
#nrow(x)
#keep <- rowSums(x)>50
#x <- x[keep,]
#dim(x)
##################


v <- voom(x, design, plot=F)

limmaVoom<-(v$E)


###################################
lcpm <- cpm(x, log=TRUE)

TMM <- cpm(x, normalized.lib.sizes=TRUE,log=T)


nrow(v$E) #17164
###################################


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

ggsave(filename ='/adddata/GO_up_regulated.png',
       device = 'png',
       dpi = 1200
)