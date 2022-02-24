library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)


query_met.hg38 <- GDCquery(
  project= c("TCGA-COAD","TCGA-READ"), 
  data.category = "DNA Methylation", 
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")

BiocManager::install("sesame")
library("sesame")
sesameDataCacheAll()
library("ChAMP")

myLoad <- data.hg38

cham.load(myLoad)

# Using sesame  http://bioconductor.org/packages/sesame/
# Please cite 10.1093/nar/gky691 and doi: 10.1093/nar/gkt090.
library(TCGAbiolinks)
proj <- c("TCGA-COAD","TCGA-READ")
query <- GDCquery(
  project = c("TCGA-COAD","TCGA-READ"),
  data.category = "Raw microarray data",
  data.type = "Raw intensities", 
  experimental.strategy = "Methylation array", 
  legacy = TRUE,
  barcode = c("TCGA-3L-AA1B-01A", "TCGA-4N-A93T-01A", "TCGA-4T-AA8H-01A", "TCGA-5M-AAT6-01A", "TCGA-5M-AATE-01A", "TCGA-A6-2671-01A", "TCGA-A6-2672-01A", "TCGA-A6-2676-01A", "TCGA-A6-2677-01B", "TCGA-A6-2678-01A", "TCGA-A6-2679-01A",
              "TCGA-A6-2680-01A", "TCGA-A6-2683-01A", "TCGA-A6-3807-01A", "TCGA-A6-3809-01A", "TCGA-A6-3810-01B", "TCGA-A6-4107-01A", "TCGA-A6-5656-01A", "TCGA-A6-5657-01A", "TCGA-A6-5659-01A", "TCGA-A6-5660-01A", "TCGA-A6-5661-01B",
              "TCGA-A6-5664-01A", "TCGA-A6-5665-01B", "TCGA-A6-5667-01A", "TCGA-A6-6140-01A", "TCGA-A6-6141-01A", "TCGA-A6-6142-01A", "TCGA-A6-6650-01B", "TCGA-A6-6651-01A", "TCGA-A6-6652-01A", "TCGA-A6-6654-01A", "TCGA-A6-6780-01B",
              "TCGA-A6-6781-01B", "TCGA-A6-6782-01A", "TCGA-A6-A567-01A", "TCGA-A6-A56B-01A", "TCGA-A6-A5ZU-01A", "TCGA-AD-6548-01A", "TCGA-AD-6895-01A", "TCGA-AD-6899-01A", "TCGA-AD-6901-01A", "TCGA-AF-2690-01A", "TCGA-AF-2691-01A",
              "TCGA-AF-2692-01A", "TCGA-AF-2693-01A", "TCGA-AF-5654-01A", "TCGA-AF-6655-01A", "TCGA-AF-6672-01A", "TCGA-AF-A56K-01A", "TCGA-AH-6903-01A", "TCGA-AM-5820-01A", "TCGA-AM-5821-01A", "TCGA-AU-3779-01A", "TCGA-AU-6004-01A",
              "TCGA-AY-4070-01A", "TCGA-AY-6386-01A", "TCGA-AY-A54L-01A", "TCGA-AY-A69D-01A", "TCGA-AY-A71X-01A", "TCGA-AZ-6601-01A", "TCGA-CM-4743-01A", "TCGA-CM-4746-01A", "TCGA-CM-4747-01A", "TCGA-CM-4751-01A", "TCGA-CM-5341-01A",
              "TCGA-CM-5348-01A", "TCGA-CM-5349-01A", "TCGA-CM-5860-01A", "TCGA-CM-5861-01A", "TCGA-CM-5863-01A", "TCGA-CM-5864-01A", "TCGA-CM-5868-01A", "TCGA-CM-6161-01A", "TCGA-CM-6162-01A", "TCGA-CM-6163-01A", "TCGA-CM-6166-01A",
              "TCGA-CM-6168-01A", "TCGA-CM-6169-01A", "TCGA-CM-6170-01A", "TCGA-CM-6172-01A", "TCGA-CM-6676-01A", "TCGA-CM-6677-01A", "TCGA-CM-6678-01A", "TCGA-CM-6679-01A", "TCGA-CM-6680-01A", "TCGA-D5-6529-01A", "TCGA-D5-6531-01A",
              "TCGA-D5-6532-01A", "TCGA-D5-6533-01A",
              "TCGA-D5-6536-01A", "TCGA-D5-6538-01A", "TCGA-D5-6539-01A", "TCGA-D5-6540-01A", "TCGA-D5-6931-01A", "TCGA-DC-4745-01A", "TCGA-DC-5869-01A", "TCGA-DC-6154-01A", "TCGA-DC-6155-01A", "TCGA-DC-6157-01A", "TCGA-DC-6160-01A", "TCGA-DC-6682-01A", "TCGA-DC-6683-01A", "TCGA-DM-A0XD-01A",
              "TCGA-DM-A1D0-01A", "TCGA-DM-A1D4-01A", "TCGA-DM-A1D7-01A", "TCGA-DM-A1D9-01A", "TCGA-DM-A1DA-01A", "TCGA-DM-A1HB-01A", "TCGA-DM-A280-01A", "TCGA-DM-A282-01A", "TCGA-DM-A285-01A", "TCGA-DM-A28A-01A", "TCGA-DM-A28C-01A",
              "TCGA-DM-A28E-01A", "TCGA-DM-A28G-01A", "TCGA-DM-A28M-01A", "TCGA-DY-A1DD-01A", "TCGA-DY-A1DF-01A", "TCGA-DY-A1H8-01A", "TCGA-EF-5831-01A", "TCGA-EI-6506-01A", "TCGA-EI-6507-01A", "TCGA-EI-6512-01A", "TCGA-EI-6882-01A",
              "TCGA-EI-6883-01A", "TCGA-EI-6885-01A", "TCGA-F4-6459-01A", "TCGA-F4-6461-01A", "TCGA-F4-6463-01A", "TCGA-F4-6569-01A", "TCGA-F4-6570-01A", "TCGA-F4-6704-01A", "TCGA-F4-6805-01A", "TCGA-F4-6806-01A", "TCGA-F4-6807-01A",
              "TCGA-F4-6809-01A", "TCGA-F4-6856-01A", "TCGA-F5-6465-01A", "TCGA-F5-6811-01A", "TCGA-F5-6814-01A", "TCGA-F5-6861-01A", "TCGA-F5-6864-01A", "TCGA-G4-6299-01A", "TCGA-G4-6302-01A", "TCGA-G4-6306-01A",
              "TCGA-G4-6309-01A", "TCGA-G4-6310-01A", "TCGA-G4-6311-01A", "TCGA-G4-6314-01A", "TCGA-G4-6315-01A", "TCGA-G4-6317-01A", "TCGA-G4-6320-01A", "TCGA-G4-6321-01A", "TCGA-G4-6323-01A", "TCGA-G4-6588-01A", "TCGA-G4-6628-01A",
              "TCGA-G5-6233-01A", "TCGA-G5-6572-01A", "TCGA-G5-6641-01A", "TCGA-NH-A50V-01A", "TCGA-NH-A5IV-01A", "TCGA-NH-A6GA-01A", "TCGA-NH-A6GC-01A", "TCGA-NH-A8F8-01A", "TCGA-QG-A5YV-01A", "TCGA-RU-A8FL-01A", "TCGA-SS-A7HO-01A",
              "TCGA-T9-A92H-01A", "TCGA-WS-AB45-01A"),
  file.type = ".idat",
  platform = "Illumina Human Methylation 450")

tryCatch(
  GDCdownload(query, method = "api", files.per.chunk = 20),
  error = function(e) GDCdownload(query, method = "client")
)

betas <- GDCprepare(query)

match.file.cases.all <- NULL
match.file.cases <- getResults(query,cols=c("cases","file_name"))
match.file.cases$project <- proj
match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)

write.csv(match.file.cases,file = "match.file.cases.csv")
cases<-read.csv(file = "IDAT/cases.csv")
cases_backup<-cases


cases_backup$sample<-substr(cases_backup$sample,1,16)
cases_backup$Sentrix_Position<-substr(cases_backup$Sentrix_Position,12,17)
cases_backup$Sentrix_ID<-substr(cases_backup$Sentrix_ID,1,10)


cases<-cases_backup

w<-BMImetaData[c(2,4)]

cases<- merge(cases,w,by='sample')

write.csv(cases,file = "cases.csv")
cases_new<-read.csv(file="./IDAT/cases.csv")
cases_new_backup<-cases_new
cases_new<-cases_new[-1]


readr::write_tsv(cases_new, path =  "./IDAT/idat_filename_case.txt")


############################ CHAMP

?champ.load
getwd()
myLoad <- champ.load(directory = "/home/dongmin/BMI methylation/IDAT",
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=1,
                     filterBeads=TRUE,
                     beadCutoff=1,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="450K")
MYLOAD<-myLoad

?champ.import
myImport  <- champ.import(directory = "/home/dongmin/BMI methylation/IDAT",
                         offset=100,
                         arraytype="450K")

?champ.filter
myLoad <- champ.filter(beta=myImport$beta,
                         M=NULL,
                         pd=myImport$pd,
                         intensity=NULL,
                         Meth=NULL,
                         UnMeth=NULL,
                         detP=NULL,
                         beadcount=NULL,
                         autoimpute=FALSE,
                         filterDetP=TRUE,
                         ProbeCutoff=0,
                         SampleCutoff=0.1,
                         detPcut=1,
                         filterBeads=TRUE,
                         beadCutoff=1,
                         filterNoCG = TRUE,
                         filterSNPs = TRUE,
                         population = NULL,
                         filterMultiHit = TRUE,
                         filterXY = TRUE,
                         fixOutlier = TRUE,
                         arraytype = "450K")


CpG.GUI(CpG=rownames(myLoad$beta),arraytype="450K")
?QC.GUI
#QC.GUI(beta=myLoad$beta,pheno=myLoad$pd$Sample_Group,arraytype="450K")
#Normalization
?champ.norm
myNorm <-  champ.norm(beta=myLoad$beta,
                      rgSet=myLoad$rgSet,
                      mset=myLoad$mset,
                      resultsDir="./CHAMP_Normalization/",
                      method="BMIQ",
                      plotBMIQ=FALSE,
                      arraytype="450K",
                      cores=10)

champ.SVD(beta=myNorm,pd=myLoad$pd)

hist(as.matrix(myNorm))

# Batch Effect Correction


?champ.runCombat
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"),logitTrans=FALSE)



#########

myDMP <- champ.DMP(beta = myCombat,
                   pheno = myLoad$pd$Sample_Group,
                   compare.group = c("High","Low"),
                   adjPVal = 1,
                   adjust.method = "BH",
                   arraytype = "450K")

myDMP$High_to_Low

?DMP.GUI
DMP.GUI(DMP=myDMP[[1]],
        beta=myCombat,
        pheno=myLoad$pd$Sample_Group,
        cutgroupnumber=2)

cg10340801

myDMR <- champ.DMR(beta=myCombat,
                             pheno=myLoad$pd$Sample_Group,
                             compare.group=NULL,
                             arraytype="450K",
                             method = "Bumphunter",
                             minProbes=7,
                             adjPvalDmr=1,
                             cores=3,
                             ## following parameters are specifically for Bumphunter method.
                             maxGap=300,
                             cutoff=NULL,
                             pickCutoff=TRUE,
                             smooth=TRUE,
                             smoothFunction=loessByCluster,
                             useWeights=FALSE,
                             permutations=NULL,
                             B=250,
                             nullMethod="bootstrap",
                             ## following parameters are specifically for probe ProbeLasso       method.
                             meanLassoRadius=375,
                             minDmrSep=1000,
                             minDmrSize=50,
                             adjPvalProbe=1,
                             Rplot=T,
                             PDFplot=T,
                             resultsDir="./CHAMP_ProbeLasso/",
                             ## following parameters are specifically for DMRcate method.
                             rmSNPCH=T,
                             fdr=0.1,
                             dist=2,
                             mafcut=0.05,
                             lambda=1000,
                             C=2)
?DMR.GUI
DMR.GUI(DMR=myDMR,
        beta=myCombat,
        pheno=myLoad$pd$Sample_Group,
        runDMP=TRUE,
        compare.group=NULL,
        arraytype="450K")

#DMR
# up-hypo : TBX20, FRZB
# down-hyper : ZNF544, HSPA1A

?champ.Block
myBlock <-  champ.Block(beta=myNorm,
                        pheno=myLoad$pd$Sample_Group,
                        arraytype="450K",
                        maxClusterGap=250000,
                        B=500,
                        bpSpan=250000,
                        minNum=10,
                        cores=5)

Block.GUI(Block=myBlock,beta=myNorm,pheno=myLoad$pd$Sample_Group,runDMP=TRUE,compare.group=NULL,arraytype="450K")




?champ.GSEA
myGSEA <-     champ.GSEA(beta=myCombat,
                         DMP=myDMP[[1]],
                         DMR=myDMR,
                         CpGlist=NULL,
                         Genelist=NULL,
                         pheno=myLoad$pd$Sample_Group,
                         method="fisher",
                         arraytype="450K",
                         Rplot=TRUE,
                         adjPval=0.1,
                         cores=1)

View(myGSEA$DMP)

?champ.CNA
myCNA <-   champ.CNA(intensity=myCombat,
                     pheno=myLoad$pd$Sample_Group,
                     control=TRUE,
                     controlGroup="champCtls",
                     sampleCNA=TRUE,
                     groupFreqPlots=TRUE,
                     Rplot=FALSE,
                     PDFplot=TRUE,
                     freqThreshold=0.3,
                     resultsDir="./CHAMP_CNA",
                     arraytype="450K")


DMP<-myDMP$Low_to_High

hist(as.matrix(DMP$logFC))

range(DMP$deltaBeta)


## Low-High 라서 DEG랑 반대로 봐야한다.
hyper<-DMP[DMP$P.Value < 0.05 & DMP$deltaBeta > 0.05,]
hypergene<-as.data.frame(hyper$gene)

hypo<-DMP[DMP$P.Val < 0.05 & DMP$deltaBeta < -0.05,]
hypogene<-as.data.frame(hypo$gene)

write.csv(hypergene,file = "hypergene.csv")
write.csv(hypogene,file = "hypogene.csv")

save.image(file="20220111.RData")

#Input genes; convert to ENTREZID 
eg = bitr(hyper$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg1 = bitr(hypo$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#GO-Enrichment analysis 

go <- enrichGO(eg$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") #keyType ='SYMBOL',pAdjustMethod = "BH",
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") # ont= "BP","CC","MF"

dotplot(go, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=5) + facet_grid(ONTOLOGY~., scale="free")
dotplot(go1, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=5) + facet_grid(ONTOLOGY~., scale="free")
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

dotplot(kg,color='pvalue')
dotplot(kg1,color='pvalue')

View(kg@result)

ggsave(filename ='/home/dongmin/TCGA-COAD,READ/BMI/GO_up_regulated.png',
       device = 'png',
       dpi = 1200
)


dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
plotMA(res)


gseaplot(kg,geneSetID = 1)

