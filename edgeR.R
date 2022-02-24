#edgeR


#edgeR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")




library(edgeR)
install.packages("glmnet", repos = "https://cran.us.r-project.org")
library(glmnet)


library(edgeR)
y <- DGEList(counts=eData,group=metaData$Group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~metaData$Group)

write.csv(design,file = "design.csv")
design<-read.csv(file = "design.csv")

y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)

topTags(qlf)
?glmQLFTest
qlf$table

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

lrt$table

res<-qlf$table



up_regulated<-res[res$PValue < 0.05 & res$logFC > 0.5,]
down_regulated<-res[res$PValue < 0.05 & res$logFC < -0.5,]



rownames(up_regulated)
rownames(down_regulated)


nrow(up_regulated)
nrow(down_regulated)


write.csv(up_regulated,file = "/home/dongmin/TCGA-COAD,READ/BMI/up_regulated.csv")
write.csv(down_regulated,file = "/home/dongmin/TCGA-COAD,READ/BMI/down_regulated.csv")

dd<-read.csv(file = "dd.csv")


#Input genes; convert to ENTREZID 
eg = bitr(rownames(up_regulated), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg1 = bitr(rownames(down_regulated), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#GO-Enrichment analysis 

go <- enrichGO(eg$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") #keyType ='SYMBOL',pAdjustMethod = "BH",
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") # ont= "BP","CC","MF"

dotplot(go, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=10, label_format = 50) + facet_grid(ONTOLOGY~., scale="free")
dotplot(go1, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=10, label_format = 50) + facet_grid(ONTOLOGY~., scale="free")
##############################
go <- enrichGO(eg$SYMBOL,keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1,pAdjustMethod = "BH") #keyType ='SYMBOL',pAdjustMethod = "BH",
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=1) # ont= "BP","CC","MF"

dotplot(go, split="ONTOLOGY",color="pvalue", font.size =10,showCategory=7) + facet_grid(ONTOLOGY~., scale="free")

View(go@result)
?EnhancedVolcano
library(EnhancedVolcano)
EnhancedVolcano(res,
                title = "Volcano plot",
                subtitle = bquote(italic("BMI>=30 vs BMI<=25")),
                lab = "",
                x = 'logFC',
                y = 'PValue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 5.0,
                #colAlpha = 1,
                #legend=c('NS','Log (base 2) fold-change','P value',
                 #        'P value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)
?ggsave

ggsave(filename ='/home/dongmin/TCGA-COAD,READ/BMI/BMIedgeR_vol.png',
       device = 'png',
       dpi = 1200,
       width=9000,
       height = 5850,
       units = "px"
)


?enrichKEGG


kg <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 1)

kg1 <- enrichKEGG(gene         = eg1$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 1)

dotplot(kg,color='pvalue', font.size =13,showCategory=10, label_format = 50)
dotplot(kg1,color='pvalue', font.size =13,showCategory=10, label_format = 50)

View(kg@result)

ggsave(filename ='/home/dongmin/TCGA-COAD,READ/BMI/GO_up_regulated.png',
       device = 'png',
       dpi = 1200
)


dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
plotMA(res)


gseaplot(kg,geneSetID = 1)




save.image(file = "20220111.RData")


test<-plotMDS(y,top=20000)


#hub genes

PT_hub <- read.csv(file = "./PT_hub.csv")


# HighUp, HighDown, LowUp, LowDown
#Input genes; convert to ENTREZID 
eg = bitr(PT_hub$Up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg1 = bitr(PT_hub$Down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg2= bitr(PT_hub$hubgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#GO-Enrichment analysis 

go <- enrichGO(eg$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=0.05,pAdjustMethod = "BH") #keyType ='SYMBOL',pAdjustMethod = "BH",
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=0.05,pAdjustMethod = "BH") # ont= "BP","CC","MF"
go2 <- enrichGO(eg2$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff=0.05,pAdjustMethod = "BH")

dotplot(go, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=10, label_format = 80) + facet_grid(ONTOLOGY~., scale="free")
dotplot(go1, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=10,label_format = 100) + facet_grid(ONTOLOGY~., scale="free")
dotplot(go2, split="ONTOLOGY",color="pvalue", font.size =13,showCategory=10,label_format = 100) + facet_grid(ONTOLOGY~., scale="free")


kg <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kg1 <- enrichKEGG(gene         = eg1$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
kg2 <- enrichKEGG(gene         = eg2$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

dotplot(kg, color='pvalue',font.size =13,showCategory=10, label_format = 50)
dotplot(kg1,color='pvalue',font.size =13,showCategory=10, label_format = 50)
dotplot(kg2,color='pvalue',font.size =13,showCategory=10, label_format = 50)

View(kg@result)

ggsave(filename ='/home/dongmin/TCGA-COAD,READ/BMI/KEGG_up-hypo.png',
       device = 'png',
       dpi = 1200,
       width = 11000,
       height = 13000,
       units = "px"
)


save.image(file = "20220111.RData")
