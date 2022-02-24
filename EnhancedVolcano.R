BiocManager::install('EnhancedVolcano',force = TRUE)

library(EnhancedVolcano)





ggsave(filename ='/home/dongmin/TCGA-COAD,READ/BMI/EnhancedVolcano.png',
       device = 'png',
       dpi = 1200
)

EnhancedVolcano(res,
                lab = rownames(res),
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
