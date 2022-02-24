


#윌콕스 단측검정. 대립가설 왼쪽 < 오늘쪽 (평균값)
wilcox.test(high_CIBERSORTx$Neutrophils,low_CIBERSORTx$Neutrophils,alternative="less")

#윌콕스 단측검정. 대립가설 왼쪽 > 오늘쪽 (평균값)
wilcox.test(high_CIBERSORTx$Neutrophils,low_CIBERSORTx$Neutrophils,alternative="greater")

##윌콕스 양측검정.
wilcox.test(high_CIBERSORTx$Neutrophils,low_CIBERSORTx$Neutrophils,alternative="two.sided")


summary(high_CIBERSORTx$T.cells.follicular.helper)

summary(low_CIBERSORTx$T.cells.follicular.helper)

