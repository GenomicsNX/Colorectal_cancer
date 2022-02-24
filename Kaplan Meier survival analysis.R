
# 5 years KM-survival analysis

right_dataset[,"OS.month"]<- (right_dataset$OS.time)/30

right_dataset$OS.month
right_dataset[,"OS.month"]<- (right_dataset$OS.time)/30


right_dataset_60month<-right_dataset[right_dataset$OS.month<=60,]

right_dataset_60month<-right_dataset[right_dataset$OS.time<=1826,]

###################################################################



sfit<-survfit(Surv(PFI.time,PFI)~right_Macrophages.M1,data=right_dataset_60month_cox)
plot(sfit)


ggsurv<-ggsurvplot(sfit,pval=TRUE, risk.table=FALSE,pval.size=10,
                   title="Progression free survival",
                   fontsize=10,
                   censor.size=7,
                   legend.labs=c("High", "Low"), legend.title="M1 Macrophages",
                   palette=c("red","dodgerblue2"),
                   font.legend= c(20, "bold", "black"),
                   ggtheme = theme_bw(),
                   xlab = "Time(days)",
                   font.main = c(20, "bold", "black"),
                   font.x = c(20, "bold", "black"),
                   font.y = c(20, "bold", "black"),
                   font.tickslab = c(10, "plain", "black"),
                   
)

ggsurv$plot <- ggsurv$plot + 
  theme(legend.title = element_text(size = 18, color = "black", face = "bold"),
        legend.text = element_text(size = 18, color = "black", face = "bold"),
        axis.text.x = element_text(size = 20, color = "black", face = "plain"),
        axis.text.y = element_text(size = 20, color = "black", face = "plain"),
        axis.title.x = element_text(size = 20, color = "black", face = "plain"),
        axis.title.y = element_text(size = 20, color = "black", face = "plain"))
ggsurv







ggsave(filename ='/home/dongmin/colorectal cancer microbiome/adddata/Right_Surv/PFS_M1macrophages.png',
       device = 'png',
       dpi = 1200,
       width= 10000,
       height = 8000,
       units = "px")



######################################################################
# Cox 회귀분석


right_dataset_60month_cox$gender
right_dataset_60month_cox<-merge(COADREAD_clinical[c(15,4,104,61)],right_dataset_60month,by="PATIENT")
right_dataset_60month_cox<- subset(right_dataset_60month_cox,lymphatic_invasion!="")
right_dataset_60month_cox<- subset(right_dataset_60month_cox,venous_invasion!="")

right_dataset_60month_cox <- right_dataset_60month_cox[!is.na(right_dataset_60month_cox$ajcc_pathologic_stage),]

right_dataset_60month_cox$Disease_stage <-ifelse(right_dataset_60month_cox$ajcc_pathologic_stage%in%
                                                   c("Stage I","Stage IA","Stage IIA", "Stage IIB", "Stage II"),"I-II","III-IV")


#COX Regression
out = coxph(Surv(OS.time,OS)~right+age_at_index+gender+Disease_stage, data=right_dataset_60month_cox)

ggforest(out, data=right_dataset_60month_cox,
         fontsize = 1.2)

summary(out)
out

right_dataset_60month[,c(1311,14)]

ggsave(filename ='/home/dongmin/colorectal cancer microbiome/adddata/riskscore/microbiome/g__Pragia.png',
       device = 'png',
       dpi = 1200,
       width= 10000,
       height = 8000,
       units = "px")

save.image(file="20220107.RData")
