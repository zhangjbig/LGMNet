
#生存分析
library(dplyr)
library(readxl)
library(survival)
library(survminer)
library(GSVA)
library(gridExtra)


markers = readRDS('Test/DEG/data/com_degs_up0.rds')
markers = markers[c(1,2,4,5)]

##CGGA  
####CGGA数据处理####
#subtype
mmc3 <- read_excel("subtype_marker.xls")
subtype_cl_row <- which(mmc3[, 1] == "A. Classical Subtype")
subtype_me_row <- which(mmc3[, 1] == "B. Mesenchymal Subtype")
subtype_pn_row <- which(mmc3[, 1] == "C. Proneural Subtype")
subtype_n_row <- which(mmc3[, 1] == "D. Neural Subtype")
# 提取 CL 数据
subtype_cl <- mmc3[(subtype_cl_row + 2):(subtype_me_row - 3), ]
subtype_me <- mmc3[(subtype_me_row + 2):(subtype_pn_row - 3), ]
subtype_pn <- mmc3[(subtype_pn_row + 2):(subtype_n_row - 3), ]
subtype_n <- mmc3[(subtype_n_row + 2):nrow(mmc3), ]

cl_marker <- unique(trimws(unlist(strsplit(Reduce(union,subtype_cl[,7]), ","))))
me_marker <- unique(trimws(unlist(strsplit(Reduce(union,subtype_me[,7]), ","))))
pn_marker <- unique(trimws(unlist(strsplit(Reduce(union,subtype_pn[,7]), ","))))
n_marker <- unique(trimws(unlist(strsplit(Reduce(union,subtype_n[,7]), ","))))

subtype_markers <- list(CL = cl_marker,
                        ME = me_marker,
                        PN = pn_marker,
                        N = n_marker)
saveRDS(subtype_markers,"subtype_markers.rds")
subtype_markers = readRDS("subtype_markers.rds")

# 确定每个样本的亚型
get_sample_subtypes <- function(gsva_data) {
  # 遍历每一行
  sample_subtypes <- apply(gsva_data, 1, function(row) {
    # 找出得分最高的列名
    max_col <- colnames(gsva_data)[which.max(row)]
    return(max_col)
  })
  return(sample_subtypes)
}

#CGGA325
#exp
tiantan325_exp <- read.table("data/CGGA/325/mRNAseq-readcounts-325gliomas.txt", header = T, sep = '\t')
rownames(tiantan325_exp) <- tiantan325_exp$gene_name
tiantan325_exp <- tiantan325_exp[2:ncol(tiantan325_exp)]
save(tiantan325_exp,file = "data/CGGA/325/tiantan325_exp.Rdata")
#clinical
tiantan325cli <- read.table("data/CGGA/325/CGGA.mRNAseq_325_clinical.20200506.txt",as.is=T,header = T,row.names=1,sep = "\t")
# subtype
tiantan325_exp_gs <- as.matrix(tiantan325_exp)
subtype_gsva <- gsva(tiantan325_exp_gs, subtype_markers, method = "gsva")
subtype_gsva_df <- as.data.frame(t(subtype_gsva))
subtype_tiantan325 <- get_sample_subtypes(subtype_gsva_df)
# 亚型加入临床信息
tiantan325cli <- cbind(tiantan325cli,subtype_tiantan325)
colnames(tiantan325cli)[13] = 'subtype'
tiantan325cli <- tiantan325cli[!is.na(tiantan325cli$PRS_type),]

#合并多变量信息：生存/age/gender/grade/IDH/subtype/1p19q/PRS_type
tiantan325_cli <- cbind(tiantan325cli[,6:7],tiantan325cli[,5],tiantan325cli[,4],tiantan325cli[,3],
                        tiantan325cli[,10],tiantan325cli[,13],tiantan325cli[,11],tiantan325cli[,1])
colnames(tiantan325_cli) <- c("OS(Days)","OS_Censor","age","gender","grade","IDH","subtype","1p19q","PRS_type")
#1为高风险0为低风险
tiantan325_cli$grade <- as.factor(recode(tiantan325_cli$grade, "WHO II" = 2, "WHO III" = 3, "WHO IV" = 4))
tiantan325_cli$IDH <- as.factor(recode(tiantan325_cli$IDH, "Wildtype" = 1, "Mutant" = 0))
tiantan325_cli$`1p19q` <- as.factor(recode(tiantan325_cli$`1p19q`, "Non-codel" = 1, "Codel" = 0))
tiantan325_cli$PRS_type <- as.factor(tiantan325_cli$PRS_type)
tiantan325_cli$gender<-trimws(tiantan325_cli$gender)
save(tiantan325_cli,file = "data/CGGA/325/tiantan325_cli.Rdata")


#CGGA693
#exp
tiantan693_exp <- read.table("data/CGGA/693/mRNAseq-readcounts-693gliomas.txt", header = T, sep = '\t')
rownames(tiantan693_exp) <- tiantan693_exp$gene_name
tiantan693_exp <- tiantan693_exp[2:ncol(tiantan693_exp)]
save(tiantan693_exp,file = "data/CGGA/693/tiantan693_exp.Rdata")
#clinical
tiantan693cli <- read.table("data/CGGA/693/CGGA.mRNAseq_693_clinical.20200506.txt",as.is=T,header = T,row.names=1,sep = "\t")
# subtype
tiantan693_exp_gs <- as.matrix(tiantan693_exp)
subtype_gsva <- gsva(tiantan693_exp_gs, subtype_markers, method = "gsva")
subtype_gsva_df <- as.data.frame(t(subtype_gsva))
subtype_tiantan693 <- get_sample_subtypes(subtype_gsva_df)
# 亚型加入临床信息
tiantan693cli <- cbind(tiantan693cli,subtype_tiantan693)
colnames(tiantan693cli)[13] = 'subtype'
tiantan693cli <- tiantan693cli[!is.na(tiantan693cli$PRS_type),]

#合并多变量信息：生存/age/gender/grade/IDH/subtype/1p19q/PRS_type
tiantan693_cli <- cbind(tiantan693cli[,6:7],tiantan693cli[,5],tiantan693cli[,4],tiantan693cli[,3],
                        tiantan693cli[,10],tiantan693cli[,13],tiantan693cli[,11],tiantan693cli[,1])
colnames(tiantan693_cli) <- c("OS(Days)","OS_Censor","age","gender","grade","IDH","subtype","1p19q","PRS_type")
#1为高风险0为低风险
tiantan693_cli$grade <- as.factor(recode(tiantan693_cli$grade, "WHO II" = 2, "WHO III" = 3, "WHO IV" = 4))
tiantan693_cli$IDH <- as.factor(recode(tiantan693_cli$IDH, "Wildtype" = 1, "Mutant" = 0))
tiantan693_cli$`1p19q` <- as.factor(recode(tiantan693_cli$`1p19q`, "Non-codel" = 1, "Codel" = 0))
tiantan693_cli$PRS_type <- as.factor(tiantan693_cli$PRS_type)
tiantan693_cli$gender<-trimws(tiantan693_cli$gender)
save(tiantan693_cli,file = "data/CGGA/693/tiantan693_cli.Rdata")


####tiantan1018####
load(file = "data/CGGA/325/tiantan325_exp.Rdata")
load(file = "data/CGGA/325/tiantan325_cli.Rdata")
load(file = "data/CGGA/693/tiantan693_exp.Rdata")
load(file = "data/CGGA/693/tiantan693_cli.Rdata")
#1018去掉NA。1014个
tiantan1014_cli <- rbind(tiantan325_cli,tiantan693_cli)
tiantan1018_exp <- cbind(tiantan325_exp,tiantan693_exp)
tiantan1014_exp <- tiantan1018_exp[,colnames(tiantan1018_exp) %in% rownames(tiantan1014_cli)]
#gsva
tiantan1014_exp_gs <- as.matrix(tiantan1014_exp)
tiantan1014_exp_gs <- log2(tiantan1014_exp_gs + 1)
tiantan.gsva.res<-gsva(tiantan1014_exp_gs,markers,method="gsva")#行基因列样本
tiantan.gsva.sc<-as.data.frame(t(tiantan.gsva.res))
tiantan.gsva.sc <- tiantan.gsva.sc[rownames(tiantan1014_cli),] #确保样本名顺序一致
write.csv(tiantan.gsva.sc, 'data/new_cgga_gsva.csv', row.names = TRUE)

#Survival Analysis
#合并生存和基因表达
tiantan1014_survival_bind <- cbind(tiantan1014_cli[,1:2],tiantan.gsva.sc[,1:4])

tiantan1014.sc <- tiantan1014_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan1014.sc$`OS(Days)`), event = as.numeric(tiantan1014.sc$OS_Censor))
data.subset=colnames(tiantan.gsva.sc)
Summtable1=data.frame()
#单变量cox
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tiantan1014.sc)[i+2]<-"clu"
  fit1.mv <- coxph(surv_object  ~clu, data = tiantan1014.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.res=as.data.frame(mv.res)
  Summtable1=rbind(Summtable1,mv.res)
  rownames(Summtable1)[i]=YY
  tiantan1014.sc <- tiantan1014_survival_bind
}
write.csv(Summtable1,"Test/survival/CGGA1014/univariate/tiantan.gsva.univariate.1014samples.csv")
#k-m
Final <- list()
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  tiantan1014.sc <- tiantan1014.sc %>% mutate(Expression.Level = ifelse(tiantan1014.sc[YY]>=median(tiantan1014.sc[,YY]), "Positive", "Negative"))
  tiantan1014.sc$Expression.Level <-factor(tiantan1014.sc$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = tiantan1014.sc)
  XX <- ggsurvplot(fit1, data = tiantan1014.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                   risk.table = T,font.main = c(12, "bold"),legend.title = "status",title=YY, legend.labs = c("Not-Enriched", "Enriched"),
                   font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
  Final[[i]] = XX
}
kmplot <- arrange_ggsurvplots(Final,print = TRUE, ncol = 1,nrow = 4,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
ggsave("Test/survival/CGGA1014/k-m/tiantan.gsva.km.1014samples.png", plot = kmplot, width = 6, height = 25, dpi = 300)

#多变量
#合并生存和基因表达
tiantan1014_survival_bind <- cbind(tiantan1014_cli[,1:7],tiantan1014_cli[,9],tiantan.gsva.sc[,1:4])
colnames(tiantan1014_survival_bind)[8] <- 'PRS_type'
tiantan1014.sc <- tiantan1014_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan1014.sc$`OS(Days)`), event = as.numeric(tiantan1014.sc$OS_Censor))
# gsva/age/gender/grade/IDH/subtype/PRS
rnames=c("gsva","age","gender","grade3","grade4","IDH","subtypeME","subtypeN","subtypePN",'Recurrent','Secondary')
Summtable2=as.data.frame(row.names =rnames,x=rnames)
#cox多变量
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tiantan1014.sc)[i+8]<-"clu"
  fit1.mv <- coxph(surv_object  ~ clu +age+gender+grade+IDH+subtype+PRS_type, data = tiantan1014.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste("Test/survival/CGGA1014/multivariate/tiantan.",data.subset[i],".multivariate.1014sample.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable2[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable2)[i]=YY
  tiantan1014.sc <- tiantan1014_survival_bind
}
write.csv(Summtable2,"Test/survival/CGGA1014/multivariate/tiantan.gsva.multivariate.1014sample.csv")


####23级####
#临床
tiantanlgg_cli <- tiantan1014_cli[tiantan1014_cli$grade %in% c(2,3),] #1014 -> 625个，625个lgg
#表达
tiantan.gsva.sc<-as.data.frame(t(tiantan.gsva.res))
tiantanlgg.gsva <- tiantan.gsva.sc[rownames(tiantanlgg_cli),]

#Survival Analysis
#合并生存和基因表达
tiantan625_survival_bind <- cbind(tiantanlgg_cli[,1:2],tiantanlgg.gsva[,1:4])

tiantan625.sc <- tiantan625_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan625.sc$`OS(Days)`), event = as.numeric(tiantan625.sc$OS_Censor))
data.subset=colnames(tiantan.gsva.sc)
Summtable1=data.frame()
#单变量cox
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tiantan625.sc)[i+2]<-"clu"
  fit1.mv <- coxph(surv_object  ~clu, data = tiantan625.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.res=as.data.frame(mv.res)
  Summtable1=rbind(Summtable1,mv.res)
  rownames(Summtable1)[i]=YY
  tiantan625.sc <- tiantan625_survival_bind
}
write.csv(Summtable1,"Test/survival/CGGA1014/univariate/tiantan.gsva.univariate.625lgg.csv")
#k-m
Final <- list()
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  tiantan625.sc <- tiantan625.sc %>% mutate(Expression.Level = ifelse(tiantan625.sc[YY]>=median(tiantan625.sc[,YY]), "Positive", "Negative"))
  tiantan625.sc$Expression.Level <-factor(tiantan625.sc$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = tiantan625.sc)
  XX <- ggsurvplot(fit1, data = tiantan625.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                   risk.table = T,font.main = c(12, "bold"),legend.title = "status",title=YY, legend.labs = c("Not-Enriched", "Enriched"),
                   font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
  Final[[i]] = XX
}
kmplot <- arrange_ggsurvplots(Final,print = TRUE, ncol = 1,nrow = 4,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
ggsave("Test/survival/CGGA1014/k-m/tiantan.gsva.km.625lgg.png", plot = kmplot, width = 6, height = 25, dpi = 300)

#多变量
#合并生存和基因表达
tiantan625_survival_bind <- cbind(tiantanlgg_cli[,1:7],tiantanlgg_cli[,9],tiantanlgg.gsva[,1:4])
colnames(tiantan625_survival_bind)[8] <- 'PRS_type'
tiantan625.sc <- tiantan625_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan625.sc$`OS(Days)`), event = as.numeric(tiantan625.sc$OS_Censor))
# gsva/age/gender/grade/IDH/subtype
rnames=c("gsva","age","gender","grade3","grade4","IDH","subtypeME","subtypeN","subtypePN",'Recurrent','Secondary')
Summtable2=as.data.frame(row.names =rnames,x=rnames)
#cox多变量
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tiantan625.sc)[i+8]<-"clu"
  fit1.mv <- coxph(surv_object  ~ clu +age+gender+grade+IDH+subtype+PRS_type, data = tiantan625.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste("Test/survival/CGGA1014/multivariate/tiantan.",data.subset[i],".multivariate.625lgg.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable2[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable2)[i]=YY
  tiantan625.sc <- tiantan625_survival_bind
}
write.csv(Summtable2,"Test/survival/CGGA1014/multivariate/tiantan.gsva.multivariate.625lgg.csv")

####测试只有4级####
#临床
tiantangbm_cli <- tiantan1014_cli[tiantan1014_cli$grade==4,] #1014 -> 389
tiantangbm_cli <- tiantangbm_cli[!is.na(tiantangbm_cli$`OS(Days)`),] #去掉NA还有374个
#表达
tiantan.gsva.sc<-as.data.frame(t(tiantan.gsva.res))
tiantangbm.gsva <- tiantan.gsva.sc[rownames(tiantangbm_cli),]

#Survival Analysis
#合并生存和基因表达
tiantan389_survival_bind <- cbind(tiantangbm_cli[,1:2],tiantangbm.gsva[,1:4])

tiantan389.sc <- tiantan389_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan389.sc$`OS(Days)`), event = as.numeric(tiantan389.sc$OS_Censor))
data.subset=colnames(tiantan.gsva.sc)
Summtable1=data.frame()
#单变量cox
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tiantan389.sc)[i+2]<-"clu"
  fit1.mv <- coxph(surv_object  ~clu, data = tiantan389.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.res=as.data.frame(mv.res)
  Summtable1=rbind(Summtable1,mv.res)
  rownames(Summtable1)[i]=YY
  tiantan389.sc <- tiantan389_survival_bind
}
write.csv(Summtable1,"Test/survival/CGGA1014/univariate/tiantan.gsva.univariate.389gbm.csv")
#k-m
Final <- list()
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  tiantan389.sc <- tiantan389.sc %>% mutate(Expression.Level = ifelse(tiantan389.sc[YY]>=median(tiantan389.sc[,YY]), "Positive", "Negative"))
  tiantan389.sc$Expression.Level <-factor(tiantan389.sc$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = tiantan389.sc)
  XX <- ggsurvplot(fit1, data = tiantan389.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                   risk.table = T,font.main = c(12, "bold"),legend.title = "status",title=YY, legend.labs = c("Not-Enriched", "Enriched"),
                   font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
  Final[[i]] = XX
}
kmplot <- arrange_ggsurvplots(Final,print = TRUE, ncol = 1,nrow = 4,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
ggsave("Test/survival/CGGA1014/k-m/tiantan.gsva.km.389gbm.png", plot = kmplot, width = 6, height = 25, dpi = 300)

#多变量
#合并生存和基因表达
tiantan389_survival_bind <- cbind(tiantangbm_cli[,1:7],tiantangbm_cli[,9],tiantangbm.gsva[,1:4])
colnames(tiantan389_survival_bind)[8] <- 'PRS_type'
tiantan389.sc <- tiantan389_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan389.sc$`OS(Days)`), event = as.numeric(tiantan389.sc$OS_Censor))
# gsva/age/gender/grade/IDH/subtype
rnames=c("gsva","age","gender","grade3","grade4","IDH","subtypeME","subtypeN","subtypePN",'Recurrent','Secondary')
Summtable2=as.data.frame(row.names =rnames,x=rnames)
#cox多变量
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tiantan389.sc)[i+8]<-"clu"
  fit1.mv <- coxph(surv_object  ~ clu +age+gender+grade+IDH+subtype+PRS_type, data = tiantan389.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste("Test/survival/CGGA1014/multivariate/tiantan.",data.subset[i],".multivariate.389gbm.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable2[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable2)[i]=YY
  tiantan389.sc <- tiantan389_survival_bind
}
write.csv(Summtable2,"Test/survival/CGGA1014/multivariate/tiantan.gsva.multivariate.389gbm.csv")


####TCGA数据处理####
#exp
load('data/TCGA/RNAseq.log2normcounts-plus-1.RData')  # 667个样本，17301个基因
tcga667_exp <- as.data.frame(exp)
colnames(tcga667_exp) <- gsub("-","_",colnames(tcga667_exp))
save(tcga667_exp,file = "data/TCGA/tcga667_exp.Rdata")

#clinical
load('data/TCGA/glioma.clinicInfo.2017.RData')
AnnotationTCGA<-cbind(clinicInfo$survival,clinicInfo$vital,clinicInfo$age,clinicInfo$gender,clinicInfo$grade,clinicInfo$idhstatus,clinicInfo$codel)
colnames(AnnotationTCGA)<-c("OS(Months)","OS_Censor","age","gender","grade","IDH","1p19q")
AnnotationTCGA<-as.data.frame(AnnotationTCGA)
rownames(AnnotationTCGA) <- gsub("-","_",rownames(AnnotationTCGA))
tcga667_cli<-AnnotationTCGA[colnames(tcga_exp),]
tcga667_cli$`OS(Months)`<-as.numeric(tcga667_cli$`OS(Months)`)
tcga667_cli$age<-as.numeric(tcga667_cli$age)
tcga667_cli$OS_Censor<-as.numeric(tcga667_cli$OS_Censor)
#1为高风险0为低风险
tcga667_cli$grade <- as.numeric(recode(tcga667_cli$grade, "G2" = 2, "G3" = 3, "G4" = 4))
tcga667_cli$IDH <- as.numeric(recode(tcga667_cli$IDH, "WT" = 1, "Mutant" = 0))
tcga667_cli$`1p19q` <- as.numeric(recode(tcga667_cli$`1p19q`, "non-codel" = 1, "codel" = 0))

tcga812_cli <- read_excel("data/TCGA/TCGA812_clinical+info.xlsx")
tcga812_cli <- as.data.frame(tcga812_cli)
rownames(tcga812_cli) <-tcga812_cli$PatientID
tcga812_cli <- tcga812_cli[,2:ncol(tcga812_cli)]
rownames(tcga812_cli) <- gsub("-","_",rownames(tcga812_cli))
tcga667_subtype_cli <- tcga812_cli[rownames(tcga667_cli),]

#合并多变量信息：生存/age/gender/grade/IDH/subtype/1p19q
tcga667_cli <- cbind(tcga667_cli[1:6],tcga667_subtype_cli$`Transcriptome Subtype`,tcga667_cli[,7])
colnames(tcga667_cli)[7:8] <- c("subtype","1p19q")
tcga667_cli$subtype[tcga667_cli$subtype %in% c("NA", "NE")] <- "N"
save(tcga667_cli,file = "data/TCGA/tcga667_cli.Rdata")


####TCGA####
load("data/TCGA/tcga667_exp.Rdata")
load("data/TCGA/tcga667_cli.Rdata")
tcga667_cli$grade <- factor(tcga667_cli$grade)

#gsva
tcga_exp_gs <- as.matrix(tcga667_exp)
tcga.gsva.res<-gsva(tcga_exp_gs,markers,method="gsva")#行基因列样本
# gsva.res <- nes(tcga_exp_gs,markers)
tcga.gsva.sc<-as.data.frame(t(tcga.gsva.res))
tcga.gsva.sc <- tcga.gsva.sc[rownames(tcga667_cli),] #确保样本名顺序一致

#Survival Analysis
#合并生存和基因表达
tcga_survival_bind <- cbind(tcga667_cli[,1:2],tcga.gsva.sc[,1:4])

tcga.sc <- tcga_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga.sc$`OS(Months)`), event = as.numeric(tcga.sc$OS_Censor))
data.subset=colnames(tcga.gsva.sc)
Summtable1=data.frame()
#单变量cox
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tcga.sc)[i+2]<-"clu"
  fit1.mv <- coxph(surv_object  ~clu, data = tcga.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.res=as.data.frame(mv.res)
  Summtable1=rbind(Summtable1,mv.res)
  rownames(Summtable1)[i]=YY
  tcga.sc <- tcga_survival_bind
}
write.csv(Summtable1,"Test/survival/TCGA667/univariate/tcga.gsva.univariate.667samples.csv")
#k-m
Final <- list()
# tcga.sc <- tcga_survival_bind
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  tcga.sc <- tcga.sc %>% mutate(Expression.Level = ifelse(tcga.sc[YY]>=median(tcga.sc[,YY]), "Positive", "Negative"))
  tcga.sc$Expression.Level <-factor(tcga.sc$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = tcga.sc)
  XX <- ggsurvplot(fit1, data = tcga.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                   risk.table = T,font.main = c(12, "bold"),legend.title = "status",title=YY, legend.labs = c("Not-Enriched", "Enriched"),
                   font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Months)")
  Final[[i]] = XX
}

kmplot <- arrange_ggsurvplots(Final,print = TRUE, ncol = 1,nrow = 5,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
ggsave("Test/survival/TCGA667/k-m/tcga.gsva.km.667samples.png", plot = kmplot, width = 6, height = 25, dpi = 300)

#多变量
#Survival Analysis
#合并生存和基因表达
tcga_survival_bind <- cbind(tcga667_cli[,1:7],tcga.gsva.sc[,1:4])

tcga.sc <- tcga_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga.sc$`OS(Months)`), event = as.numeric(tcga.sc$OS_Censor))
data.subset <- colnames(tcga.gsva.sc)
# gsva/age/gender/grade/IDH/MGMT/1p19q/KPS/subtype
rnames=c("gsva","age","gender","grade3","grade4","IDH","subtypeME","subtypeN","subtypePN") #
Summtable2=as.data.frame(row.names =rnames,x=rnames)
#cox多变量
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tcga.sc)[i+7]<-"clu"
  
  fit1.mv <- coxph(surv_object  ~clu+age+gender+grade+IDH+subtype, data = tcga.sc) #
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste("Test/survival/TCGA667/multivariate/tcga.",data.subset[i],".multivariate.667sample.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable2[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable2)[i]=YY
  tcga.sc <- tcga_survival_bind
}
write.csv(Summtable2,"Test/survival/TCGA667/multivariate/tcga.gsva.multivariate.667sample.csv")

####只留23级LGG####
tcgalgg_cli =tcga667_cli[which(tcga667_cli$grade %in% c(2,3)),] #lgg454个

tcga.gsva.sc<-as.data.frame(t(tcga.gsva.res))
tcgalgg.gsva <- tcga.gsva.sc[rownames(tcgalgg_cli),] #确保样本名顺序一致

#Survival Analysis
#合并生存和基因表达
tcga454_survival_bind <- cbind(tcgalgg_cli[,1:2],tcgalgg.gsva[,1:4])

tcga454.sc <- tcga454_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga454.sc$`OS(Months)`), event = as.numeric(tcga454.sc$OS_Censor))
data.subset=colnames(tcgalgg.gsva)
Summtable1=data.frame()
#单变量cox
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tcga454.sc)[i+2]<-"clu"
  fit1.mv <- coxph(surv_object  ~clu, data = tcga454.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.res=as.data.frame(mv.res)
  Summtable1=rbind(Summtable1,mv.res)
  rownames(Summtable1)[i]=YY
  tcga454.sc <- tcga454_survival_bind
}
write.csv(Summtable1,"Test/survival/TCGA667/univariate/tcga.lgg.univariate.454lgg.csv")
#k-m
Final <- list()
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  tcga454.sc <- tcga454.sc %>% mutate(Expression.Level = ifelse(tcga454.sc[YY]>=median(tcga454.sc[,YY]), "Positive", "Negative"))
  tcga454.sc$Expression.Level <-factor(tcga454.sc$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = tcga454.sc)
  XX <- ggsurvplot(fit1, data = tcga454.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                   risk.table = T,font.main = c(12, "bold"),legend.title = "status",title=YY, legend.labs = c("Not-Enriched", "Enriched"),
                   font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
  Final[[i]] = XX
}
kmplot <- arrange_ggsurvplots(Final,print = TRUE, ncol = 1,nrow = 5,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
ggsave("Test/survival/tcga667/k-m/tcga.lgg.km.454lgg.png", plot = kmplot, width = 6, height = 25, dpi = 300)


#多变量
#Survival Analysis
#合并生存和基因表达
tcga454_survival_bind <- cbind(tcgalgg_cli[,1:7],tcgalgg.gsva[,1:4])

tcga454.sc <- tcga454_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga454.sc$`OS(Months)`), event = as.numeric(tcga454.sc$OS_Censor))
# gsva/age/gender/grade/IDH/subtype
rnames=c("gsva","age","gender","grade3","grade4","IDH","subtypeME","subtypeN","subtypePN")
Summtable2=as.data.frame(row.names =rnames,x=rnames)
#cox多变量
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tcga454.sc)[i+7]<-"clu"
  fit1.mv <- coxph(surv_object  ~ clu +age+gender+grade+IDH+subtype, data = tcga454.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste("Test/survival/TCGA667/multivariate/tcga.lgg.",data.subset[i],".multivariate.454lgg.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable2[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable2)[i]=YY
  tcga454.sc <- tcga454_survival_bind
}
write.csv(Summtable2,"Test/survival/TCGA667/multivariate/tcga.lgg.multivariate.454lgg.csv")

####只留4级GBM####
tcgagbm_cli =tcga667_cli[which(tcga667_cli$grade %in% 4),]  #153个
tcga.gsva.sc<-as.data.frame(t(tcga.gsva.res))
tcgagbm.gsva <- tcga.gsva.sc[rownames(tcgagbm_cli),] #确保样本名顺序一致

#Survival Analysis
#合并生存和基因表达
tcga153_survival_bind <- cbind(tcgagbm_cli[,1:2],tcgagbm.gsva[,1:4])

tcga153.sc <- tcga153_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga153.sc$`OS(Months)`), event = as.numeric(tcga153.sc$OS_Censor))
data.subset=colnames(tcgagbm.gsva)
Summtable1=data.frame()
#单变量cox
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tcga153.sc)[i+2]<-"clu"
  fit1.mv <- coxph(surv_object  ~clu, data = tcga153.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.res=as.data.frame(mv.res)
  Summtable1=rbind(Summtable1,mv.res)
  rownames(Summtable1)[i]=YY
  tcga153.sc <- tcga153_survival_bind
}
write.csv(Summtable1,"Test/survival/TCGA667/univariate/tcga.gbm.univariate.153gbm.csv")
#k-m
Final <- list()
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  Expression.Level=paste(YY,"Levels",sep=".")
  tcga153.sc <- tcga153.sc %>% mutate(Expression.Level = ifelse(tcga153.sc[YY]>=median(tcga153.sc[,YY]), "Positive", "Negative"))
  tcga153.sc$Expression.Level <-factor(tcga153.sc$Expression.Level)
  fit1 <- survfit(surv_object ~Expression.Level, data = tcga153.sc)
  XX <- ggsurvplot(fit1, data = tcga153.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                   risk.table = T,font.main = c(12, "bold"),legend.title = "status",title=YY, legend.labs = c("Not-Enriched", "Enriched"),
                   font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
  Final[[i]] = XX
}
kmplot <- arrange_ggsurvplots(Final,print = TRUE, ncol = 1,nrow = 5,surv.plot.height = NULL,risk.table.height = NULL,ncensor.plot.height = NULL)
ggsave("Test/survival/tcga667/k-m/tcga.gbm.km.153gbm.png", plot = kmplot, width = 6, height = 25, dpi = 300)


#多变量
#Survival Analysis
#合并生存和基因表达
tcga153_survival_bind <- cbind(tcgagbm_cli[,1:7],tcgagbm.gsva[,1:4])

tcga153.sc <- tcga153_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga153.sc$`OS(Months)`), event = as.numeric(tcga153.sc$OS_Censor))
# gsva/age/gender/grade/IDH/subtype
rnames=c("gsva","age","gender","grade3","grade4","IDH","subtypeME","subtypeN","subtypePN")
Summtable2=as.data.frame(row.names =rnames,x=rnames)
#cox多变量
for( i in 1: length(data.subset)){
  YY= data.subset[i]
  colnames(tcga153.sc)[i+7]<-"clu"
  fit1.mv <- coxph(surv_object  ~ clu +age+gender+grade+IDH+subtype, data = tcga153.sc)
  mv.res <- summary(fit1.mv)$coefficients
  mv.fname <- paste("Test/survival/TCGA667/multivariate/tcga.gbm.",data.subset[i],".multivariate.153gbm.csv",sep="")
  write.csv(mv.res,mv.fname)
  mv.res=as.data.frame(mv.res)
  Summtable2[i]=mv.res$`Pr(>|z|)`
  colnames(Summtable2)[i]=YY
  tcga153.sc <- tcga153_survival_bind
}
write.csv(Summtable2,"Test/survival/TCGA667/multivariate/tcga.gbm.multivariate.153gbm.csv")

####加入影像####
####TRAIN####
load("data/CGGA/Discovery.RData")
##全部样本
# #临床
# tiantan1014_cli
# #表达
# tiantan1014_exp
# tiantan.gsva.sc
#和影像样本取交集
tiantan.img.cli <- tiantan1014_cli[rownames(tiantan1014_cli) %in% tiantan$ID,] #167个
tiantan.img.gsva <- tiantan.gsva.sc[rownames(tiantan.gsva.sc) %in% tiantan$ID,]
#cut_off值
tiantan1014.median.1 <- median(tiantan.gsva.sc[,'cluster1'])
#验证生存差异是否显著
tiantan167_survival_bind <- cbind(tiantan.img.cli[,1:2],tiantan.img.gsva[,1])
colnames(tiantan167_survival_bind)[3] <- 'gsva'
tiantan167.sc <- tiantan167_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan167.sc$`OS(Days)`), event = as.numeric(tiantan167.sc$OS_Censor))
#km曲线
tiantan167.sc <- tiantan167.sc %>% mutate('Expression.Level' = ifelse(tiantan167.sc$gsva>=tiantan1014.median.1, "Positive", "Negative"))
tiantan167.sc$Expression.Level <-factor(tiantan167.sc$Expression.Level)
fit1 <- survfit(surv_object ~Expression.Level, data = tiantan167.sc)
XX1 <- ggsurvplot(fit1, data = tiantan167.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                  risk.table = T,font.main = c(12, "bold"),legend.title = "status",title='cluster1_img', legend.labs = c("Not-Enriched", "Enriched"),
                  font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
ggsave(filename = "survival/CGGA1014/k-m/tiantan.gsva.km.167img.png", 
       plot = arrangeGrob(XX1$plot,XX1$table,nrow=2,heights = c(3, 1)), 
       width = 5, height = 5, dpi = 300)



#只留23级LGG
# #临床
# tiantanlgg_cli
# #表达
# tiantanlgg.gsva
#和影像样本取交集，tiantan109个
tiantan.lgg.img.cli <- tiantanlgg_cli[rownames(tiantanlgg_cli) %in% tiantan$ID,] #109个
tiantan.lgg.img.gsva <- tiantanlgg.gsva[rownames(tiantanlgg.gsva) %in% tiantan$ID,]
#cut_off值
tiantan625.median.1 <- median(tiantanlgg.gsva[,'cluster1'])
#验证生存差异是否显著
tiantan109_survival_bind <- cbind(tiantan.lgg.img.cli[,1:2],tiantan.lgg.img.gsva[,1])
colnames(tiantan109_survival_bind)[3] <- 'gsva'
tiantan109.sc <- tiantan109_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tiantan109.sc$`OS(Days)`), event = as.numeric(tiantan109.sc$OS_Censor))
tiantan109.sc <- tiantan109.sc %>% mutate('Expression.Level' = ifelse(tiantan109.sc$gsva>=tiantan.median.1, "Positive", "Negative"))
tiantan109.sc$Expression.Level <-factor(tiantan109.sc$Expression.Level)
fit1 <- survfit(surv_object ~Expression.Level, data = tiantan109.sc)
XX1 <- ggsurvplot(fit1, data = tiantan109.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                  risk.table = T,font.main = c(12, "bold"),legend.title = "status",title='cluster1_LGG_img', legend.labs = c("Not-Enriched", "Enriched"),
                  font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
ggsave(filename = "survival/CGGA1014/k-m/tiantan.gsva.km.109imglgg.png", 
       plot = arrangeGrob(XX1$plot,XX1$table,nrow=2,heights = c(3, 1)), 
       width = 5, height = 5, dpi = 300)

# # 把样本高风险和低风险单提出来
# tiantan.high <- rownames(tiantan109.sc[which(tiantan109.sc$Expression.Level == "Positive"), ])
# tiantan.low <- rownames(tiantan109.sc[which(tiantan109.sc$Expression.Level == "Negative"), ])


####TEST####
load("data/TCGA/TCGA.RData")
#处理样本名
# 定义一个函数来处理单个 ID
process_id <- function(id) {
  parts <- strsplit(id, "_")[[1]]
  if (nchar(parts[2]) == 1) {
    parts[2] <- sprintf("%02d", as.numeric(parts[2]))
  }
  if (nchar(parts[3]) < 4) {
    parts[3] <- sprintf("%04d", as.numeric(parts[3]))
  }
  new_id <- paste(parts, collapse = "_")
  return(new_id)
}
TCGA$ID <- sapply(TCGA$ID, process_id)
save(TCGA,file = "data/TCGA/TCGA.RData")
##全部样本
# #临床
# tcga667_cli
# #表达
# tcga667_exp
# tcga.gsva.sc
#和影像样本取交集
tcga.img.cli <- tcga667_cli[rownames(tcga667_cli) %in% TCGA$ID,] #131个
save(tcga.img.cli,file='data/TCGA/tcga_img_cli131.Rdata')
tcga.img.gsva <- tcga.gsva.sc[rownames(tcga.gsva.sc) %in% TCGA$ID,]
#cut_off值
tcga667.median.1 <- median(tcga.gsva.sc[,'cluster1'])
#验证生存差异是否显著
tcga131_survival_bind <- cbind(tcga.img.cli[,1:2],tcga.img.gsva[,1])
colnames(tcga131_survival_bind)[3] <- 'gsva'
tcga131.sc <- tcga131_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga131.sc$`OS(Months)`), event = as.numeric(tcga131.sc$OS_Censor))
#自己cutoff
tcga131.sc <- tcga131.sc %>% mutate('Expression.Level' = ifelse(tcga131.sc$gsva>=tcga667.median.1, "Positive", "Negative"))
tcga131.sc$Expression.Level <-factor(tcga131.sc$Expression.Level)
fit1 <- survfit(surv_object ~Expression.Level, data = tcga131.sc)
XX1 <- ggsurvplot(fit1, data = tcga131.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                  risk.table = T,font.main = c(12, "bold"),legend.title = "status",title='cluster1_img_self_cutoff', legend.labs = c("Not-Enriched", "Enriched"),
                  font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
ggsave(filename = "survival/TCGA667/k-m/tcga_131img_self_667cutoff.png", 
       plot = arrangeGrob(XX1$plot,XX1$table,nrow=2,heights = c(3, 1)), 
       width = 5, height = 5, dpi = 300)
#train的cutoff
tcga131.sc <- tcga131_survival_bind
tcga131.sc <- tcga131.sc %>% mutate('Expression.Level' = ifelse(tcga131.sc$gsva>=tiantan1014.median.1, "Positive", "Negative"))
tcga131.sc$Expression.Level <-factor(tcga131.sc$Expression.Level)
fit1 <- survfit(surv_object ~Expression.Level, data = tcga131.sc)
XX2 <- ggsurvplot(fit1, data = tcga131.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                  risk.table = T,font.main = c(12, "bold"),legend.title = "status",title='cluster1_img_train_cutoff', legend.labs = c("Not-Enriched", "Enriched"),
                  font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
ggsave(filename = "survival/TCGA667/k-m/tcga_131img_train1014_cutoff.png", 
       plot = arrangeGrob(XX2$plot,XX2$table,nrow=2,heights = c(3, 1)), 
       width = 5, height = 5, dpi = 300)



#只留23级LGG
# #临床
# tcgalgg_cli
# #表达
# tcgalgg.gsva
#和影像样本取交集
tcgalgg.img.cli <- tcgalgg_cli[rownames(tcgalgg_cli) %in% TCGA$ID,] #93个
tcgalgg.img.gsva <- tcgalgg.gsva[rownames(tcgalgg.gsva) %in% TCGA$ID,]
#cut_off值
tcga454.median.1 <- median(tcgalgg.gsva[,'cluster1'])
#验证生存差异是否显著
tcga93_survival_bind <- cbind(tcgalgg.img.cli[,1:2],tcgalgg.img.gsva[,1])
colnames(tcga93_survival_bind)[3] <- 'gsva'
tcga93.sc <- tcga93_survival_bind
#生存事件
surv_object <- Surv(time = as.numeric(tcga93.sc$`OS(Months)`), event = as.numeric(tcga93.sc$OS_Censor))
#自己cutoff
tcga93.sc <- tcga93.sc %>% mutate('Expression.Level' = ifelse(tcga93.sc$gsva>=tcga454.median.1, "Positive", "Negative"))
tcga93.sc$Expression.Level <-factor(tcga93.sc$Expression.Level)
fit1 <- survfit(surv_object ~Expression.Level, data = tcga93.sc)
XX3 <- ggsurvplot(fit1, data = tcga93.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                  risk.table = T,font.main = c(12, "bold"),legend.title = "status",title='cluster1_imgLGG_self_cutoff', legend.labs = c("Not-Enriched", "Enriched"),
                  font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
ggsave(filename = "survival/TCGA667/k-m/tcga_93img_self_454cutoff.png", 
       plot = arrangeGrob(XX3$plot,XX3$table,nrow=2,heights = c(3, 1)), 
       width = 5, height = 5, dpi = 300)
#train的cutoff
tcga93.sc <- tcga93_survival_bind
tcga93.sc <- tcga93.sc %>% mutate('Expression.Level' = ifelse(tcga93.sc$gsva>=tiantan625.median.1, "Positive", "Negative"))
tcga93.sc$Expression.Level <-factor(tcga93.sc$Expression.Level)
fit1 <- survfit(surv_object ~Expression.Level, data = tcga93.sc)
XX4 <- ggsurvplot(fit1, data = tcga93.sc, pval = TRUE,pval.cExpression.Levelrd = c(50, 1),pval.size=6,legend="top",palette = c("Blue","Red"),
                  risk.table = T,font.main = c(12, "bold"),legend.title = "status",title='cluster1_imgLGG_train_cutoff', legend.labs = c("Not-Enriched", "Enriched"),
                  font.x = c(12, "bold"),font.y = c(12, "bold"),font.tickslab = c(12, "plain"),font.legend = c(12),xlab="Time (Days)")
ggsave(filename = "survival/TCGA667/k-m/tcga_93img_train625_cutoff.png", 
       plot = arrangeGrob(XX4$plot,XX4$table,nrow=2,heights = c(3, 1)), 
       width = 5, height = 5, dpi = 300)
####保存分组数据####
save(tiantan167.sc,tiantan109.sc,tcga131.sc,tcga93.sc, file = 'survival/tiantan_tcga_sc.Rdata')




####保存训练样本####
##全部样本
#train
train167.level <- tiantan167.sc[,3:4]
train167.level$ID <- rownames(train167.level)
colnames(train167.level)<- c('level','label','ID')
train167.level$label = as.factor(recode(train167.level$label, 'Positive'=1, 'Negative'=0))
train167.img <- tiantan[,which(colnames(tiantan) == 'exponential_firstorder_10Percentile'):ncol(tiantan)]
train167.img <- train167.img[,-c(265,269,270,355)] #去除特征值相同的列
train167.img <- cbind(tiantan$ID,train167.img)
colnames(train167.img)[1] = 'ID'
train167 <- merge(train167.level, train167.img, by = "ID")
write.csv(train167,'data/img/train167.csv',row.names = FALSE)

#test
test131.level <- tcga131.sc[,3:4]
test131.level$ID <- rownames(test131.level)
colnames(test131.level)<- c('level','label','ID')
test131.level$label = as.factor(recode(test131.level$label, 'Positive'=1, 'Negative'=0))
test131.img <- TCGA[,which(colnames(TCGA) == 'exponential_firstorder_10Percentile'):ncol(tiantan)]
test131.img <- test131.img[,-c(265,269,270,355)] #去除特征值相同的列
test131.img <- cbind(TCGA$ID,test131.img)
colnames(test131.img)[1] = 'ID'
test131 <- merge(test131.level, test131.img, by = "ID")
write.csv(test131,'data/img/test131.csv', row.names = FALSE)


##LGG
#train
train109.level <- tiantan109.sc[,3:4]
train109.level$ID <- rownames(train109.level)
colnames(train109.level)<- c('level','label','ID')
train109.level$label = as.factor(recode(train109.level$label, 'Positive'=1, 'Negative'=0))
train109.img <- tiantan[,which(colnames(tiantan) == 'exponential_firstorder_10Percentile'):ncol(tiantan)]
train109.img <- train109.img[,-c(265,269,270,355)] #去除特征值相同的列
train109.img <- cbind(tiantan$ID,train109.img)
colnames(train109.img)[1] = 'ID'
train109 <- merge(train109.level, train109.img, by = "ID")
write.csv(train109,'data/img/train109.csv',row.names = FALSE)

#test
test93.level <- tcga93.sc[,3:4]
test93.level$ID <- rownames(test93.level)
colnames(test93.level)<- c('level','label','ID')
test93.level$label = as.factor(recode(test93.level$label, 'Positive'=1, 'Negative'=0))
test93.img <- TCGA[,which(colnames(TCGA) == 'exponential_firstorder_10Percentile'):ncol(tiantan)]
test93.img <- test93.img[,-c(265,269,270,355)] #去除特征值相同的列
test93.img <- cbind(TCGA$ID,test93.img)
colnames(test93.img)[1] = 'ID'
test93 <- merge(test93.level, test93.img, by = "ID")
write.csv(test93,'data/img/test93.csv', row.names = FALSE)
