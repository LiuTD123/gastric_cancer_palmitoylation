rm(list = ls())
foldpath <- paste0("D:/workdir/12stadb/08borutax")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
set.seed(125)
library(glmnet)
library(dplyr)
library(ggplot2)
library(tidyverse)
# ---------------------------------------机器学习部分
# 
library(xgboost)
library(randomForest)
library(caret)
library(DALEX)
library(stats)
library(e1071)
#  # # ------------------Boruta -------

library(Boruta)
train_expr <-readRDS("..\\00data\\tcga_count.rds")
train_group <- readRDS("..\\00data\\tcga_group.rds")
# group <- select(group,-1)
sig_expore <- read.csv('../03relatgene/hub_gene.csv')

hub_gene <- sig_expore$x

data <- train_expr[hub_gene,]
data <- na.omit(data)
group <- train_group

# hub_gene <- read.table("/data/nas1/zhuxuying/04.YQNN-10301-9/08_Expression/hub.csv",sep = ",",header = T)

# names(group)[2] <- "gene"
identical(colnames(data),group$sample)

data <- as.data.frame(t(data))
# data <- data[,hub_gene]
data$sample <- rownames(data)

dat <- merge(data,group,by = "sample")

# dat <- dat[,-10]
dat <- dat[order(dat$group,decreasing = F),]
# names(dat)[1] <- "sample"

# names(dat)[13] <- "group"
rownames(dat) <- dat$sample
dat1 <- dat[-1]

dat1$group <- factor(dat1$group,levels = c('Tumor','Normal'))

num <- length(colnames(dat1))-1
dat1[,1:num] <- as.data.frame(lapply(dat1[,1:num],as.numeric))
dat1 <- na.omit(dat1)

res.Boruta<-Boruta(x=dat1[,2: ncol(dat1)-1], y=as.factor(dat1[,ncol(dat1)]), pValue=0.01, mcAdj=T,
                   maxRuns=1000)
Boruta<-attStats(res.Boruta) #给出Boruta算法的结果
write.csv(Boruta,"09.Boruta.csv") #,head = T
table(res.Boruta$finalDecision)
# Tentative Confirmed  Rejected 
# 9        63        39 

Boruta <- Boruta[order(Boruta$medianImp, decreasing = T),]
boruta_geneids<-Boruta[Boruta$decision=='Confirmed',]%>%rownames(.)
boruta_geneids
# [1] "ANP32B"   "ECT2"     "CD44"     "NASP"     "RNASEH1"  "TMPO"     "CDC45"    "LIN9"     "BARD1"    "SLC25A17" "CCDC77"   "LMNB1"   
# [13] "CENPH"    "MAD2L1"   "MTHFD2"   "TTLL4"    "RRM2"     "TCF19"    "AP3M2"    "IPO9"     "MSH6"     "PSMC3IP"  "CENPO"    "GPR19"   
# [25] "PSRC1"    "HELLS"    "ZNF473"   "CPOX"     "GTPBP4"   "RFT1"     "MTHFD1"   "EZH2"     "PRR11"    "ZWILCH"   "KIF23"    "RAD54L"  
# [37] "FIGNL1"   "NUDT5"    "TTK"      "EARS2"    "KNTC1"    "PRIM1"    "PARPBP"   "DNA2"     "OIP5"     "NDC80"    "NAT10"    "SPAG5"   
# [49] "RRP1B"    "WDHD1"    "CHEK2"    "KDM1A"    "KNSTRN"   "NUP155"  

##定义一个函数提取每个变量对应的重要性值。

library(dplyr)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  # sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
  #   summarise(median=median(Importance)) %>% arrange(median)
  # sortedVariable <- as.vector(sortedVariable$Variable)
  # 
  # 
  # boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(res.Boruta)

head(boruta.variable.imp)
boruta.variable.imp <- boruta.variable.imp[]

# 绘制Boruta算法运行过程中各个变量的重要性得分的变化 （绿色是重要的变量，红色是不重要的变量，蓝色是影子变量，黄色是Tentative变量）
# 用来查看是否有必要增加迭代的次数以便再次确认Tentative变量中是否有一部分为有意义的特征变量。
# ，黄色变量部分随着迭代还是有部分可能高于最高值，可以继续尝试增加迭代次数。

pdf("boruta_imphistory.pdf",w=5,h=4)
Boruta::plotImpHistory(res.Boruta)
dev.off()

# 提取重要变量
boruta.finalVars <- data.frame(Item=getSelectedAttributes(res.Boruta, withTentative = F), Type="Boruta")

# devtools::install_github("Tong-Chen/YSX")
# library(YSX)
# devtools::install_github("Tong-Chen/ImageGP")
library(ImageGP)

p <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                legend_variable = "finalDecision", legend_variable_order = c("Tentative", "Confirmed", "Rejected", "shadowMax", "shadowMean", "shadowMin"),
                xtics_angle = 90,
                xvariable_order = c(rownames(Boruta)))

# ggsave(p,file = "09.Boruta.pdf",w = 18, h = 8,family='Times')
# ggsave(p,file = "09.Boruta.png",w = 18, h = 8,family='Times')
pdf("09.Boruta.pdf",w = 18, h = 8,family='Times')
p
dev.off()

write.table(boruta_geneids,"10.boruta_geneids.csv",sep = "\t",row.names = F,col.names = T,quote = F)
