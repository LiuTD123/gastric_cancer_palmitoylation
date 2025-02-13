rm(list = ls())

foldpath <- paste0("D:/workdir/12stadb/05model")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE))
  devtools::install_github("souravc83/fastAdaboost")

# if (!requireNamespace("Mime", quietly = TRUE))
#   devtools::install_github("l-magnificence/Mime")

library(Mime1)
library(readxl)
# BiocManager::install("mice")
library(mice)

# -------OS----------
# ----------------训练集--------------------
data_train <- readRDS("../00data/tcga_fpkm.rds")
# --------------
# library(tibble)
# library(magrittr)
# # 加载edgeR包
# library(edgeR)
# # 假设你有一个计数矩阵counts，行名为基因，列名为样本
# # 例如：counts <- matrix(c... , nrow=..., dimnames=list(genes, samples))
# 
# # 创建DGEList对象
# y <- DGEList(counts=data_train)
# 
# # 计算文库大小
# y$samples$lib.size <- colSums(y$counts)
# 
# # 计算归一化因子（TMM方法）
# y <- calcNormFactors(y)
# 
# # 计算FPKM
# y_fpkm <- cpm(y, normalized.lib.sizes = TRUE)
# y_fpkm <- y_fpkm / 1000 * 1e9 / y$samples$lib.size
# 
# # 查看结果
# print(head(y_fpkm))
# data_train <- y_fpkm
# -------------
data_surv <- readRDS("../00data/survival_dat.rds")

data_train <- as.data.frame(t(data_train))
data_train$sample <- rownames(data_train)

traindata <- merge(data_surv,data_train,by = "sample")
colnames(traindata)[1:3] <- c("ID","OS.time","OS") 
traindata$OS.time <- round(traindata$OS.time*365)
traindata$OS <- as.numeric(traindata$OS)

# ----------------验证集--------------------
dir ="../00geodata"
samples=list.files(dir,pattern = "^GSE.*RData$")

index <- c()
for (sample in samples){
  if (nchar(sample)<17){
    index <- c(index,sample)
  }
}

for (i in index){
  load(paste0(dir,"/",i))
}
# 
# # -------group# ------------------------------------------
# # load("./exampledata/Example.cohort.Rdata") # 生存数据与基因表达信息
# # load("./exampledata/genelist.Rdata")
# # list_train_vali_Data[[5]][1:5,1:5]$OS
# #                 ID    OS.time OS   MT-CO1   MT-CO3
# #60  TCGA.DH.A66B.01 1281.65322  0 13.77340 13.67931
# #234 TCGA.HT.7607.01   96.19915  1 14.96535 14.31857
# #42  TCGA.DB.A64Q.01  182.37755  0 13.90659 13.65321
# #126 TCGA.DU.8167.01  471.97707  0 14.90695 14.59776
# #237 TCGA.HT.7610.01 1709.53901  0 15.22784 14.62756
# # 其中list_train_vali_Data是含有2个数据集的列表，每个数据集的第一列为ID ，2-3列为生存信息（OS.time ，OS） ，后面为基因表达量。
# range(GSE84437$OS.time)
# GSE84437 <- na.omit(GSE84437)
# 
# # -----------------
# GSE84437_exp <- GSE84437[,c(1,4:19588)]
# GSE84437_surv <- GSE84437[,c(1:3)]
# 
# rownames(GSE84437_exp) <-GSE84437_exp$ID
# GSE84437_exp <- GSE84437_exp[,-1]
# 
# y <- DGEList(counts=GSE84437_exp)
# 
# # 计算文库大小
# y$samples$lib.size <- colSums(y$counts)
# 
# # 计算归一化因子（TMM方法）
# y <- calcNormFactors(y)
# 
# # 计算FPKM
# y_fpkm <- cpm(y, normalized.lib.sizes = TRUE)
# y_fpkm <- y_fpkm / 1000 * 1e9 / y$samples$lib.size
# 
# # 查看结果
# print(head(y_fpkm))
# GSE84437_exp<- as.data.frame(y_fpkm)
# 
# GSE84437_exp$ID <- rownames(GSE84437_exp)
# GSE84437 <- merge(GSE84437_surv,GSE84437_exp,by = "ID")
# -------------------------
list_train_vali_Data <- list("Train" = traindata,
                             "GSE13861" = GSE13861,
                             "GSE26942" = GSE26942,
                             "GSE28541" = GSE28541,
                             "GSE29272" = GSE29272,
                             "GSE66229" = GSE66229
                             # "GSE84437" = GSE84437,
                             # "GSE15459" = GSE15459
)

save(list_train_vali_Data, file = "list_train_vali_Data.RData")

load("train_vali_data.RData")
# 二 构建预后模型
# 1. 构建101机器学习模型组合
genes <- read.csv("../08borutax/10.boruta_geneids.csv")
genelist <- genes$x

genelist <- Reduce(intersect,list(colnames(GSE13861),
                                  colnames(GSE15459),
                                  colnames(GSE26942),
                                  colnames(GSE28541),
                                  colnames(GSE29272),
                                  # colnames(GSE66229),
                                  # colnames(GSE84437),
                                  genelist))

save(genelist, file = "intergenes.RData")
# # # # 删除小于100天生存期的样本
for (i in 1:length(list_train_vali_Data)){
  list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][list_train_vali_Data[[i]]$OS.time >90,]
  list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][list_train_vali_Data[[i]]$OS.time <365*10,]
  list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][,c("ID","OS.time","OS",genelist)]
  list_train_vali_Data[[i]] <- na.omit(list_train_vali_Data[[i]])
}

save(list_train_vali_Data, file = "list_train_vali_Data_deal.RData")
# trace("ML.Dev.Prog.Sig",edit = T, where=asNamespace("Mime1"))

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Train,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.1,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize =5,
                       seed = 5201314)

save(res, file = "res_all_ML.Dev.RData")

load("res_all_ML.Dev1.RData")

# 提取系数
ml <- res$ml.res
# ml <- ml$`StepCox[forward] + RSF`
modelgenes <- ml$`RSF + GBM`$fit$var.names

save(ml,file = "ml_details.RData")
write.csv(modelgenes,file = "modelgenes.csv")

# ML.Dev.Prog.Sig() 可选 all, single 和 double三种模式. all 为所有10种算法 以及 组合 . 
# single 为用10种算法中的一种.
# double 为两种算法的组合，一般情况下使用 all 模式.
# 默认情况下 unicox.filter.for.candi 为 T , 会先对训练集进行单因素cox分析，unicox_p_cutoff 显著的基因会用于构建预后模型.

# 如果使用自己数据的时候，需要注意：
# （1）替换自己数据注意前三列的要求，且将多个数据集以列表形式存储
# （2）分析之前最好先确认 所有数据集中是否 有基因集列表中的所有基因 ，减少报错。
# （3）种子数确定好，会有一些小的影响 。

data2 <- data2 %>% 
  dplyr::select(ID , OS.time , OS, genelist)
# 通过View(ML.Dev.Prog.Sig) 检查函数，设置unicox.filter.for.candi = T 后会先做单因素cox分析

# 2. C-index 展示
# 示例数据list_train_vali_Data 为2个数据集的list，结果图中队列为2个，最后两列为Cindex的均值，这也就是机器学习模型组合文献中的主图。
pdf("cindex_dis.pdf",w=9,h=18)
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[1],
               order = names(list_train_vali_Data),
               width = 0.35
)
dev.off()

cindex <- res$Cindex.res
write.csv(cindex,file = "cindex.csv",quote = F)

# 3. 查看指定模型的结果
# 假设我们选择第一个模型（StepCox[forward] + plsRcox） ，可以单独查看该模型下各个数据集的cindex表现

method_select = "StepCox[forward] + plsRcox"

pdf("cindex_dis_select.pdf",w=5,h=5)
cindex_dis_select(res,
                  model= method_select,
                  order= names(list_train_vali_Data))
dev.off()

# 也可以查看该模型下各个数据集的KM曲线情况
survplot <- vector("list",3) 
pdf("km_.pdf",w=8,h=7)
for (i in c(1:3)) {
  survplot[[i]]<-rs_sur(res, model_name = method_select,
                        dataset = names(list_train_vali_Data)[i],
                        # color=c("blue","green"),
                        median.line = "hv",
                        # cutoff = 0.5,
                        conf.int = T,
                        xlab="Day",
                        pval.coord=c(1000,0.9)
  )
  
}
dev.off()

# rs_sur(res, model_name = method_select,
#        dataset = names(list_train_vali_Data)[6],
#        # color=c("blue","green"),
#        median.line = "hv",
#        cutoff = 0.5,
#        conf.int = T,
#        xlab="Day",
#        pval.coord=c(1000,0.9)
# )

pdf("km_combine.pdf",w=16,h=5)
aplot::plot_list(gglist=survplot,ncol=3)
dev.off()

# 提取模型RS结果
# 这里有个很重要的点是要提取指定模型下的RS结果，然后就可以根据自己的需求重新绘制KM 以及 独立预后分析，森林图，列线图等其他分析了。
# 结果都在res中，根据str(res)知道对应的信息，提取即可

head(res$riskscore$`StepCox[forward] + RSF`[[1]])
head(res$riskscore$`StepCox[forward] + RSF`[[2]])

trainriskscore <- res$riskscore[[method_select]][[1]]
testriskscore <- res$riskscore[[method_select]][[2]]
externalriskscore <- res$riskscore[[method_select]][[3]]

save(trainriskscore,testriskscore,externalriskscore,file = "riskscore_superpc.RData")
# 提取特定模型的基因
ml <- res$ml.res
modelgenes <- ml$`StepCox[forward] + RSF`$terms
genes <- colnames(modelgenes)

save(modelgenes,genes,file = "modelgenes.RData")

# ---------提取风险评分-------
method_select = "StepCox[forward] + plsRcox"
rs <- res$riskscore[[method_select]]

save(rs,file = paste0("riskscore_",method_select,".RData"))

# -----------------
# 4. AUC结果
library(survivalROC)
library(survminer)
library(survival)
library(patchwork)
library(timeROC)
# ----------timeROC----------------
dataspan = seq(360,180*10,180)

for (i in 1:length(rs)){
  if (!nrow(list_train_vali_Data[[i]]) == nrow(rs[[i]])){
    sample <- intersect(list_train_vali_Data[[i]]$ID,rs[[i]]$ID)
    list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][list_train_vali_Data[[i]]$ID %in% sample,]
    rs[[i]] <- rs[[i]][rs[[i]]$ID %in% sample,]
  }
}
ROC_merge <- list()

for (i in c(1:6)){
  ROC_merge[[i]] <- timeROC(T=rs[[i]]$OS.time,   
                            delta=rs[[i]]$OS,   
                            marker=rs[[i]]$RS,   
                            cause=1,                #阳性结局指标数值
                            weighting="marginal",   #计算方法，默认为marginal
                            times=dataspan,       #时间点，选取1年，3年和5年的生存率
                            iid=TRUE)
}

# ----------绘图----------------
# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

col=c("#EB4B17", "#2775AB", '#4C8045',"#D8D155","#E0367A","#91612D")
# --------时间推移----------
pdf(paste0("timeauc_",method_select,".pdf"),w=7,h=5)
par(mai=c(1,1,1,2))
plot(dataspan, ROC_merge[[1]]$AUC, lwd=3, type = "l", col = col[1], 
     ylim = c(0.3, 1),
     # xlim = c(360,1200),
     xlab = "Time (days)", ylab = "AUC", main = paste0("AUC"), bty = "l", xaxt = "n")
axis(1, dataspan)
for (i in c(2:6)){
  lines(dataspan, lwd=3, ROC_merge[[i]]$AUC, col = col[i])
}

# legendtxt <- c()
# for (i in 1:6){
#   legendtxt <- c(legendtxt,paste0(names(rs[i])," AUC=",round(mean(ROC_merge[[i]][["AUC"]],na.rm=T),3),
#                                   " CI(",round(mean(confint(ROC_merge[[i]], level = 0.95)$CI_AUC[,1])/100,2),
#                                   "-",round(mean(confint(ROC_merge[[i]], level = 0.95)$CI_AUC[,2])/100,2),")\n"))
# }
legendtxt <- c()
for (i in c(1:6)){
  legendtxt <- c(legendtxt,paste0(names(rs[i])," AUC=",round(mean(ROC_merge[[i]][["AUC"]],na.rm=T),3)))
}
legend(x= 2000,y= 0.8,xpd = TRUE,
       legendtxt,
       col= col,
       lty=1, lwd=3, cex = 0.8)

dev.off()

# -----------------------------------------------------------------
# 计算每个模型的1年，3年，5年 的 auc值 ，并可视化所有模型的1年auc结果
cutoff_1 <- 1
cutoff_2 <- 3
cutoff_3 <- 5
cutoff_4 <- 7
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_1,
                             auc_cal_method="KM")
all.auc.2y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_2,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_3,
                             auc_cal_method="KM")
all.auc.4y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_4,
                             auc_cal_method="KM")
save.image("auc_over.RData")

pdf(file = paste0("roc-",cutoff_1,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

pdf(paste0("roc-",cutoff_2,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.2y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

pdf(paste0("roc-",cutoff_3,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.3y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

pdf(paste0("roc-",cutoff_4,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.4y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

# 绘制选定模型下的auc曲线
roc_vis(all.auc.5y,
        model_name = method_select,
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)

pdf("auc_histon.pdf",w=7,h=5)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name= method_select,
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(cutoff_1,cutoff_2,cutoff_3))
dev.off()

# 5. 模型比较
# 该包还提供了和之前文献报道的预后模型比较的函数，当然只提供了胶质瘤的。
# 那如果你做的是其他癌种呢？可以通过查看函数了解是怎样的输入形式，然后就做对应的替换后就可以分析

cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,
                                             type.sig = c('Glioma','LGG','GBM'),
                                             list_input_data = list_train_vali_Data)

pdf("model_compare_dataset5_rsf.pdf",w=9,h=16)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="StepCox[forward] + RSF",
            dataset=names(list_train_vali_Data))
dev.off()
