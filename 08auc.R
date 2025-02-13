rm(list = ls())

foldpath <- paste0("D:/workdir/12stadb/08auc")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

library(survivalROC)
library(survminer)
library(survival)
library(patchwork)
library(timeROC)

# ------------风险评分----------
load("../05model/list_train_vali_Data.RData")
genes <- read.csv("../08borutax/10.boruta_geneids.csv")
genelist <- genes$x

genelist2 <- Reduce(intersect, lapply(list_train_vali_Data, 
                                     function(df) intersect(colnames(df), genelist)))


# # # # 删除小于100天生存期的样本
for (i in 1:length(list_train_vali_Data)){
  list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][list_train_vali_Data[[i]]$OS.time >90,]
  list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][list_train_vali_Data[[i]]$OS.time <365*10,]
  list_train_vali_Data[[i]] <- list_train_vali_Data[[i]][,c("ID","OS.time","OS",genelist2)]
  list_train_vali_Data[[i]] <- na.omit(list_train_vali_Data[[i]])
}

# cox-------
library(ggthemes)
library(survivalROC)
library(survminer)
library(grid)
library(ggvenn)
library(pROC)
library(xgboost)
library(randomForest)
library(caret)
library(DALEX)
library(stats)
library(e1071)
library(glmnet)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(scales)
library(survminer)
library(glmnet)
library(readxl)
library(readr)
library(tidyverse)
library(plyr)
library(glmnet)
library(survminer)
library(gridExtra)
library(patchwork)
library(forestplot)
library(survival)
library(dplyr)
library(plyr)
library(tibble)
library(tidyr)

# unicox ------------------------------------------------------------------

diff_expr_clinical <- list_train_vali_Data[[1]]
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum

covariates <- colnames_sum[-which(colnames_sum %in% c("ID","OS", "OS.time"))]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))
univ_models <- lapply(univ_formulas,
                      function(x) coxph(x, data = list_train_vali_Data[[1]]))

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
write.table(res, file = "03.uni_cox_OS_res.txt", sep = "\t", quote = F, row.names = T)

# 选择单因素P值 0.05
res_results<- na.omit(res[which(as.numeric(res$p.value) < 0.05),])
write.table(res_results, file = "04.uni_cox_OS_sig0.05.txt", sep = "\t", quote = F, row.names = T)
dim(res_results)
# [1] 14 2

res_results_all <- res_results
# 单因素cox森林图----
res_results<- res_results_all[order(res_results_all$p.value),]
# res_results <- res_results[-5,]
# res_results <- res_results[-3,]
# res_results <- res_results[-9,]
res_results_plot <- res_results
unix_res <- tidyr::separate(res_results_plot, "HR (95% CI for HR)", into = c("HR", "HR.95L", "HR.95H"), sep = " ") %>% 
  tidyr::separate("HR.95L", into = c("HR.95L", "HR.95H"), sep = "\\-") %>% 
  tibble::rownames_to_column(var="GeneName")
unix_res$HR.95L <- gsub("\\(", "", unix_res$HR.95L)
unix_res$HR.95H <- gsub("\\)", "", unix_res$HR.95H)
unix_res[, 2:ncol(unix_res)] <- as.data.frame(apply(unix_res[, 2:ncol(unix_res)], 2, as.numeric))
unix_res <- unix_res[order(unix_res$p.value),]

hz <- paste(round(unix_res$HR,4),
            "(",round(unix_res$HR.95L,3),
            "-",round(unix_res$HR.95H,3),")",sep = "")
tabletext <- cbind(c(NA,"GeneName", unix_res$GeneName),
                   c(NA,"P value", ifelse(unix_res$p.value<0.001,
                                          "< 0.001",
                                          round(unix_res$p.value,4))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
# 12

pdf(file = "13.univariate_cox_forest.pdf", height = 7, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE, rep(FALSE, 54)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,unix_res$HR),
           lower=c(NA,NA,unix_res$HR.95L), #95%置信区间下限
           upper=c(NA,NA,unix_res$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2, 3), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=1.2, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface = "bold", fontfamily="Times"),
                          summary = gpar(cex=1.2, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.2, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily="Times")),
           xlab="Hazard Ratio",
           cex.lab=2,
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()
# ---------------------------------------------------------------------
# Lasso -------------------------------------------------------------------
x_all <- subset(diff_expr_clinical, select = -c(ID, OS, OS.time))
x_all <- x_all[,rownames(res_results)]
y_all <- subset(diff_expr_clinical, select = c(OS, OS.time))

fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS),
              family = "cox")
plot(fit, xvar = "lambda",label = TRUE, las=1)

png(filename = "05.lasso_model.png", height = 450, width = 600)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "05.lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()

set.seed(99999)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),
                  nfold=10,
                  family = "cox")
plot(cvfit, las =1)

pdf(file = "06.lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()

# ------多因素----------
train_data <- list_train_vali_Data[[1]][,c("OS.time","OS",unix_res$GeneName)]
train_data <- na.omit(train_data)
cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(colnames(train_data)[3:26], collapse = "+")))
cox_more <- coxph(cox_data, data = train_data)

cox_zph <- cox.zph(cox_more)
names(cox_more$coefficients) <- rownames(cox_zph$table)[1:2]

mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
  }

multi_os_clinical <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res2 <- rbind(res) %>% as.data.frame()

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_os_clinical.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

pdf(file = "07.multi_cox.pdf", family = "Times", height = 5, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           # is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,2),TRUE,rep(FALSE,4)),
           is.summary = c(rep(FALSE,7),TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),
           # is.summary = c(rep(FALSE,16),TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           # xticks = c(0, 1, 2, 4, 6, 8), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily = "Times"),
                          xlab=gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate",
           clip = c(0,7)) # 垂直于x轴的网格线，对应每个刻度
dev.off()
# # 1.Lasso模型   ----------------
# # library(lance）
# ## 01.获取数据集 -----------------------------------------------------------  
# ### 03.LASSO
train_data <- list_train_vali_Data[[1]][,c("OS.time","OS",unix_res$GeneName)]
x_all <- subset(train_data, select = -c(OS.time, OS))
y_all <- subset(train_data, select = c(OS.time, OS))


##  04.拟合模型 ----------------------------------------------------------------
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS),
              family = "cox")


## 05.交叉验证 -----------------------------------------------------------------
set.seed(99999)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=10,
                  family = "cox")
##画图

coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1] 0.02545705

# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
length(coef.min)
# 14
coef.min
# GHR     0.25419692
# ASPA    .         
# BGN     0.05611489
# ADH1B   .         
# FAP     0.03458269
# CHRDL1  .         
# COL3A1  .         
# TIAM1   0.11867804
# EFS     .         
# CADM3   .         
# THBS2   .         
# CLSPN   .         
# KHK    -0.15566338
# STIL    .         

# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
# "GHR"   "BGN"   "FAP"   "TIAM1" "KHK"

saveRDS(lasso_geneids, file= "lasso_geneids.rds")

saveRDS(coef.min, file="coef.min.rds")


lasso_geneids <- data.frame(symbol=lasso_geneids)
write.csv(lasso_geneids, "01.lasso_genes.csv",quote = F,row.names = F)

x <- coef(fit)
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1]
# extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1]
# compute L1 norm

#图片美化
head(tmp)
col <- rep(c('#ff9898','#dedb8e','#99e0ab','#D94F04','#007172',"#E9967A","#FA8072",
             '#025259','#c49d93','#aec6e8','#F2C6C2','#86A69D',"black",  "#FFA07A"),2)
pdf("02.lasso_model_name.pdf",height = 5, width = 7,family='Times')
ggplot(tmp,aes(log(lambda),value,color = coef)) +
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=1) +
  xlab("Log(Lambda)") +
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+
  scale_color_manual(values= col)+

  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  annotate('text',x = -4.3,y=0.18,
           label=paste0('Optimal Lambda =', cvfit$lambda.min),
           color='black')+
  guides(col=guide_legend(ncol = 2))
dev.off()
# 
# #--------------------------- 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = cvfit$lambda.min)  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1] 0.02545705
df.coef = cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame
write.table(df.coef, "Lasso_Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)  #6个

save.image("lasso_over.Rdata")
load("lasso_over.Rdata")

# ----------timeROC----------------
# ---------融合----------
rs <- list_train_vali_Data

coxsig <- lasso_geneids$symbol

rs <- lapply(rs, 
             function (x) x[,c("ID","OS.time","OS",coxsig)])

cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(coxsig, collapse = "+")))

cox_model_all <- coxph(cox_data,data = list_train_vali_Data[[1]])
coef <- cox_model_all$coefficients
write.csv(cox_model_all$coefficients,file = "coef.csv")

risk_scores <- predict(cox_model_all, newdata = list_train_vali_Data[[1]], type = "risk")
rs[[1]]$RS <- risk_scores

for (i in 2:6){
  cox_model_all <- coxph(cox_data,data = list_train_vali_Data[[i]])
  risk_scores <- predict(cox_model_all, newdata = list_train_vali_Data[[i]], type = "risk")
  rs[[i]]$RS <- risk_scores
}

save(rs, file = "riskscore.RData")
# ---------------------------
load("riskscore.RData")
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
  ROC_merge[[i]] <- timeROC(T=list_train_vali_Data[[i]]$OS.time,   
                            delta=list_train_vali_Data[[i]]$OS,   
                            marker=rs[[i]]$RS,   
                            cause=1,                #阳性结局指标数值
                            weighting="marginal",   #计算方法，默认为marginal
                            times=dataspan,       #时间点，选取1年，3年和5年的生存率
                            iid=TRUE)
}

save(ROC_merge,file = "roc_merge.RData")
# ----------绘图----------------
# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

col=c("#EB4B17", "#2775AB", '#4C8045',"#D8D155","#E0367A","#91612D","#E0867B","#35112D")
# --------时间推移----------
pdf("timeauc.pdf",w=8,h=5)
par(mai=c(1,1,1,3))
plot(dataspan, ROC_merge[[1]]$AUC, lwd=3, type = "l", col = col[1], 
     ylim = c(0.5, 0.9),
     xlim = c(360,1800),
     xlab = "Time (days)", ylab = "AUC", main = paste0("AUC"), bty = "l", xaxt = "n")
axis(1, dataspan)
for (i in c(1:6)){
  lines(dataspan, lwd=3, ROC_merge[[i]]$AUC, col = col[i])
}

legendtxt <- c()
for (i in 1:6){
  legendtxt <- c(legendtxt,paste0(names(rs[i])," AUC=",round(mean(ROC_merge[[i]][["AUC"]],na.rm=T),3),
                                  " CI(",round(mean(confint(ROC_merge[[i]], level = 0.95)$CI_AUC[,1])/100,2),
                                  "-",round(mean(confint(ROC_merge[[i]], level = 0.95)$CI_AUC[,2])/100,2),")\n"))
}
# legendtxt <- c()
# for (i in c(1:6)){
#   legendtxt <- c(legendtxt,paste0(names(rs[i])," AUC=",round(median(ROC_merge[[i]][["AUC"]],na.rm=T),3)))
# }
legend(x= 2000,y= 0.8,xpd = TRUE,
       legendtxt,
       col= col,
       lty=1, lwd=3, cex = 0.8)

dev.off()
