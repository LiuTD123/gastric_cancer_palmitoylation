rm(list = ls())

foldpath <- "D:/workdir/12stadb/10diag"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(glmnet)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(xgboost)
library(randomForest)
library(caret)
library(DALEX)
library(stats)
library(e1071)
library(readxl)
library(gplots)
library(tableone)

# ----数据集整理---------
# ----------stad-----------
data_train <- readRDS("../00data/tcga_fpkm.rds")
tcga_group <-readRDS("..\\00data\\tcga_group.rds")
genelist <- readRDS("../08auc/lasso_geneids.rds")

data_train <- as.data.frame(t(data_train[genelist,]))
data_train$sample <- rownames(data_train)
data_train <- merge(tcga_group,data_train, by = "sample")
table(data_train$group)
# -----geo------
dir ="../00geodata"
samples=list.files(dir,pattern = "^GSE.*RData$")

# ------------------
gseid <- "GSE13861"
load(paste0("../00geodata/dat_",gseid,".RData"))
exp <- as.data.frame(t(exp_symbol[genelist,]))
exp$sample <- rownames(exp)

load(paste0("../00geodata/",gseid,"_cliall.RData"))
table(cli$characteristics_ch1)
cli <- cli[,c("geo_accession","characteristics_ch1")]
colnames(cli) <- c("sample","group")

cli <- cli[cli$group %in% c("gastric adenocarcinoma","normal surrounding gastric tissue"),]
cli$group <- ifelse(cli$group == "gastric adenocarcinoma","Tumor","Normal")

table(cli$group)
GSE13861 <- merge(cli,exp,by = "sample")
table(GSE13861$group)
# ------------
gseid <- "GSE26942"

load(paste0("../00geodata/dat_",gseid,".RData"))
exp <- as.data.frame(t(exp_symbol[genelist,]))
exp$sample <- rownames(exp)

load(paste0("../00geodata/",gseid,"_cliall.RData"))
table(cli$source_name_ch1)
cli <- cli[,c("geo_accession","source_name_ch1")]
colnames(cli) <- c("sample","group")

cli$group <- ifelse(cli$group == "Gastric tumor tissue","Tumor","Normal")

table(cli$group)
GSE26942 <- merge(cli,exp,by = "sample")
table(GSE26942$group)
# ------------------
gseid <- "GSE29272"
load(paste0("../00geodata/dat_",gseid,".RData"))
exp <- as.data.frame(t(exp_symbol[genelist,]))
exp$sample <- rownames(exp)

load(paste0("../00geodata/",gseid,"_cliall.RData"))
table(cli$characteristics_ch1)
cli <- cli[,c("geo_accession","characteristics_ch1")]
colnames(cli) <- c("sample","group")

cli$group <- ifelse(cli$group == "tissue: adjacent tissue normal gastric glands","Normal","Tumor")

table(cli$group)
GSE29272 <- merge(cli,exp,by = "sample")
table(GSE29272$group)
# ---------------------------------
gseid <- "GSE66229"

load(paste0("../00geodata/dat_",gseid,".RData"))
exp <- as.data.frame(t(exp_symbol[genelist,]))
exp$sample <- rownames(exp)

load(paste0("../00geodata/",gseid,"_cliall.RData"))
table(cli$characteristics_ch1)
cli <- cli[,c("geo_accession","characteristics_ch1")]
colnames(cli) <- c("sample","group")

cli$group <- ifelse(cli$group == "tissue: Gastric tumor","Tumor","Normal")

table(cli$group)
GSE66229 <- merge(cli,exp,by = "sample")
table(GSE66229$group)
# -------stad整理----------
stad_diagnosticdata <- list("Train" = data_train,
                            "GSE13861" = GSE13861,
                            "GSE26942" = GSE26942)

save(stad_diagnosticdata, file = "stad_diagnosticdata.RData")
table(stad_diagnosticdata[[1]]$group)

# -------------
load("./stad_diagnosticdata.RData")

for (i in 1:3){
  stad_diagnosticdata[[i]]$group <- ifelse(stad_diagnosticdata[[i]]$group == "Tumor",1,0)
  stad_diagnosticdata[[i]] <- stad_diagnosticdata[[i]][,-1]
}

library(dplyr)

data <- stad_diagnosticdata[[1]]

# ----------------------------------------------------------------
# load("data.Rdata")
# lasso <- read.csv("../../03lasso/model1/03.lasso.csv")
# input <- input[,c("TotalIgE",lasso$symbol)]
# testset <- testset[,c("TotalIgE",lasso$symbol)]
input <- data

ctrl <- trainControl(method="cv", number=500)
classif_knn <- train(group~., data = input,
                     method = "knn", trControl=ctrl, tuneLength=10)

classif_pls <- train(group~., data = input,
                     method = "pls",validation = 'CV')

# 支持向量机模型
classif_svm <- train(group~., data = input, 
                     method = "svmRadial")

# 随机森林模型
classif_rf <- train(group~., data = input, 
                    method = "rf",
                    ntree = 300)
# 广义线性模型
classif_glm <- train(group~., data = input,
                     method = 'glm')


# # 极限梯度提升模型
x = model.matrix(group~.,input)
model_martix_train<-model.matrix(group~., input)

data_train <- xgb.DMatrix(x, label =as.numeric(as.factor(input$group)))

params <- list(
  objective = "reg:squarederror"
)
# 
classif_xgboost <- xgb.train(params, data_train, nrounds = 100)
# 
# save(classif_xgboost,classif_glm,classif_rf,classif_svm,file="fourModel_classif.RData")
explainer_knn<-explain(classif_knn,label = "KNN",
                       data = input,
                       y = input$group)

explainer_pls<-explain(classif_pls,label = "PLS",
                       data = input,
                       y = input$group)

explainer_svm<-explain(classif_svm,label = "SVM",
                       data = input,
                       y = input$group)

explainer_rf<-explain(classif_rf,label = "RF",
                      data = input,
                      y = input$group)

explainer_glm<-explain(classif_glm,label = "GLM",
                       data = input,
                       y = input$group)

# -------------------------------xgboost
predict_logit <- function(model,x){
  raw_x <-predict(model,x)
  exp(raw_x)/(1+exp(raw_x))
}

logit <- function(x){
  exp(x)/(1+exp(x))
}

explainer_xgboost<-explain(classif_xgboost,
                           label = "xgboost",
                           data = x,
                           y = as.numeric(input$group),
                           predict_function = predict_logit,
                           link = logit
)

# save(explainer_xgboost,explainer_glm,explainer_rf,explainer_svm,file="fourModel_explainer.RData")

# model performance
# per_knn<-model_performance(explainer_knn)
library(pROC)
per_pls<-model_performance(explainer_pls, measure = "auc")
per_knn<-model_performance(explainer_knn, measure = "auc")
per_svm<-model_performance(explainer_svm, measure = "auc")
per_rf<-model_performance(explainer_rf, measure = "auc")
per_glm<-model_performance(explainer_glm, measure = "auc")
per_xgboost<-model_performance(explainer_xgboost, measure = "auc")

auc_pls <- auc(per_pls$residuals$observed,per_pls$residuals$predicted)
auc_glm <- auc(per_glm$residuals$observed,per_glm$residuals$predicted)
auc_rf <- auc(per_rf$residuals$observed,per_rf$residuals$predicted)
auc_svm<- auc(per_svm$residuals$observed,per_svm$residuals$predicted)
auc_knn <- auc(per_knn$residuals$observed,per_knn$residuals$predicted)
auc_xgboost <- auc(per_xgboost$residuals$observed,per_xgboost$residuals$predicted)

# ROC曲线
psvm <- plot(per_svm, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_svm, 3)))

pknn <- plot(per_knn, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_knn, 3)))

prf <- plot(per_rf, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_rf, 3)))

pglm <- plot(per_glm, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_glm, 3)))

ppls <- plot(per_pls, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_pls, 3)))

pxgboost <- plot(per_xgboost, geom = "roc")+
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_xgboost, 3)))

library(cowplot)
pdf("Train_ML_all_ROC2.pdf",h = 8, w=14)
plot_grid(pknn, pglm, ppls, prf, psvm, pxgboost, ncol=3)
dev.off()

# ------------ROC曲线放大一个图------
plotlist <- list(pknn, pglm, ppls, prf, psvm, pxgboost)

method <- c("KNN","GLM","PLS","RF","SVM","XGBoost")
col=c("#EB4B17", "#2775AB", "#91612D",'#4C8045',"#D8D155","#E0867B","#35112D","#E0367A")

auclist <- c(auc_knn,auc_glm,auc_pls,auc_rf,auc_svm,auc_xgboost)

pdf("Train_diagnostic_auc.pdf",w = 6, h = 6)
for (i in 1:length(plotlist)){
  if (i == 1){
    plot(x = plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr, lwd=3, type = "l", col = col[i],
         xlab = "False positive rate", ylab = "True positive rate",
         main = paste0("Gastric Cancer Diagnostic"),
         bty = "l", xaxt = "n")
  } else {
    lines(plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr,
          lwd=3, col = col[i])
  }
}
lines(x=c(0,1),
      y=c(0,1),
      lwd=3, col = "gray")
legend("bottomright",
       legend=c(paste0(method[1]," AUC=",round(auclist[1],3),"\n"),
                paste0(method[2]," AUC=",round(auclist[2],3),"\n"),
                paste0(method[3]," AUC=",round(auclist[3],3),"\n"),
                paste0(method[4]," AUC=",round(auclist[4],3),"\n"),
                paste0(method[5]," AUC=",round(auclist[5],3),"\n"),
                paste0(method[6]," AUC=",round(auclist[6],3)
                )),
       col=col,
       lwd=2,
       title="Curve")
dev.off()

# ===============================
# ------------------------------------------------------

# plot model performance
pdf("04.machine_learning_residuals_line2.pdf",height = 4,width = 6)
plot(per_knn,per_glm,per_pls,per_rf,per_svm,per_xgboost)
dev.off()

pdf("04.machine_learning_residuals_box2.pdf",height = 4,width = 6)
plot(per_knn,per_glm,per_pls,per_rf,per_svm,per_xgboost,geom = "boxplot")
dev.off()

# importance
importance_knn<-variable_importance(
  explainer_knn,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_pls<-variable_importance(
  explainer_pls,
  loss_function = loss_root_mean_square
)
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_xgboost<-variable_importance(
  explainer_xgboost,
  loss_function = loss_root_mean_square
)
write.csv(importance_knn,file = "importance_knn.csv")
write.csv(importance_glm,file = "importance_glm.csv")
write.csv(importance_pls,file = "importance_pls.csv")
write.csv(importance_rf,file = "importance_rf.csv")
write.csv(importance_xgboost,file = "importance_xgboost.csv")
write.csv(importance_svm,file = "importance_svm.csv")

pdf("04.machine_learning_importance2.pdf",height = 11,width = 22)
# par(mfrow = c(2,3), xpd=TRUE)
# plot(importance_knn,importance_glm,importance_pls,importance_rf,importance_svm,importance_xgboost)
knn <- plot(importance_knn)
glm <- plot(importance_glm)
pls <- plot(importance_pls)
rf <- plot(importance_rf)
svm <- plot(importance_svm)
xgboost <- plot(importance_xgboost)
# par(2,2,2,2)
plot_grid(knn, glm, pls, rf, svm, xgboost, ncol=3)
dev.off()

# save.image("modeltrained.RData")

# ---------------验证-----------------
# ====================================
for (i in 2:3){
  datasets <- names(stad_diagnosticdata[i])
  
  testset <- stad_diagnosticdata[[i]]
  
  explainer_knn<-explain(classif_knn,label = "KNN",
                         data = testset,
                         y = testset$group)
  
  explainer_pls<-explain(classif_pls,label = "PLS",
                         data = testset,
                         y = testset$group)
  
  explainer_svm<-explain(classif_svm,label = "SVM",
                         data = testset,
                         y = testset$group)
  
  explainer_rf<-explain(classif_rf,label = "RF",
                        data = testset,
                        y = testset$group)
  
  explainer_glm<-explain(classif_glm,label = "GLM",
                         data = testset,
                         y = testset$group)
  
  # -------------------------------xgboost
  x = model.matrix(group~.,testset)

  params <- list(
    objective = "reg:squarederror"
  )
  # 
  # classif_xgboost <- xgb.train(params, data_train, nrounds = 100)
  predict_logit <- function(model,x){
    raw_x <-predict(model,x)
    exp(raw_x)/(1+exp(raw_x))
  }
  
  logit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  explainer_xgboost<-explain(classif_xgboost,
                             label = "xgboost",
                             data = x,
                             y = as.numeric(testset$group),
                             predict_function = predict_logit,
                             link = logit
  )
  
  
  library(pROC)
  per_pls<-model_performance(explainer_pls, measure = "auc")
  per_knn<-model_performance(explainer_knn, measure = "auc")
  per_svm<-model_performance(explainer_svm, measure = "auc")
  per_rf<-model_performance(explainer_rf, measure = "auc")
  per_glm<-model_performance(explainer_glm, measure = "auc")
  per_xgboost<-model_performance(explainer_xgboost, measure = "auc")
  
  auc_pls <- auc(per_pls$residuals$observed,per_pls$residuals$predicted)
  auc_glm <- auc(per_glm$residuals$observed,per_glm$residuals$predicted)
  auc_rf <- auc(per_rf$residuals$observed,per_rf$residuals$predicted)
  auc_svm<- auc(per_svm$residuals$observed,per_svm$residuals$predicted)
  auc_knn <- auc(per_knn$residuals$observed,per_knn$residuals$predicted)
  auc_xgboost <- auc(per_xgboost$residuals$observed,per_xgboost$residuals$predicted)
  
  # ROC曲线
  psvm <- plot(per_svm, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_svm, 3)))
  
  pknn <- plot(per_knn, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_knn, 3)))
  
  prf <- plot(per_rf, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_rf, 3)))
  
  pglm <- plot(per_glm, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_glm, 3)))
  
  ppls <- plot(per_pls, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_pls, 3)))
  
  pxgboost <- plot(per_xgboost, geom = "roc")+
    annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_xgboost, 3)))
  
  # ---------画图-----------
  plotlist <- list(pknn, pglm, ppls, prf, psvm, pxgboost)
  method <- c("KNN","GLM","PLS","RF","SVM","XGBoost")
  col=c("#EB4B17", "#2775AB", "#91612D",'#4C8045',"#D8D155","#E0867B","#35112D","#E0367A")
  
  auclist <- c(auc_knn,auc_glm,auc_pls,auc_rf,auc_svm,auc_xgboost)
  
  pdf(paste0(datasets,"_diagnostic_auc.pdf"),w = 6, h = 6)
  for (i in 1:length(plotlist)){
    if (i == 1){
      plot(x = plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr, lwd=3, type = "l", col = col[i],
           xlab = "False positive rate", ylab = "True positive rate",
           main = paste0(datasets," Diagnostic"),
           bty = "l", xaxt = "n")
    } else {
      lines(plotlist[[i]]$data$fpr, y = plotlist[[i]]$data$tpr,
            lwd=3, col = col[i])
    }
  }
  lines(x=c(0,1),
        y=c(0,1),
        lwd=3, col = "gray")
  legend("bottomright",
         legend=c(paste0(method[1]," AUC=",round(auclist[1],3),"\n"),
                  paste0(method[2]," AUC=",round(auclist[2],3),"\n"),
                  paste0(method[3]," AUC=",round(auclist[3],3),"\n"),
                  paste0(method[4]," AUC=",round(auclist[4],3),"\n"),
                  paste0(method[5]," AUC=",round(auclist[5],3),"\n"),
                  paste0(method[6]," AUC=",round(auclist[6],3)
                  )),
         col=col,
         lwd=2,
         title="Curve")
  dev.off()
}
