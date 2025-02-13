options(timeout = Inf)
rm(list = ls())
foldpath <- paste0("D:/workdir/12stadb/07cox")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

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

fpkm <- readRDS("../00data/tcga_fpkm.rds")
group <- readRDS("../00data/tcga_group.rds")
# keygene <- read.csv("/data/nas1/liuhouyan/01-BJTC-577-2/04_candidate_gene/DEGs-PRGs.csv")
# load("/data/nas1/liuhouyan/01-BJTC-577-2/15_wgcna/cluster_and_expression.RData")
keygene <- read.csv("../08borutax/10.boruta_geneids.csv")

# ppi_top20 <- read.csv("/data/nas1/liuhouyan/01-BJTC-577-2/04_candidate_gene/ppi/10.Venn_Con.csv")

survival_dat <- readRDS("../00data/survival_dat.rds")
colnames(survival_dat) <- c("id","OS.time","OS")
# survival_dat$id <- substr(survival_dat$id, 0, 15)

fpkm_1 <- na.omit(fpkm[keygene$x,]) %>%t %>% as.data.frame() %>% rownames_to_column(var = "id")
# fpkm_1 <- fpkm[ppi_top20$x,] %>%t %>% as.data.frame() %>% rownames_to_column(var = "id")
fpkm_1$id <- gsub('\\.','\\-',fpkm_1$id)

fpkm_tem <- merge(fpkm_1,survival_dat,by = "id")

rownames(fpkm_tem) <- fpkm_tem$id
fpkm_tem <- fpkm_tem[,-1]
# fpkm_tem <- fpkm_tem[,-11]
df_merge <- fpkm_tem

# unicox ------------------------------------------------------------------

diff_expr_clinical <- df_merge
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum

covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))
univ_models <- lapply(univ_formulas,
                      function(x) coxph(x, data = diff_expr_clinical))

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
# [1] 5 2

res_results_all <- res_results
res_results1 <- res_results_all[! rownames(res_results_all)%in%c("CHST8","HS6ST3","EXTL1"),]
res_results2 <- res_results_all[c("CHST8","HS6ST3","EXTL1"),]
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

pdf(file = "13.univariate_cox_forest.pdf", height = 9, width = 10, onefile = F)
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
           xticks = c(0, 1, 2), #横坐标刻度
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

save.image("unicox.Rdata")
# ---------------------------------------------------------------------
# # PH检验

for (gene in rownames(res_results)){
  # gene <- cox.zph(univ_models$gene)
  # print(gene)
  print(paste0(gene,"<-","cox.zph(univ_models$",gene,")"))
  print(paste0(gene, "<- ggcoxzph(",gene, ",caption = '",gene,"')"))
  print(paste0("ggsave(","'./ph/ggcoxzph_",gene,".pdf'",",arrangeGrob(grobs = ",gene,"))"))
  print(paste0("ggsave(","'./ph/ggcoxzph_",gene,".png'",",arrangeGrob(grobs = ",gene,"))"))
}

CCNF<-cox.zph(univ_models$CCNF)
CCNF<- ggcoxzph(CCNF,caption = 'CCNF')
ggsave('./ph/ggcoxzph_CCNF.pdf',arrangeGrob(grobs = CCNF))
ggsave('./ph/ggcoxzph_CCNF.png',arrangeGrob(grobs = CCNF))
NACAD<-cox.zph(univ_models$NACAD)
NACAD<- ggcoxzph(NACAD,caption = 'NACAD')
ggsave('./ph/ggcoxzph_NACAD.pdf',arrangeGrob(grobs = NACAD))
ggsave('./ph/ggcoxzph_NACAD.png',arrangeGrob(grobs = NACAD))
CILP<-cox.zph(univ_models$CILP)
CILP<- ggcoxzph(CILP,caption = 'CILP')
ggsave('./ph/ggcoxzph_CILP.pdf',arrangeGrob(grobs = CILP))
ggsave('./ph/ggcoxzph_CILP.png',arrangeGrob(grobs = CILP))
ZNF316<-cox.zph(univ_models$ZNF316)
ZNF316<- ggcoxzph(ZNF316,caption = 'ZNF316')
ggsave('./ph/ggcoxzph_ZNF316.pdf',arrangeGrob(grobs = ZNF316))
ggsave('./ph/ggcoxzph_ZNF316.png',arrangeGrob(grobs = ZNF316))
SNAPC4<-cox.zph(univ_models$SNAPC4)
SNAPC4<- ggcoxzph(SNAPC4,caption = 'SNAPC4')
ggsave('./ph/ggcoxzph_SNAPC4.pdf',arrangeGrob(grobs = SNAPC4))
ggsave('./ph/ggcoxzph_SNAPC4.png',arrangeGrob(grobs = SNAPC4))
SIGLEC1<-cox.zph(univ_models$SIGLEC1)
SIGLEC1<- ggcoxzph(SIGLEC1,caption = 'SIGLEC1')
ggsave('./ph/ggcoxzph_SIGLEC1.pdf',arrangeGrob(grobs = SIGLEC1))
ggsave('./ph/ggcoxzph_SIGLEC1.png',arrangeGrob(grobs = SIGLEC1))
CERCAM<-cox.zph(univ_models$CERCAM)
CERCAM<- ggcoxzph(CERCAM,caption = 'CERCAM')
ggsave('./ph/ggcoxzph_CERCAM.pdf',arrangeGrob(grobs = CERCAM))
ggsave('./ph/ggcoxzph_CERCAM.png',arrangeGrob(grobs = CERCAM))
ZMIZ2<-cox.zph(univ_models$ZMIZ2)
ZMIZ2<- ggcoxzph(ZMIZ2,caption = 'ZMIZ2')
ggsave('./ph/ggcoxzph_ZMIZ2.pdf',arrangeGrob(grobs = ZMIZ2))
ggsave('./ph/ggcoxzph_ZMIZ2.png',arrangeGrob(grobs = ZMIZ2))
PRR11<-cox.zph(univ_models$PRR11)
PRR11<- ggcoxzph(PRR11,caption = 'PRR11')
ggsave('./ph/ggcoxzph_PRR11.pdf',arrangeGrob(grobs = PRR11))
ggsave('./ph/ggcoxzph_PRR11.png',arrangeGrob(grobs = PRR11))

a <- paste(rownames(res_results),collapse = ",")

ggsave("ggcoxzph_all.pdf",w = 14,h = 12,
       arrangeGrob(grobs = c(CCNF,NACAD,CILP,ZNF316,SNAPC4,SIGLEC1,CERCAM,ZMIZ2,PRR11)))

# FOS<-cox.zph(univ_models$FOS)
# FOS<- ggcoxzph(FOS,title = "FOS")
# ggsave("ggcoxzph_FOS.pdf",arrangeGrob(grobs = FOS))
# ggsave("ggcoxzph_FOS.png",arrangeGrob(grobs = FOS))
# 
# STK40<-cox.zph(univ_models$STK40)
# STK40<- ggcoxzph(STK40,title = "STK40")
# ggsave("ggcoxzph_STK40.pdf",arrangeGrob(grobs = STK40))
# ggsave("ggcoxzph_STK40.png",arrangeGrob(grobs = STK40))

# ----------------------------
# PH检验
x <- df_merge[,c(101,100)]
uniSigExp = df_merge
uniSigExp = uniSigExp[,rownames(res_results)]
uniSigExp <- merge(uniSigExp,x, by = "row.names")
rownames(uniSigExp) <- uniSigExp$Row.names
uniSigExp <- uniSigExp[,-1]
dat <- uniSigExp
outPH=data.frame()
for(i in colnames(dat[,1:11])){
  cox <- coxph(Surv(OS.time, OS) ~ dat[,i], data = dat)
  test.ph <- cox.zph(cox)
  #coxP=test.ph$coefficients[,"Pr(>|z|)"]
  outPH=rbind(outPH,
              cbind(id=i,
                    p=test.ph$table[1,"p"])
  )
}

sigPH=outPH[as.numeric(as.vector(outPH$p))>0.05,]  # 46
write.table(sigPH,file="02.PH.Sig.txt",sep="\t",row.names=F,quote=F)
# -----------------------------------------

# Lasso -------------------------------------------------------------------

# diff_expr_clinical <- read.table("./02.TCGA-BLCA_logtpm_cli.xls", header = T)
# res_results_0.05 <- read.table("./04.uni_cox_OS_sig0.05.txt", header = T, sep = "\t", check.names = F)
x_all <- subset(diff_expr_clinical, select = -c(OS, OS.time))
# x_all <- x_all[,rownames(res_results)]
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

set.seed(555555555)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),
                  nfold=10,
                  family = "cox")
plot(cvfit, las =1)

png(filename = "06.lasso_verify.png", height = 450, width = 600)
plot(cvfit, las =1)
dev.off()
pdf(file = "06.lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()

# 1.Lasso模型   ----------------
# library(lance）

fpkm <- readRDS("..\\00data\\tcga_fpkm.rds")
group <- readRDS("../00data/tcga_group.rds")
## 01.获取数据集 -----------------------------------------------------------  
dat <- fpkm
# dat <- log2(dat+1)
keygene <- read.csv("../08borutax//10.boruta_geneids.csv")

survival <- readRDS("../00data/survival_dat.rds")

gene <- res_results


## 02.合并生存数据
survival_dat <- t(dat[rownames(gene),colnames(dat) %in% survival$sample])
train_dat <- survival_dat %>% data.frame()
train_dat$sample <- rownames(train_dat)
train_dat <- merge(survival,train_dat,by='sample')
rownames(train_dat) <- train_dat$sample
train_dat<-train_dat[,c(-1)]
colnames(train_dat)


### 03.LASSO

train_data <- train_dat
x_all <- subset(train_data, select = -c(fustat, futime))
y_all <- subset(train_data, select = c(fustat, futime))


##  04.拟合模型 ----------------------------------------------------------------
fit <- glmnet(as.matrix(x_all), Surv(y_all$futime,y_all$fustat), 
              family = "cox") 


## 05.交叉验证 -----------------------------------------------------------------
set.seed(555555)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$futime,y_all$fustat),nfold=10,
                  family = "cox") 
##画图

coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1]  0.03441655

# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
length(coef.min)
# 19
# coef.min
# 19 x 1 sparse Matrix of class "dgCMatrix"
# 1
# GHR       0.218422910
# FNDC1     0.058541256
# CHRDL1    .          
# ADH1B     .          
# ASPA      .          
# BGN       .          
# COL3A1    .          
# COL3A1.1  .          
# EFS       .          
# FAP       .          
# LGI4      .          
# CADM3     .          
# STIL      .          
# THBS2     .          
# KIF18B    .          
# KHK      -0.101646739
# ANKRD35   .          
# TIAM1     0.007329178
# KCNA5     .  

# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids

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

pdf("02.lasso_model_name.pdf",height = 5, width = 8,family='Times')
ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Log(Lambda)") + 
  ylab('Coefficients')+ 
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values= rep(c('#ff9898','#dedb8e','#99e0ab','#D94F04','#007172',"#E9967A","#FA8072",
                               '#025259','#c49d93','#aec6e8','#F2C6C2','#86A69D',"black",  "#FFA07A"),2))+
  
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

#--------------------------- 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = cvfit$lambda.min)  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1] 0.03441655
df.coef = cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame
write.table(df.coef, "Lasso_Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)  #6个

save.image("lasso_over.Rdata")
load("lasso_over.Rdata")

