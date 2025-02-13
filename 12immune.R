rm(list = ls())

foldpath <- "D:/workdir/12stadb/12immune"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
#####免疫浸润分析
library(psych)
library(ggcorrplot)
library(corrplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(reshape2)
library(magrittr)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(immunedeconv)
library(tidyr)
library(IOBR)
library(reshape2)
library(rstatix)
library(patchwork)
library(CIBERSORT)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos = rforge, dependencies = TRUE)
# install.packages("D:/workdir/Rpackage/estimate_1.0.13.tar.gz")

# 输入参考
# OvarianCancerExpr <- system.file("extdata", "sample_input.txt", package = "estimate")
# exampleinput <- read.table(OvarianCancerExpr)[1:4,1:4]
# s516      s518      s519      s520
# C9orf152  4.881540  4.575656  3.739469  3.695996
# ELMO2     7.298054  7.555440  7.533202  7.382355
# CREB3L1   5.569164  5.700406  5.959730  5.770007
# RPS11    13.389937 13.848820 13.642862 13.654622

library(estimate)
data_train <- readRDS("../00data/tcga_fpkm.rds")

head(data_train)[1:4,1:4]
# data_train$gene <- rownames(data_train)
# data_train <- data_train[,c(449,1:448)]

write.table(data_train,file = "comexp.txt",sep = "\t", quote =F)
input_file_dir <- './comexp.txt'
output_file_dir <- './genes.gct'
output_estimate <- './estimate_score.gct'

which(duplicated(data_train$symbol))

data_train$symbol <- rownames(data_train)
data_train <- data_train[,c(449,1:448)]
data_train <- na.omit(data_train)

# 根据表达谱生成gct
filterCommonGenes(input.f = input_file_dir, 
                  output.f = "genes.gct", 
                  id = "GeneSymbol")

estimateScore(input.ds = output_file_dir,
              output.ds = "estimate_score.gct", 
              platform = "affymetrix")

scores <- read.table("estimate_score.gct",skip = 2,header = T)
rownames(scores) <- scores[,1]
scores <- t(scores[,4:ncol(scores)])

# plotPurity函数可以根据保存好的文件来挑选对应的样本进行可视化：
# plotPurity(scores = "estimate_score.gct", samples = "GSM1192715", platform = "affymetrix")

lassogene <- colnames(rs[[1]])[c(4:8)]

# -------------------------
# =========================
load("../08auc/riskscore.RData")
risk <- rs[[1]]
cut <- survminer::surv_cutpoint(risk, #数据集
                                minprop = 0.25,
                                time = "OS.time", #生存时间
                                event = "OS", #生存状态
                                variables = "RS"  #需要计算的数据列名
)
cut

risk$group <- ifelse(risk$RS > summary(cut)[1,1],'High','Low' )

TrainGroup <- group <- risk[,c(1,10),drop = F]
save(risk, file = "../08auc/riskgroup.RData")

load("../08auc/riskgroup.RData")
# --------------------
data_train <- as.data.frame(t(data_train[,-1]))
exprSet <- data.frame(scores,
                      GHR = data_train[,"GHR"],
                      BGN = data_train[,"BGN"],
                      FAP = data_train[,"FAP"],
                      TIAM1 = data_train[,"TIAM1"],
                      KHK = data_train[,"KHK"]
)

library(ggstatsplot)

# ------------StromalScore----------------
p1 <- ggscatterstats(exprSet,x = StromalScore,y = GHR,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p2 <- ggscatterstats(exprSet,x = StromalScore,y = BGN,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p3 <- ggscatterstats(exprSet,x = StromalScore,y = FAP,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p4 <- ggscatterstats(exprSet,x = StromalScore,y = TIAM1,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p5 <- ggscatterstats(exprSet,x = StromalScore,y = KHK,marginal = T,
                     type = "nonparametric",title = "log2TPM")
# ------------ImmuneScore----------------
p1 <- ggscatterstats(exprSet,x = ImmuneScore,y = GHR,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p2 <- ggscatterstats(exprSet,x = ImmuneScore,y = BGN,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p3 <- ggscatterstats(exprSet,x = ImmuneScore,y = FAP,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p4 <- ggscatterstats(exprSet,x = ImmuneScore,y = TIAM1,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p5 <- ggscatterstats(exprSet,x = ImmuneScore,y = KHK,marginal = T,
                     type = "nonparametric",title = "log2TPM")

# ----------ESTIMATEScore-----------------
p1 <- ggscatterstats(exprSet,x = ESTIMATEScore,y = GHR,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p2 <- ggscatterstats(exprSet,x = ESTIMATEScore,y = BGN,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p3 <- ggscatterstats(exprSet,x = ESTIMATEScore,y = FAP,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p4 <- ggscatterstats(exprSet,x = ESTIMATEScore,y = TIAM1,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p5 <- ggscatterstats(exprSet,x = ESTIMATEScore,y = KHK,marginal = T,
                     type = "nonparametric",title = "log2TPM")

# ----------TumorPurity-----------------
p1 <- ggscatterstats(exprSet,x = TumorPurity,y = GHR,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p2 <- ggscatterstats(exprSet,x = TumorPurity,y = BGN,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p3 <- ggscatterstats(exprSet,x = TumorPurity,y = FAP,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p4 <- ggscatterstats(exprSet,x = TumorPurity,y = TIAM1,marginal = T,
                     type = "nonparametric",title = "log2TPM")

p5 <- ggscatterstats(exprSet,x = TumorPurity,y = KHK,marginal = T,
                     type = "nonparametric",title = "log2TPM")

pdf("Stromalscore.pdf",w = 23,h=3.5)
pdf("Immunescore.pdf",w = 23,h=3.5)
pdf("ESTIMATEScore.pdf",w = 23,h=3.5)
pdf("TumorPurity.pdf",w = 23,h=3.5)

p1+p2+p3+p4+p5 + plot_layout(ncol = 5)
dev.off()

# # ----------------------------cibersort评估细胞丰度-免疫细胞比例--------------------------
a = data_train

k = !duplicated(rownames(a));table(k)
exp = a[k,]

rownames(exp) = unique(rownames(a))
colnames(exp) = str_remove(colnames(exp),"TPM")

exp[1:4,1:4]
exp2 = as.data.frame(exp)
# exp2 = rownames_to_column(exp2)
exp2 <- as.data.frame(t(exp))
write.table(exp2,file = "exp.txt",row.names = F,quote = F,sep = "\t")

# 内置数据LM22.txt，记录了22种免疫细胞的基因表达特征数据。
lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")

exp2 <- na.omit(exp2)
TME.results = cibersort(lm22f,
                        "comexp.txt",
                        QN = F
)

TME.results <- TME.results[-1,]
save(TME.results,file = "TME.result.RData")
TME.results[1:4,1:4]
#             B cells naive B cells memory Plasma cells T cells CD8
# GSM404005    0.00000000    0.118715823    0.5064640 0.002298742
# GSM404006    0.02075989    0.000000000    0.7821873 0.007487980
# GSM404007    0.00000000    0.224421334    0.2433118 0.001972307
# GSM404008    0.03454889    0.002089211    0.5739203 0.000000000

re <- TME.results[,-(23:25)]

# 堆积柱状图
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
Group = str_sub(colnames(exp),1,str_length(colnames(exp)))

# 肿瘤正常--------
tcga_group <-readRDS("..\\00data\\tcga_group.rds")
group <- tcga_group
# 高低风险--------
group <- TrainGroup
re2 <- re[TrainGroup$ID,]
exp2 <- exp2[,TrainGroup$ID]
# ---------
dat <- re2 %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(group = group$group) %>%
  gather(key = Cell_type,value = Proportion,-Sample,-group) %>%
  arrange(group)

dat$Sample = factor(dat$Sample,ordered = T,levels = unique(dat$Sample)) #定横坐标顺序

write.csv(dat,file = "cell_prop_risk.csv")


# 先把group排序，然后将sample设为了因子，确定排序后的顺序为水平，所以两图的顺序是对应的。
dat2 = data.frame(a = 1:ncol(exp2),
                  b = 1,
                  group = sort(group$group))

p1 = ggplot(dat2,aes(x = a, y = b)) +
  geom_tile(aes(fill = group)) +
  scale_fill_manual(values = mypalette(22)[1:length(unique(Group))]) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")

p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) +
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

library(patchwork)

p <- p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")

ggsave(filename = '04.histon_plot_risk.pdf',p,w=15,h=8)
# 
# # #------------------------- 箱线图--------------------------
# # # 全是0的行去掉
k = colSums(re)>0;table(k)
# #
# # ## k
# # ## FALSE  TRUE
# # ##     1    21
# #
re2 = re2[,k]
library(tinyarray)

# group <- t(group)
# factor(group[2,])
# #

p <- draw_boxplot(t(re2)%>%as.data.frame(),factor(group$group,labels = c("Low","High")),
                  drop = F,
                  method = "wilcox.test",
                  color = mypalette(length(unique(Group))))+
  scale_fill_manual(values = c("#5F9EA0","#FF7F24"))+
  labs(x = "Cell Type", y = "Estimated Proportion")
p
ggsave(filename = '05.box_plot2_risk.pdf',p,w=10,h=5)

# 细胞相关性----
tiics_result <- re2
# %>% lc.tableToNum
#tiics_result <- tiics_result[, group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
#tiics_result <- tiics_result[rownames(expr),]

# diff <- read.csv('./02.cell_diff_wilcox_test.csv')
# tiics_result <- tiics_result[,diff$immune_cell]

diffcell <- c("Mast cells resting",
              "Monocytes",
              "Macrophages M0",
              "B cells naive",
              "Macrophages M2",
              "T cells CD4 memory activated",
              "Eosinophils",
              "Mast cells activated",
              "Neutrophils",
              "T cells follicular helper",
              "Dendritic cells resting"
)
res3 <- tiics_result[diffcell,]
res3 <- t(res3)%>%as.data.frame()

# #过滤掉表达为0的
# res3 <- res3[,which(colSums(res3) > 0)]
# res3 <- res3[,order(colnames(res3),decreasing = F)]
library(ggcorrplot)
cor_data <- cor(res3,method="spearman")
corp <- cor_pmat(res3)

write.csv(cor_data,'03.cell_cor_r.csv',quote=F)
write.csv(corp,'03.cell_cor_p.csv',quote=F)
env.cor <- round(cor((res3),method="spearman"), 3)
# env.p <-round(cor_pmat((gene_dat),method = "spearman"),3) 
cor_p <- WGCNA::corPvalueStudent(env.cor,nrow(res3))

pdf("cor.pdf", width = 7, height = 7)
cor.plot<-corrplot(corr =env.cor,p.mat = cor_p,type="upper",
                   col = colorRampPalette(c("blue", "white", "red"))(50),
                   tl.pos="lt",tl.col="black", 
                   insig = "label_sig", sig.level = c(.001,.01, .05),
                   pch.cex=1,pch.col = "black",order = "AOE")
cor.plot<-corrplot(corr = env.cor,type="lower",add=TRUE,method="number",
                   col = colorRampPalette(c("blue", "white", "red"))(50),
                   tl.pos="n",tl.col="black",tl.cex=1.2,
                   diag=FALSE, cl.pos="n",pch.col = "black",
                   number.cex = 0.7,order = "AOE")
dev.off()

#  生物标志物与关键免疫细胞相关性----
exprlog <- exp2
# group <- comgroup
# multi <- modelgenes

gene <- lassogene

expr <- t(exprlog) %>% as.data.frame()
expr <- expr[,gene]
colnames(group) <- c('sample', 'group')

tiics_result <- re2
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[, group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[rownames(expr),]

diff <- diffcell
tiics_result <- tiics_result[,diff]


tem <- intersect(rownames(expr), rownames(tiics_result))
expr <- expr[tem,]
tiics_result <- tiics_result[tem,]
identical(rownames(expr), rownames(tiics_result))
cor_r <- cor(expr,tiics_result,method = "spearman") 
#cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(dat_exp_diff)))
library(psych)
d <- corr.test(expr,tiics_result,use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"04.correlation_cor.csv")

# 相关性热图带显著性----
data <- read.csv('04.correlation_cor.csv',row.names = 1,check.names = F)
data <- data %>%
  mutate(text = case_when( #设置label，并加入判断，当P值符合特定条件就显示"\n"外加特定数量的*号
    Pvalue <= 0.001 ~ "***", #P<0.001就显示回车加三个星号
    between(Pvalue, 0.001, 0.01) ~ "**", #P为0.001-0.01 显示回车加两个*号
    between(Pvalue, 0.01, 0.05) ~ "*",  #P为0.01-0.05 显示回车加一个星号
    T ~ ""))
# data <- data %>%
#   mutate(text = case_when(  # 一定要 get 到 case_when() 函数奥秘
#     Pvalue > 0 ~ paste(round(Pvalue, 2), "+"), # round() 只保留两位小数
#     Pvalue < 0 ~ paste(round(Pvalue, 2), "-")))
p <- 
  ggplot(data, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") + # 这里可以用 windowns 小工具 takecolor 取色，看中哪个文章就吸哪个文章
  # 比如这篇 https://www.nature.com/articles/nmeth.1902 
  geom_text(aes(label = text),col ="black",size = 5) +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold")) + #调整y轴文字
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) +   # 修改 legend 内容
  scale_x_discrete(position = "top") #
p

ggsave(file=paste0('correlation_biomarker_','.pdf'), height = 6, width = 7, p)

