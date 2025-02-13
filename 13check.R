options(timeout = Inf)

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)

rm(list = ls())
foldpath <- "D:/workdir/12stadb/13check"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
# check----
check <- read.csv("../../basedata/key_checkpoint_gene.csv", header=T)
load("../08auc/riskgroup.RData")
TrainExp <- readRDS("../00data/tcga_fpkm.rds")

TrainExp <- TrainExp[check$Symbol,]
TrainExp <- na.omit(TrainExp)
TrainExp <- t(TrainExp) %>% data.frame()
TrainExp$sample <- rownames(TrainExp)

rownames(risk) <- risk$ID
TrainGroup <-group <- risk[,10,drop = F] %>% rownames_to_column(.,var = "sample")

aa <- merge(TrainExp,TrainGroup)
aa <- gather(data = aa,key = gene,value = exp,-c("sample","group")) 

# -------------------做一个差异基因列表
library(rstatix)
stat_res <- aa %>% 
  group_by(gene) %>% 
  wilcox_test(exp ~ group) %>% 
  adjust_pvalue(method = "BH") %>% 
  add_significance("p")
stat_res
write.csv(stat_res, file = "wilcoxon_res_diffgene.csv", row.names = F, quote = F)

deffgenes <- c()

for (i in 1:nrow(stat_res)){
  p <- stat_res[i, 8]
  if (p<0.05){
    deffgenes<-append(deffgenes, stat_res[i,1])
  }
}

stat_res_diff <- subset(stat_res,p<0.05)
write.csv(stat_res, file = "wilcoxon_res.csv", row.names = F, quote = F)
write.csv(stat_res_diff, file = "wilcoxon_res_diffgene.csv", row.names = F, quote = F)
# ------------------------------绘图

pdf("01.check.pdf",w=12,h=6)
ggplot(aa, aes(x=gene, y=exp, fill = group)) +
  scale_fill_manual(values = c("#FF7F24","#5F9EA0")) + #设置颜色
  ylim(0,12)+
  geom_boxplot()+
  stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
  # annotate(geom = "text", x = violin_dat1$gene, y = 13, size = 5,#(修改显著问题)
  #          label =as.character( violin_dat1$P.Value.signif)) +
  xlab("") + 
  ylab("Expression") + 
  labs(title = "")+ 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
dev.off()

# ----------差异免疫检查点和风险评分基因相关性分析------------
library(dplyr)
library(psych)
library(Hmisc)

TrainExp <- readRDS("../00data/tcga_fpkm.rds")
check <- read.csv("wilcoxon_res_diffgene.csv", header=T)

load("../08auc/riskgroup.RData")
riskScore <- risk[,c(1,9)]

TrainExp <- subset(TrainExp, rownames(TrainExp) %in% check$gene)
# TrainExp <- TrainExp[check$genes,]
TrainExp <- na.omit(TrainExp)

intersect_result <- intersect(riskScore$ID,colnames(TrainExp))

TrainExp <- TrainExp[,intersect_result]


# 相关性分析
TrainExp <- t(TrainExp) %>% data.frame()
identical(rownames(TrainExp), riskScore$ID)
# riskScore <- riskScore[,1,drop = F] %>% rownames_to_column(.,var = "sample")
riskScore <- data.frame(riskScore)

riskScore <- na.omit(riskScore)


library(ggplot2)
library(magrittr)
library(tibble)
library(reshape2)
library(dplyr)
library(psych)
library(Hmisc)

TrainExp <- readRDS("../00data/tcga_fpkm.rds")
multi <- read.csv("../08auc/01.lasso_genes.csv")

intersect_gene <- multi$symbol

# TrainExp <- t(TrainExp) %>% data.frame()

checkExp <- TrainExp[check$gene,]
checkExp <- na.omit(checkExp)
preExp <- TrainExp[multi$symbol,]
preExp <- na.omit(preExp)

# =========整合到一张图============
# 相关性分析
checkExp <- t(checkExp) %>% data.frame()
checkExp <- checkExp[riskScore$ID,]
preExp <- t(preExp) %>% data.frame()

preExp$ID <- rownames(preExp)
preExp <- merge(preExp,riskScore,by = "ID")
rownames(preExp) <- gsub("\\-",".",preExp$ID)
preExp <- preExp[,-1]

cor_r <- cor(preExp,checkExp,method = "spearman") 

d <- corr.test(preExp,checkExp,use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"04.pre_check_correlation_cor.csv")


# 相关性热图带显著性----
data <- read.csv('04.pre_check_correlation_cor.csv',row.names = 1,check.names = F)
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

data$gene <- factor(data$gene,levels = c("BGN","FAP","GHR","KHK","TIAM1","RS"),
                    labels = c("BGN","FAP","GHR","KHK","TIAM1","RS"))
p <- ggplot(data, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 0.5) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
  geom_text(aes(label = text), col = "black", size = 5, hjust = 0.5, vjust = 0.8) +  # 添加 hjust 和 vjust 参数
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, face = "bold"),
    legend.title = element_text(hjust = 0.5), 
    legend.text = element_text(hjust = 0.5, vjust = 0.5)
  ) +
  labs(
    fill = paste0(
      " * p < 0.05\n",
      " ** p < 0.01\n",
      " *** p < 0.001\n",
      "Correlation"
    )
  ) +
  scale_x_discrete(position = "top")
p
# ggsave(file=paste0('pre_correlation','.png'), height = 6, width = 5,p)
ggsave(file=paste0('pre_correlation','.pdf'), height = 8, width = 5, p)
