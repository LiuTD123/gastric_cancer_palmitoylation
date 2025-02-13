rm(list = ls())
foldpath <- "D:/workdir/12stadb/17expr"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(pROC)
library(glmnet)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(tidyr)
# -------------------------------------------训练集表达量验证
check <- read.csv("..\\08auc/01.lasso_genes.csv")
TrainExp <- readRDS("..\\00data\\tcga_fpkm.rds")
group <- readRDS("..\\00data\\tcga_group.rds")

TrainExp <- TrainExp[check$symbol,]
TrainExp <- na.omit(TrainExp)
TrainExp <- t(TrainExp) %>% data.frame()
TrainExp$sample <- rownames(TrainExp)

risk <- group
rownames(risk) <- risk$sample

TrainGroup <- group <- risk[,2,drop = F] %>% rownames_to_column(.,var = "sample")

aa <- merge(TrainExp,TrainGroup)
aa <- gather(data = aa,key = gene,value = exp,-c("sample","group")) 

# -------------------做一个差异基因列表
# library(rstatix)
# stat_res <- aa %>% 
#   group_by(gene) %>% 
#   wilcox_test(exp ~ group) %>% 
#   adjust_pvalue(method = "BH") %>% 
#   add_significance("p")
# stat_res
# # write.csv(stat_res, file = "wilcoxon_res_diffgene.csv", row.names = F, quote = F)
# 
# deffgenes <- c()
# 
# for (i in 1:nrow(stat_res)){
#   p <- stat_res[i, 8]
#   if (p<0.05){
#     deffgenes<-append(deffgenes, stat_res[i,1])
#   }
# }
# 
# stat_res_diff <- subset(stat_res,p<0.05)
# write.csv(stat_res_diff, file = "wilcoxon_res_diffgene.csv", row.names = F, quote = F)
# ------------------------------绘图

violin_dat1 <- read.table("..\\01DEG\\01.DEG_ALL.txt",sep = "\t",header = T)
violin_dat1$gene <- rownames(violin_dat1)
violin_dat1 <- violin_dat1[check$symbol,]

data <- violin_dat1 %>%
  mutate(text = case_when( #设置label，并加入判断，当P值符合特定条件就显示"\n"外加特定数量的*号
    pvalue <= 0.001 ~ "***", #P<0.001就显示回车加三个星号
    between(pvalue, 0.001, 0.01) ~ "**", #P为0.001-0.01 显示回车加两个*号
    between(pvalue, 0.01, 0.05) ~ "*",  #P为0.01-0.05 显示回车加一个星号
    T ~ ""))

pdf("01.train_gene.pdf",w=6,h=4)
ggplot(aa, aes(x=gene, y=exp, fill = group)) +
  scale_fill_manual(values = c("#5F9EA0","#FF7F24")) + #设置颜色
  # ylim(3,15)+
  geom_boxplot()+
  # stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
  annotate(geom = "text", x = data$gene, y = 13, size = 5,#(修改显著问题)
           label =as.character(data$text)) +
  xlab("") + 
  ylab("Expression") + 
  labs(title = "")+ 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
  theme(axis.text.x=element_text(face = "bold", color="black",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
dev.off()

# ------------------------------绘图

pdf("01.veri_gene.pdf",w=16,h=6)
ggplot(veriaa, aes(x=gene, y=exp, fill = group)) +
  scale_fill_manual(values = c("#5F9EA0","#FF7F24")) + #设置颜色
  # ylim(3,15)+
  geom_boxplot()+
  stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
  annotate(geom = "text", x = violin_dat1$gene, y = 13, size = 5,#(修改显著问题)
           label =as.character( violin_dat1$P.Value.signif)) +
  xlab("") + 
  ylab("Expression") + 
  labs(title = "")+ 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
dev.off()

png("01.veri_gene.png",w=1600,h=600)
ggplot(veriaa, aes(x=gene, y=exp, fill = group)) +
  scale_fill_manual(values = c("#5F9EA0","#FF7F24")) + #设置颜色
  ylim(3,15)+
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