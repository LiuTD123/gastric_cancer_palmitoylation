options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/12stadb/15medicinsens"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

# BiocManager::install(c('sva', 'car', 'genefilter', 'preprocessCore', 'ridge'))
# install.packages("../../pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(devtools)
library(Ipaper)
library(rstatix)
library(magrittr)
library(parallel)
library(car)
library(pRRophetic)

load("../08auc/riskgroup.RData")
# riskscore$risk <- ifelse(riskscore$riskScore>median(riskscore$riskScore), "High", "Low")

riskscore <- risk
rownames(riskscore) <-riskscore$ID

# model_coef <- read.table("/data/nas1/liuhouyan/01-BJTC-577-2/05_model_build/Lasso_Coefficients.xls", header = T)
# model_gene <- model_coef$gene
# length(model_gene)# 4

train_fpkm <- readRDS("../00data/tcga_fpkm.rds")
exprData <- train_fpkm[, riskscore$ID]
exprData <- as.matrix(log2(exprData + 1))

ic50 <- data.frame(riskscore$sample)
a <- data.frame(row.names=riskscore$sample,riskscore$group)
colnames(a)<-'group'
drug<-read.table('../../basedata/Medicinal_Sensity_drugs.txt',sep='\t',header=F)
# exprData <- as.data.frame(exprData)
# exprData <- na.omit(exprData)

# -------------修改----------------
cnt<-1
# 在calcPhenotype函数的第6，8，10，12和66行class函数后面加个[1]
trace(calcPhenotype, edit = T)
# summarizeGenesByMean函数的第19行，trace(summarizeGenesByMean, edit = T)
trace(summarizeGenesByMean, edit = T)
while (cnt < 139) {
  print(cnt)
  
  predictedPtype <- pRRopheticPredict(as.matrix(exprData), drug[cnt,],selection=1)
  
  Tipifarnib<-data.frame(predictedPtype)
  
  colnames(Tipifarnib)<-drug[cnt,]
  
  a<-cbind(a,Tipifarnib)
  
  cnt = cnt + 1
}

write.csv(a,'01.IC50.csv')
b<-a

# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}

c<-removeColsAllNa(b)

na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
# View(x)
dim(x)# 353 139

risk <- riskscore
high_group <- rownames(risk)[which(risk$group=="High")]
low_group <- rownames(risk)[which(risk$group=="Low")]

# x$sample <- rownames(x)
# dat <- t(x)%>%as.data.frame()
# dat <- dat[-c(1,140),]

dat <- x %>% as.data.frame() %>% tibble::rownames_to_column(var = "ID") %>% 
  tidyr::gather(., drug, value, -c("group", "ID"))

rTable <- dat %>% 
  group_by(drug) %>% 
  wilcox_test(value ~ group) %>% 
  adjust_pvalue(method = "BH") %>% 
  add_significance("p.adj")

rTable$sig <- ifelse(rTable$p.adj < 0.05,
                     ifelse(rTable$p.adj < 0.01, 
                            ifelse(rTable$p.adj < 0.001,
                                   ifelse(rTable$p.adj < 0.0001,
                                          paste(rTable$drug, "****",  sep = ""),
                                          paste(rTable$drug, "***", sep = "")),
                                   paste(rTable$drug, "**", sep = "")),
                            paste(rTable$drug, "*",  sep = "")), 
                     rTable$drug)

write.csv(rTable,file = "02.drugs_wilcox_res.csv")

rTable2 <- rTable[which(rTable$p.adj < 0.05),]
rTable2 <- rTable2[order(rTable2$p.adj),]
dim(rTable2)
# 104  11

dat2 <- subset(x, select = -group)
dat2 <- dat2[which(colnames(dat2) %in% rTable2$drug)]
dat3 <- exprData %>% t %>% as.data.frame()
pre_gene <- read.csv("../08auc/01.lasso_genes.csv")
dat3 <- dat3[,pre_gene$symbol]

cor_res <- cor(dat2, dat3, method = "spearman")
cor_p <- WGCNA::corPvalueStudent(cor_res, nrow(dat2))
cor_padj <- p.adjust(cor_p, method = "BH")
cor_padj <- matrix(cor_padj, nrow = dim(cor_p)[1])
rownames(cor_padj) <- rownames(cor_p)
colnames(cor_padj) <- colnames(cor_p)

cor_res2 <- cor_res %>% as.data.frame %>% tibble::rownames_to_column(var = "drug") %>% 
  tidyr::gather(., gene, Cor, -drug)
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "drug") %>% 
  tidyr::gather(., gene, Pvalue, -drug)
cor_p3 <- cor_padj %>% as.data.frame %>% tibble::rownames_to_column(var = "drug") %>% 
  tidyr::gather(., gene, Padj, -drug)
cor_res3 <- cbind(cor_res2, cor_p2, cor_p3)[,c(1,2,3,6,9)]
write.table(cor_res3, file = "03.drug_modelgene_cor.xls", sep = "\t", quote = F, row.names = F)


cor_plot <- cor_res3[which(cor_res3$drug %in% rTable2$drug[1:10]),]

ggheatmap <- ggplot(cor_plot, aes(gene, drug, fill = Cor))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 30, hjust = 1),
        axis.text.y = element_text(size = 30)) +
  coord_flip()
ggheatmap +
  geom_text(aes(gene, drug, label = paste0(signif(Cor,2),"\n(",signif(Padj,2),")")), color = "black", size =10, family = "Times") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    text = element_text(face = "bold", family = "Times"))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggsave(filename = "04.cor_heatmap.png", height = 10, width = 20, bg = "white")
dev.off()
ggsave(filename = "04.cor_heatmap.pdf", height = 10, width = 20, bg = "white")


# -------------------------------箱线图------------------------------------------

calc_res <- read.csv('01.IC50.csv',row.names = 1)
calc_res <- calc_res[,-1]

risk_score <- riskscore

calc_res <- calc_res[rownames(risk_score),]


identical(rownames(risk_score),rownames(calc_res))

# risk_score <- risk_score[,4]
# risk_score <- as.data.frame(risk_score)

drug <- colnames(calc_res)
high_group <- rownames(risk_score)[which(risk_score$group=="High")]
low_group <- rownames(risk_score)[which(risk_score$group=="Low")]
calc_res$risk <- risk_score$group

rTable <- read.csv("02.drugs_wilcox_res.csv")

rTable2 <- rTable[which(rTable$p.adj < 0.05),]
rTable2 <- rTable2[order(rTable2$p.adj),]
write.csv(rTable2,'02.drugs_wilcox_diff_res.csv')

rTable3 <-rTable2[c(1:10),]
# rTable3 <-rTable2

# 箱线图----

data<-calc_res[,rTable3$drug]
data<-t(data)
cell.data <- data.frame(Immune_Cell=rownames(data), 
                        data, 
                        pvalue=rTable3$p.adj)
plot.cell <- cell.data[which(cell.data$pvalue<0.05),]
#plot.cell<-cell.data
diff_tiics <- rownames(plot.cell)

violin_dat <- gather(plot.cell, 
                     key=Group,
                     value=score, 
                     -c("Immune_Cell","pvalue")
)

risk_score$sample <- rownames(risk_score)
train_group <-risk_score[,c(11,10)]

# train_group$sample <- gsub("-",".",train_group$sample)
group_case<-train_group[train_group$group=='High',]
group_case<-as.character(group_case$sample)
group_control<-train_group[train_group$group=='Low',]
group_control<-as.character(group_control$sample)


# # 替换字符
# for (i in group_case){
#   b <- sub("\\-", ".", i)
#   group_case1 <- c(group_case1,b)
# }
violin_dat$Group<-gsub("\\.","-",violin_dat$Group)

violin_dat$Group <- ifelse(violin_dat$Group %in% group_control,
                           "Low", "High") 
# violin_dat$Group <- ifelse(violin_dat$Group == group_case,
#                            "Low", "High")
violin_dat$Group <- factor(violin_dat$Group, levels = c("Low", "High"))
# violin_dat <- violin_dat[,-2]
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) + ylim(-10,12)+
  
  # geom_violin(trim=T,color=alpha('black',alpha = 0.5)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  # stat_boxplot(geom="errorbar",
  #              width=0.1,
  #              position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = 21,
               outlier.fill = "black",
               outlier.size = 0.5,
  )+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  # geom_point(aes(fill = Group),
  #            size = 0.05,
  #            position = position_dodge(0.9))+
  scale_fill_manual(values= c("#45a9b8","#f76a56"))+ #设置填充的颜色
  labs(title="", x="", y = "IC50",size=20) +
  # stat_compare_means(data = violin_dat,
  #                    method = "wilcox.test", #没有直接用差异分析的结果，但检验方法是一样的
  #                    mapping = aes(group = Group),
  #                    label ="p.signif",
  #                    hide.ns = F) +
  annotate(geom = "text", x = rTable3$drug, y = 11, size = 4,#(修改显著问题)
           label =as.character( rTable3$p.adj.signif)) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        text=element_text(family = 'Times'),
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs
# ggsave('03.Box.pdf',boxplot_diff_TIICs,w=10,h=6)
ggsave('03.Box.png',boxplot_diff_TIICs,w=10,h=6)
