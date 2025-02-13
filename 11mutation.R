rm(list = ls())

foldpath <- "D:/workdir/12stadb/11mutation"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

library(TCGAmutations)
library(TCGAbiolinks)
library(maftools)
# BiocManager::install("PoisonAlien/TCGAmutations")

# --------------------------选择癌症类型-----------------
maf <- TCGAmutations::tcga_load(study = "STAD")

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

riskscore <- TrainGroup
# table(riskscore$risk)
High.sample <- riskscore$ID[which(riskscore$group == "High")]
Low.sample <- riskscore$ID[which(riskscore$group == "Low")]
maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
colnames(riskscore)[1] <- 'sample'
#maf_sample <- merge(maf_sample, riskscore, by= "sample")
maf_sample <- merge(maf_sample, riskscore, by = 'sample')
sample <- subset(maf_sample)$barcode

genes <- read.csv("../08auc/01.lasso_genes.csv")
genes <- genes$symbol
maf <- subsetMaf(maf, tsb = sample,
                 genes = genes) #只选择模型基因

maf@gene.summary$Hugo_Symbol 
# 可视化-----
##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample, ]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
pdf(file = "01.oncoplot.high_model.pdf", height = 6, width = 8)
oncoplot(maf = maf.high, top = 20,
         # fontSize = 0.9,
         SampleNamefontSize = 1.3,
         legendFontSize = 2,
         barcode_mar = 2,
         titleText = "High Risk"
)

dev.off()

low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
pdf(file = "02.oncoplot.low_model.pdf", height = 6, width = 8)
oncoplot(maf = maf.low, top = 20,
         # fontSize = 0.9,
         SampleNamefontSize = 1.3,
         legendFontSize = 2,
         barcode_mar = 2,
         titleText = "Low Risk"
)
dev.off()

# 
# # MAF整体结果图-------------
pdf(file = "03.MAF_total_res_model.pdf", height = 8, width = 10)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
class(maf_sample)
# 
# # TMB-----
tmb_table_wt_log = tmb(maf = maf)
head(tmb_table_wt_log)

TMB <- tmb_table_wt_log[,c(1,4)]
colnames(TMB) <- c('sample',"TMB")
TMB <- data.frame(TMB)
TMB$sample <- as.character(TMB$sample)
TMB$sample <- stringr::str_sub(TMB$sample, 1, 15)
TMB <- merge(TMB, riskscore, by = 'sample')

pdf("04.TMB_model.pdf",w=8,h=6)
ggplot(TMB, aes(x=risk, y=TMB, fill = risk)) +
  scale_fill_manual(values = c("red","#20B2AA")) + #设置颜色
  ylim(-0.7,2.3)+
  geom_boxplot()+
  stat_compare_means(aes(group = risk),method = 'wilcox.test',label = "p.signif",cex=4)+
  # annotate(geom = "text", x = violin_dat1$gene, y = 4250, size = 5,#(修改显著问题)
  #          label =as.character( violin_dat1$P.Value.signif)) +
  xlab("") +
  ylab("TMB value") +
  labs(title = "")+
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle = 70表示横轴标题倾斜70度
dev.off()