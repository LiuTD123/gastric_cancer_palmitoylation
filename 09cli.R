options(timeout = Inf)
rm(list = ls())

foldpath <- "D:/workdir/12stadb/09cli"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}

setwd(foldpath)
library(gghalves)
library(tidyverse)
library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(corrplot)
library(ggcorrplot)
library(corrplot)
library(graphics)
library(WGCNA)
library(rstatix)
library(IOBR)
library(reshape2)
library(stats)
library(dplyr)
library(psych)
library(Hmisc)
library(tidyr)
library(survivalROC)
library(survminer)
library(regplot)
library(survival)
library(rms)
library(forestplot)
library(tidyr)
library(survival)
library(patchwork)
library(magrittr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggpubr)
library(timeROC)
# 临床信息----
tcga_group <-readRDS("..\\00data\\tcga_group.rds")
clinical <- read.delim2("../00data/TCGA-STAD.clinical.tsv.gz", header = T, check.names = F)
head(clinical)
colnames(clinical)
# --------------------------------------------------------------------
# gender.demographic
# race.demographic
# age_at_index.demographic
# ajcc_pathologic_stage.diagnoses
# ajcc_pathologic_t.diagnoses
# ajcc_pathologic_m.diagnoses
# ajcc_pathologic_n.diagnoses

# table(clinical$tumor_stage.diagnoses)
# --------------------------------------------------------------------
luad_clinical <- data.frame(sample=str_sub(clinical$sample,1,16),
                            gender=clinical$gender.demographic,
                            age=clinical$age_at_index.demographic,
                            race=clinical$race.demographic,
                            mstage = clinical$ajcc_pathologic_m.diagnoses,
                            nstage = clinical$ajcc_pathologic_n.diagnoses,
                            tstage = clinical$ajcc_pathologic_t.diagnoses,
                            stage = clinical$ajcc_pathologic_stage.diagnoses
) %>% dplyr::distinct(sample, .keep_all = T)

clinical <- merge(tcga_group,luad_clinical,by="sample")
write.csv(clinical,file = "clidata.csv")
# 构建模型中的riskscor

load("../08auc/riskscore.RData")
riskScore <- rs[[1]]
colnames(riskScore)[1] <- "sample"
luad_cli_all <- merge(luad_clinical, riskScore, by = "sample")
luad_cli_all[luad_cli_all==""] <- NA

write.table(luad_cli_all, file = "luad_cli_riskscore.txt", sep = "\t", quote = F, row.names = F)
# 独立预后----

cli_dat <- read.table("luad_cli_riskscore.txt", header = T, sep = "\t")
# survival_dat <- readRDS("../00_data/survival_dat.rds") 
# colnames(survival_dat) <-  c("sample","OS.time", "OS")
# survival_dat$sample <- gsub("[A-Z]$", "", survival_dat$sample)
df_all <- cli_dat
# df_all <- merge(cli_dat, survival_dat, by = "sample")

# _____________________________________肿瘤临床信息分类和整理_______________________________________________-----
table(df_all$gender)
# female   male 
# 124    229 
range(df_all$age)
# 31 90

df_all$gender <- factor(df_all$gender, levels = c("male", "female"), label = c("Male", "Female"))

table(df_all$race)
# asian                 black or african american native hawaiian or other pacific islander 
# 77                                        11                                         1 
# not reported                                     white 
# 39                                       225
df_all$race <- factor(df_all$race, levels = c("asian",
                                              "black or african american",
                                              "native hawaiian or other pacific islander",
                                              "not reported",
                                              "white"), 
                      label = c("Asian", 
                                "Black",
                                "White",
                                "Other",
                                "White"))

table(df_all$mstage)
df_all$mstage <- factor(df_all$mstage, levels = c("M0",
                                                  "M1",
                                                  "MX"),
                        label = c("M0",
                                  "M1",
                                  "MX"))

table(df_all$nstage)
df_all$nstage <- factor(df_all$nstage, levels = c("N0",
                                                  "N1",
                                                  "N2",
                                                  "N3","N3a","N3b",
                                                  "NX"),
                        label = c("N0",
                                  "N1",
                                  "N2",
                                  "N3","N3","N3",
                                  "NX"))

table(df_all$tstage)
df_all$tstage <- factor(df_all$tstage, levels = c("T1","T1a","T1b",
                                                  "T2","T2a","T2b",
                                                  "T3",
                                                  "T4","T4a","T4b",
                                                  "TX"),
                        label = c("T1","T1","T1",
                                  "T2","T2","T2",
                                  "T3",
                                  "T4","T4","T4",
                                  "TX"))

table(df_all$stage)
df_all$stage <- factor(df_all$stage, levels = c("Stage I","Stage IA","Stage IB",
                                                  "Stage II","Stage IIA","Stage IIB",
                                                  "Stage III","Stage IIIA","Stage IIIB","Stage IIIC",
                                                  "Stage IV"),
                        label = c("Stage I","Stage I","Stage I",
                                  "Stage II","Stage II","Stage II",
                                  "Stage III","Stage III","Stage III","Stage III",
                                  "Stage IV"))

df_all$age <- as.numeric(df_all$age)

# ---------限制性立方样图确定年龄分级---------
dd <- datadist(df_all)
options(datadist='dd') 
for (knot in 3:10) {
  fit <- cph(Surv(OS.time,OS) ~ rcs(age,knot) , x=TRUE, y=TRUE,data= df_all)
  tmp <- extractAIC(fit)
  if(knot==3){AIC=tmp[2];nk1=3}
  if(tmp[2]<AIC){AIC=tmp[2];nk1=knot}
}
nk1 #3
fit <- cph(Surv(OS.time,OS) ~ rcs(age,nk1) , x=TRUE, y=TRUE,data=df_all)
cox.zph(fit,"rank")
anova(fit)
#               Wald Statistics          Response: Surv(time, etype == 1) 
#
# Factor     Chi-Square d.f. P     
# age        6.59       2    0.0371
# Nonlinear 1.24       1    0.2659
# TOTAL      6.59       2    0.0371
HR<-Predict(fit, age,fun=exp)
head(HR)

pdf("age_sep.pdf",w = 6, h =5)
ggplot()+
  geom_line(data=HR, aes(age,yhat),
            linetype="solid",size=1,alpha = 0.7,colour="#0070b9")+
  geom_ribbon(data=HR, 
              aes(age,ymin = lower, ymax = upper),
              alpha = 0.1,fill="#0070b9")+
  theme_classic()+
  geom_hline(yintercept=1, linetype=2,size=1)+
  geom_vline(xintercept=60.97918,size=1,color = '#E0367A')+#查表HR=1对应的age
  annotate('text',x = 65,y=0.8,
           label=paste0('Age = 61'),
           color='black', ) +
  # geom_vline(xintercept=65.26131,size=1,color = '#d40e8c')+
  labs(title = "Cancer Risk", x="Age", y="HR (95%CI)") 
dev.off()

df_all$age <- cut(df_all$age, breaks = c(-Inf, 61,  Inf), labels = c("61<", "≥61"))
table(df_all$age)
# 61< ≥61 
# 124 226
# ---------------------------------箱线图--------------------------------------------------------
gender <- ggboxplot(df_all,
                    x = "gender", y = "RS", 
                    color = "gender", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#A73030FF"),
                    add = "jitter", shape = "gender",
                    ggtheme = theme_bw()) +
  stat_compare_means(family = "Times", label.x = 1.3, label = "p.signif") +
  labs(title = "Gender", x = "") +
  # ylim(0,5) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
gender
ggsave(filename = "01.gender_riskScore_dis.pdf", width = 4, height = 4, gender)

# ----------------------------------------------------------------------
table(df_all$age)
# 61< ≥61 
# 124 226

age <- ggboxplot(df_all %>% dplyr::filter(!is.na(age)),
                 x = "age", y = "RS", 
                 color = "age", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#A73030FF"),
                 add = "jitter", shape = "age",
                 ggtheme = theme_bw()) +
  stat_compare_means(family = "Times", label.x = 1.3, label = "p.signif") +
  labs(title = "Age", x = "") +
  # ylim(0,5) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
age
ggsave(filename = "02.age_riskScore_dis.pdf", width = 4, height = 4, age)

# ---------------------------tstage-----------------------------------
table(df_all$tstage)
# my_comparisons <- list(c("T1",   "T2"),
#                        c("T1",   "T3"),
#                        c("T1",   "T4"),
#                        c("T1",   "TX"),
#                        c("T2","T3"),
#                        c("T2","T4"),
#                        c("T2","TX"),
#                        c("T3","T4"),
#                        c("T3","TX"),
#                        c("T4","TX"))

sig <- compare_means(RS ~ tstage, data = df_all,
              symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), 
                                 symbols = c("***", "*", "ns")))

my_comparisons <- list()
for (i in 1:nrow(sig)){
  my_comparisons[[i]] <- c(sig$group1[i],sig$group2[i]) 
}

t_stage <- ggboxplot(data = subset(df_all, !is.na(tstage)),
                     x = "tstage", y = "RS", 
                     color = "tstage", palette = c("#00AFBB", "#E7B800", "#FC4E07","#FFA07A","#E0367A"),
                     add = "jitter", shape = "tstage",
                     ggtheme = theme_bw()) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")+
  # stat_compare_means(family = "Times", label.x = 2.5) +
  labs(title = "T Stage", x = "") +
  # ylim(-1,2) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
t_stage
ggsave(filename = "04.t_stage_riskScore_dis.pdf", width = 5, height = 5, t_stage)

# ----------------------------nstage--------------------------------------
table(df_all$nstage)

sig <- compare_means(RS ~ nstage, data = df_all,
              symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), 
                                 symbols = c("***", "*", "ns")))

my_comparisons <- list()
for (i in 1:nrow(sig)){
  my_comparisons[[i]] <- c(sig$group1[i],sig$group2[i]) 
}

n_stage <- ggboxplot(data = subset(df_all, !is.na(nstage)),
                     x = "nstage", y = "RS", 
                     color = "nstage", palette = c("#00AFBB", "#E7B800", "#FC4E07","#FFA07A","#E0367A"),
                     add = "jitter", shape = "nstage",
                     ggtheme = theme_bw()) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif", hide.ns = F)+
  # stat_compare_means(family = "Times", label.x = 2.5) +
  labs(title = "N Stage", x = "") +
  # ylim(-1,2) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
n_stage
ggsave(filename = "05.n_stage_riskScore_dis.pdf", width = 7, height = 5, n_stage)

# ----------------------------mstage--------------------------------------
table(df_all$mstage)

sig <- compare_means(RS ~ mstage, data = df_all,
              symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), 
                                 symbols = c("***", "*", "ns")))

my_comparisons <- list()
for (i in 1:nrow(sig)){
  my_comparisons[[i]] <- c(sig$group1[i],sig$group2[i]) 
}

m_stage <- ggboxplot(data = subset(df_all, !is.na(mstage)),
                     x = "mstage", y = "RS", 
                     color = "mstage", palette = c("#00AFBB", "#E7B800", "#FC4E07","#FFA07A"),
                     add = "jitter", shape = "mstage",
                     ggtheme = theme_bw()) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")+
  # stat_compare_means(family = "Times", label.x = 2.5) +
  labs(title = "M Stage", x = "") +
  # ylim(-1,2) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
m_stage
ggsave(filename = "05.m_stage_riskScore_dis.pdf", width = 7, height = 5, m_stage)

# ----------------------------stage--------------------------------------
table(df_all$stage)

sig <- compare_means(RS ~ stage, data = df_all,
                     symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), 
                                        symbols = c("***", "*", "ns")))

my_comparisons <- list()
for (i in 1:nrow(sig)){
  my_comparisons[[i]] <- c(sig$group1[i],sig$group2[i]) 
}

stage <- ggboxplot(data = subset(df_all, !is.na(stage)),
                     x = "stage", y = "RS", 
                     color = "stage", palette = c("#00AFBB", "#E7B800", "#FC4E07","#FFA07A"),
                     add = "jitter", shape = "stage",
                     ggtheme = theme_bw()) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")+
  # stat_compare_means(family = "Times", label.x = 2.5) +
  labs(title = "Stage", x = "") +
  # ylim(-1,2) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
stage
ggsave(filename = "05.stage_riskScore_dis.pdf", width = 7, height = 5, stage)
# ---------------------------------race------------------------------------
# ----------------------------race dis--------------------------------------
table(df_all$race)

sig <- compare_means(RS ~ race, data = df_all,
              symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), 
                                 symbols = c("***", "*", "ns")))

my_comparisons <- list()
for (i in 1:nrow(sig)){
  my_comparisons[[i]] <- c(sig$group1[i],sig$group2[i]) 
}

race <- ggboxplot(data = subset(df_all, !is.na(race)),
                  x = "race", y = "RS", 
                  color = "race", palette = c("#00AFBB", "#E7B800", "#FC4E07","#FFA07A"),
                  add = "jitter", shape = "race",
                  ggtheme = theme_bw()) +
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")+
  # stat_compare_means(family = "Times", label.x = 2.5) +
  labs(title = "Race", x = "") +
  # ylim(-1,2) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", family = "Times"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", family = "Times"),
        axis.title = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 15, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"))
race
ggsave(filename = "05.race_riskScore_dis.pdf", width = 7, height = 5, race)
ggsave(filename = "05.race_riskScore_dis.png", width = 7, height = 5, race)
# -------------------------------------------------------------------------
all_cli_plot <- gender + age  + race + t_stage + n_stage + m_stage + stage + plot_layout(nrow = 2) 
# plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(family = "Times", face = "bold", size = 20)) 
all_cli_plot
ggsave(filename = "06.all_riskScore_dis.pdf", width = 15, height = 10, all_cli_plot)
# ------------------------------------------------------------------------------

# -----------------------------------------------临床特征整理-------------------
res.risk <- coxph(Surv(time = OS.time, event = OS) ~ RS, data = df_all) %>% summary
res.risk <- c(res.risk$conf.int[-2], res.risk$coefficients[5])

res.age <- coxph(Surv(time = OS.time, event = OS) ~ age, data = df_all) %>% summary
res.age <- c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender <- coxph(Surv(time = OS.time, event = OS) ~ gender, data = df_all) %>% summary
res.gender <- c(res.gender$conf.int[-2], res.gender$coefficients[5])

res.t_stage <- coxph(Surv(time = OS.time, event = OS) ~ tstage, data = df_all) %>% summary
res.t_stage <- cbind(res.t_stage$conf.int[,-2], res.t_stage$coefficients[,5])

res.n_stage <- coxph(Surv(time = OS.time, event = OS) ~ nstage, data = df_all) %>% summary
res.n_stage <- cbind(res.n_stage$conf.int[,-2], res.n_stage$coefficients[,5])

res.m_stage <- coxph(Surv(time = OS.time, event = OS) ~ mstage, data = df_all) %>% summary
res.m_stage <- cbind(res.m_stage$conf.int[,-2], res.m_stage$coefficients[,5])

res.stage <- coxph(Surv(time = OS.time, event = OS) ~ stage, data = df_all) %>% summary
res.stage <- cbind(res.stage$conf.int[,-2], res.stage$coefficients[,5])

res.race <- coxph(Surv(time = OS.time, event = OS) ~ race, data = df_all) %>% summary
res.race <- cbind(res.race$conf.int[,-2], res.race$coefficients[,5])

res.ref <- c(1, 1, 1, NA)
res <- rbind(res.risk, 
             res.age, 
             res.gender, 
             res.ref, res.stage,
             res.ref, res.t_stage, 
             res.ref, res.n_stage, 
             res.ref, res.m_stage,
             res.ref, res.race) %>% as.data.frame()

table(df_all$nstage)
res

res$Indicators <- c("riskScore", "Age", "Gender", 
                    "Stage\n(Stage I Reference)",
                    "II",
                    "III",
                    "IV",
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3",
                    "T4",
                    "TX",
                    "N Stage\n(N0 Reference)", 
                    "N1",
                    "N2",
                    "N3",
                    "NX",
                    "M Stage\n(M0 Reference)",
                    "M1",
                    "MX",
                    "Race\n(Asian Reference)",
                    "Black",
                    "White",
                    "Other"
)
colnames(res) <- c("hr","low","up","pv","Indicator")
res$p <- signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] <- NA
res$Indicators <- factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_prog_forest.xls", sep = "\t", quote = F)
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


# 单因素----

pdf(file = "07.univariate_cox_prog_forest.pdf", family = "Times", height = 11, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,
                          rep(FALSE,3),
                          TRUE,rep(FALSE,3), #Stage
                          TRUE,rep(FALSE,4), #T
                          TRUE,rep(FALSE,4), #N
                          TRUE,rep(FALSE,2), #M
                          TRUE,rep(FALSE,3)), #Race
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           # xticks = c(0, 1, 2, 4, 6, 16), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
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
           title = "Univariate",
           clip = c(0,7)) # 垂直于x轴的网格线，对应每个刻度

dev.off()

save.image("unicox.Rdata")

# 多因素----
features <- c("RS", "age","stage","tstage","nstage","mstage")

cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(features, collapse = "+")))
cox_more <- coxph(cox_data, data = df_all)

cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
names(cox_more$coefficients) <- c("riskScore",
                                  "Age",
                                  "Stage II",
                                  "Stage III",
                                  "Stage IV",
                                  "T2",
                                  "T3",
                                  "T4",
                                  "TX",
                                  "N1",
                                  "N2",
                                  "N3",
                                  "NX",
                                  "M1",
                                  "MX")

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

multi <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1,2,16,3:5,17,6:8,18,10:13,19,14:15),]

res2$Indicator <- c(res2$Indicator[1:2],
                    "Stage\n(Stage I Reference)",
                    "II",
                    "III",
                    "IV",
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3",
                    "T4",
                    "N Stage\n(N0 Reference)", 
                    "N1",
                    "N2",
                    "N3",
                    "NX",
                    "M Stage\n(M0 Reference)",
                    "M1",
                    "MX")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_os_CP.xls", sep = "\t", quote = F)
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

# -----------多因素画图---------------
pdf(file = "07.multi_cox.pdf", family = "Times", height = 12, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,
                          rep(FALSE,2),
                          TRUE,rep(FALSE,3),
                          TRUE,rep(FALSE,3),
                          TRUE,rep(FALSE,4),
                          TRUE,rep(FALSE,2)),
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

# ----------------------------
# PH检验
load("unicox.Rdata")

features <- c("RS", "age","stage","tstage","nstage","mstage")

x <- df_all[,c(17,16)]
uniSigExp = df_all[,c("OS.time","OS",features)]
rownames(uniSigExp) = df_all$sample

dat <- uniSigExp
outPH=data.frame()
for(i in colnames(dat[,3:8])){
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

# ----------------------------PH画图-----------------------
cox_RS <- coxph(Surv(OS.time, OS) ~ RS, data = dat)
RS <- cox.zph(cox_RS)
RS <- ggcoxzph(RS,caption = 'RS', font.main = 15)

cox_age <- coxph(Surv(OS.time, OS) ~ age, data = dat)
Age <- cox.zph(cox_age)
Age <- ggcoxzph(Age,caption = 'Age', font.main = 15)

cox_stage <- coxph(Surv(OS.time, OS) ~ stage, data = dat)
Stage <- cox.zph(cox_stage)
Stage <- ggcoxzph(Stage,caption = 'Stage', font.main = 15)

cox_tstage <- coxph(Surv(OS.time, OS) ~ tstage, data = dat)
Tstage <- cox.zph(cox_tstage)
Tstage <- ggcoxzph(Tstage,caption = 'Tstage', font.main = 15)

cox_nstage <- coxph(Surv(OS.time, OS) ~ nstage, data = dat)
Nstage <- cox.zph(cox_nstage)
Nstage <- ggcoxzph(Nstage,caption = 'Nstage', font.main = 15)

cox_mstage <- coxph(Surv(OS.time, OS) ~ mstage, data = dat)
Mstage <- cox.zph(cox_mstage)
Mstage <- ggcoxzph(Mstage,caption = 'Mstage', font.main = 15)

library(gridExtra)
ggsave("ggcoxzph_all.pdf",w = 9,h = 8,
       arrangeGrob(grobs = c(RS,Age,Stage,Tstage,Nstage,Mstage)))
# -----------------------------------------

# ------------------Nomogram-列线图1----------------------
df_all2 <- df_all

df_all2 <- df_all2[,c("OS.time","OS",features)]

ddist <- datadist(df_all2)
# df_all2$OS.time <- df_all2$OS.time*365
options(datadist = "ddist")
# cox_data2 <- as.formula(paste0('Surv(OS.time, OS)~', paste(c("riskScore","Stage","`T Stage`", "`N Stage`"), collapse = "+")))
cox_data2 <- as.formula(paste0('Surv(OS.time, OS)~', paste(features, collapse = "+")))
res.cox <- psm(cox_data2, data = df_all2, dist = "lognormal")
surv <- Survival(res.cox)
# df_all2$OS.time <- df_all2$OS.time *365
function(x) surv(365*1, x)
function(x) surv(365*3, x)
function(x) surv(365*5, x)
# 
# 报错Error in reformulate(attr(termobj, "term.labels")[-dropx], response = if (keep.response) termobj[[2L]],  : 
#                          'termlabels'必需是长度至少为一的字节矢量
nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)),
                    funlabel = c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    lp = F,
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99)
)

pdf(file = "05.nomogram_line_points.pdf", height = 9, width = 18)
par(family = "Times")
plot(nom.cox, cex.axis  = 1.2, cex.var = 1.5)
dev.off()

# --------------------regplot--列线图2----------
# 安装并加载所需的R包
# install.packages("regplot")
library("regplot")
mod2 <- coxph(formula =  as.formula(paste("Surv(OS.time,OS) ~ ", paste(features,collapse = "+"))), data = df_all2)
## Call:
## coxph(formula = as.formula(paste("Surv(time,status) ~ ", paste(colnames(pbc)[c(5, 
##     11, 14, 20)], collapse = "+"))), data = pbc)
## 
##             coef exp(coef)  se(coef)      z       p
## age    -0.069487  0.932872  0.025029 -2.776 0.00550
## bili    0.219639  1.245627  0.073620  2.983 0.00285
## copper  0.005304  1.005318  0.001977  2.683 0.00729
## stage   0.746101  2.108763  0.324986  2.296 0.02169
## 
## Likelihood ratio test=31.07  on 4 df, p=2.956e-06
## n= 186, number of events= 19 
##    (71 observations deleted due to missingness)

regplot(mod2,
        failtime = c(1*365,3*365,5*365), # 定义生存的概率标度
        plots = c("violin","boxes"), # 连续性变量形状，可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"；分类变量的形状，可选"no plot" "boxes" "bars" "spikes"
        points = T, # 截距项显示为0-100
        prfail = T)

# 展示其中一位患者的生存情况
regplot(mod2,
        observation=df_all2[1,], #用哪行观测
        obscol = "#326db1",
        failtime = c(1*365,3*365,5*365), 
        plots = c("violin","boxes"),
        droplines = T, # 是否画竖线
        points = T,
        title = "nomogram", # 更换标题
        # odds = T, # 是否显示OR值
        showP = T, # 是否显示变量的显著性标记（默认：T）
        rank = "sd", # 根据sd给变量排序
        # interval="confidence", # 展示可信区间
        clickable = F, # 是否可以交互（默认：F）
        prfail = T)

# -----------------校准曲-线1------------------
# 安装并加载所需的R包
features <- c("RS", "age","stage","tstage","nstage","mstage")
# # 建模并完成计算
set.seed(55555)
summary(df_all2)
# df_all2 <- na.omit(df_all2)

f1 <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(features,collapse = "+"))),
          data = df_all2, x = T,y = T,surv = T, time.inc = 365) # time.inc参数和calibrate的u参数后接天数
cal1 <- calibrate(f1, cmethod = "KM", method = "boot", u = 365, m = 100, B = 1000) # m，分组到平均包含 m 个受试者的区间；

f2 <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(features,collapse = "+"))),
          data = df_all2, x = T,y = T,surv = T, time.inc = 365*3) # time.inc参数和calibrate的u参数后接天数
cal2 <- calibrate(f2, cmethod="KM", method="boot", u = 365*3, m = 100, B = 1000)

f3 <- cph(formula =  as.formula(paste("Surv(OS.time, OS) ~ ",paste(features,collapse = "+"))),
          data = df_all2, x = T,y = T,surv = T, time.inc = 365*5) # time.inc参数和calibrate的u参数后接天数
cal3 <- calibrate(f3, cmethod="KM", method="boot", u = 365*5, m = 100, B = 1000)

pdf(file = "06.nomogram_predicted_1-5.pdf", family = "Times", height = 8, width = 8)
par(mar=c(5,4,2,3),cex=1.5,family="Times")
plot(cal1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',#便签
     ylab='Actual 1-5 year Overall Survival (proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal2,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',#便签
     ylab='Actual 1-5 year Overall Survival (proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Overall Survival',#便签
     ylab='Actual 1-5 year Overall Survival (proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()


# KM曲线----
features <- c("age","stage","tstage","nstage","mstage")
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

cox_data <- paste("Surv(OS.time, OS) ~ ",paste(features,collapse = "+"))
df_all2 <- df_all
res.mul = coxph(formula =  as.formula(cox_data), data = df_all2)
df_all$pred = predict(res.mul, newdata = df_all2, type = "lp")

# 中位数
# df_all$pred.group = ifelse(df_all$pred > median(df_all$pred, na.rm = T), "High", "Low")

# 最佳截断值
cox_data <- survminer::surv_cutpoint(df_all, #数据集
                                     minprop = 0.25,
                                     time = "OS.time", #生存时间
                                     event = "OS", #生存状态
                                     variables = "pred"  #需要计算的数据列名
)
cox_data
# cutpoint statistic
# pred 2.277692  5.631209
df_all$pred.group = ifelse(df_all$pred > cox_data$cutpoint$cutpoint, "High", "Low")

df_all$pred.group = factor(df_all$pred.group, levels = c("High", "Low"), labels = c("High risk", "Low risk"))
surv.fit = survfit(Surv(time = OS.time, event = OS) ~ pred.group, data = df_all)
palette_colors <- c("#f16c23","#2b6a99")

title <- "Clinical Factors"
train_km <- ggsurvplot(surv.fit,
                       pval = TRUE, 
                       pval.method = T,
                       conf.int = F,
                       legend.labs=c("High risk","Low risk" ),
                       legend.title="Risk Score",
                       xlab = "Overall Survival",
                       title=paste0(str_to_title(title)," KM"),
                       font.main = c(15,"bold"),
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       linetype = "strata", 
                       surv.median.line = "hv",
                       font.family = "Times",
                       risk.table.y.text.col = T,
                       risk.table.y.text = T,
                       risk.table.height = 0.35,
                       ggtheme = theme_bw(), 
                       palette = palette_colors)
train_km$plot <- train_km$plot + labs(
  title    = paste0("Survival Curves in ",title),
  subtitle = "Based on Kaplan-Meier Estimates"
)
train_km$table <- train_km$table + labs(
  caption  = paste0("Created with ",title," Data")
)
train_km <- customize_labels(
  train_km,
  font.title    = c(16, "bold"),
  font.subtitle = c(15, "bold.italic"),
  font.caption  = c(14, "plain", "orange"),
  font.x        = c(14, "bold.italic"),
  font.y        = c(14, "bold.italic"),
  font.xtickslab = c(12, "plain")
)
pvalue <- stringr::str_extract(train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]],
                               "\\d.*")
train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(pvalue)))

pdf_file_name <- paste0("clinical_km.pdf")
pdf(pdf_file_name, family = "Times", height = 6, width = 6, onefile = F)
print(train_km)
dev.off()

save.image("km.over.RData")
# ------------------------------ROC-------------------------------------------------
load("km.over.RData")

rt = subset(df_all2, select = c(OS, OS.time, RS))
rt$OS.time <- rt$OS.time / 365
colnames(rt) <- c("fustat","futime","riskscore")
ROC <- rt
cutoff_1 <- 1
cutoff_2 <- 3
cutoff_3 <- 5
year_1= survivalROC(Stime=ROC$futime,##生存时间
                    status=ROC$fustat,## 终止事件
                    marker = ROC$riskscore, ## marker value
                    predict.time = cutoff_1,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_2= survivalROC(Stime=ROC$futime,##生存时间
                    status=ROC$fustat,## 终止事件
                    marker = ROC$riskscore, ## marker value
                    predict.time = cutoff_2,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_3= survivalROC(Stime=ROC$futime,##生存时间
                    status=ROC$fustat,## 终止事件
                    marker = ROC$riskscore, ## marker value
                    predict.time = cutoff_3,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
if(T){
  pdf(file = paste0("01.ModelROC_train.pdf"),width = 8,height = 8)
  par(mar = c(5,5,5,2))
  plot(year_1$FP, year_1$TP,
       type="l",col="red",xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Fraction",
       ylab="True Positive Fraction",
       main="Moldel ROC",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5
  )
  abline(0,1,col="gray",lty=2)
  lines(year_2$FP, year_2$TP, type="l",col="#EB4B17",xlim=c(0,1), ylim=c(0,1))
  lines(year_1$FP, year_1$TP, type="l",col="#2775AB",xlim=c(0,1), ylim=c(0,1))
  lines(year_3$FP, year_3$TP, type="l",col="#4C8045",xlim=c(0,1), ylim=c(0,1))
  legend(0.6,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,3)),
                   paste("AUC of 3 year =",round(year_2$AUC,3)),
                   paste("AUC of 5 year =",round(year_3$AUC,3))),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,col=c("#2775AB","#EB4B17",'#4C8045'),
         bty = "n",# bty框的类型
         seg.len=1,cex=1.2)#
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
}

# ---------------特征ROC----------------------
rt_all <- df_all2
# # 创建Surv对象

rt_all <- na.omit(rt_all)

surv_obj_all <- Surv(time = rt_all$OS.time, event = rt_all$OS)
# 
# # 拟合Cox模型
cox_model_all <- coxph(surv_obj_all ~
                         RS+gender+race+
                         age+mstage+nstage+
                         tstage+stage, 
                       data = rt_all)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = rt_all, type = "risk")

ROC_all <- rt_all
# # 将风险评分添加到原始数据框中
ROC_all$riskScore <- risk_scores

# ROC_all$OS.time <- ROC_all$OS.time / 365
# --------------------------------------------------------------------------------
ROC_all = subset(ROC_all, select = c(OS, OS.time, RS, riskScore))
colnames(ROC_all) <- c("fustat","futime","RS", "mergeRS")

# -------------时间推移ROC----------
dataspan = seq(360,180*10,180)

load("../08auc/roc_merge.RData")

ROC_cligene <- timeROC(T=ROC_all$futime,   
                     delta=ROC_all$fustat,   
                     marker=ROC_all$mergeRS,   
                     cause=1,                #阳性结局指标数值
                     weighting="marginal",   #计算方法，默认为marginal
                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)

# ROC_gene <- timeROC(T=ROC_all$futime,   
#                     delta=ROC_all$fustat,   
#                     marker=ROC_all$RS,   
#                     cause=1,                #阳性结局指标数值
#                     weighting="marginal",   #计算方法，默认为marginal
#                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
#                     iid=TRUE)
# ----------绘图----------------
# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

col=c("#EB4B17", "#2775AB", '#4C8045',"#D8D155","#E0367A","#91612D","#E0867B","#35112D")
# --------时间推移----------
pdf("timeauc.pdf",w=8,h=5)
par(mai=c(1,1,1,3))
plot(dataspan, ROC_cligene$AUC, lwd=3, type = "l", col = col[1], 
     ylim = c(0.5, 0.8),
     # xlim = c(360,2000),
     xlab = "Time (Years)", ylab = "AUC", main = paste0("AUC"), bty = "l", xaxt = "n")
axis(1, dataspan)

lines(dataspan, lwd=3, ROC_merge[[1]]$AUC, col = col[2])

legendtxt <- c()
legendtxt <- c(paste0("Merged AUC=",round(mean(ROC_cligene[["AUC"]],na.rm=T),3),
                      " CI(",round(mean(confint(ROC_cligene, level = 0.95)$CI_AUC[,1])/100,2),
                      "-",round(mean(confint(ROC_cligene, level = 0.95)$CI_AUC[,2])/100,2),")\n"),
               paste0("Gene AUC=",round(mean(ROC_merge[[1]][["AUC"]],na.rm=T),3),
                                  " CI(",round(mean(confint(ROC_merge[[1]], level = 0.95)$CI_AUC[,1])/100,2),
                                  "-",round(mean(confint(ROC_merge[[1]], level = 0.95)$CI_AUC[,2])/100,2),")\n"))

# legendtxt <- c()
# for (i in c(1:6)){
#   legendtxt <- c(legendtxt,paste0(names(rs[i])," AUC=",round(median(ROC_merge[[i]][["AUC"]],na.rm=T),3)))
# }
legend(x= 1900,y= 0.7,xpd = TRUE,
       legendtxt,
       col= col,
       lty=1, lwd=3, cex = 0.8)

dev.off()
# -------------ROC曲线--------------
ROC_all$OS.time <- ROC_all$OS.time / 365
cutoff_1 <- 1
cutoff_2 <- 3
cutoff_3 <- 5

ROC <- ROC_all
ROC1 <- timeROC(T=ROC$futime,
                delta=ROC$fustat,
                marker=ROC$riskscore,
                cause=1,                #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(cutoff_1,cutoff_2,cutoff_3),       #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
ci = confint(ROC1, level = 0.95)$CI_AUC

title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)

pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC1,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC1,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC1,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC1[["AUC"]][1],3)),
         paste0("AUC at ",cutoff_2," year: ",round(ROC1[["AUC"]][2],3)),
         paste0("AUC at ",cutoff_3," year: ",round(ROC1[["AUC"]][3],3))),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
dev.off()

# -----------------单临床特征ROC-----------------------
df_all2 <- df_all[,c("OS","OS.time",features)]
df_all2$OS.time <- df_all2$OS.time/365

ROC.risk <- timeROC(T=df_all2$OS.time,
                    delta=df_all2$OS,   
                    marker=df_all2$RS,   
                    cause=1,                
                    weighting="marginal",   
                    times=c(cutoff_1,cutoff_2,cutoff_3),   
                    iid=TRUE)

df_all2$age <- as.numeric(df_all2$age)
ROC.age <- timeROC(T=df_all2$OS.time,
                   delta=df_all2$OS,   
                   marker=df_all2$age,   
                   cause=1,                
                   weighting="marginal",   
                   times=c(cutoff_1,cutoff_2,cutoff_3),   
                   iid=TRUE)

df_all2$nstage <- as.numeric(df_all2$nstage)
ROC.nstage <- timeROC(T=df_all2$OS.time,
                      delta=df_all2$OS,   
                      marker=df_all2$nstage,   
                      cause=1,                
                      weighting="marginal",   
                      times=c(cutoff_1,cutoff_2,cutoff_3),   
                      iid=TRUE)

df_all2$tstage <- as.numeric(df_all2$tstage)
ROC.tstage <- timeROC(T=df_all2$OS.time,
                               delta=df_all2$OS,   
                               marker=df_all2$tstage,   
                               cause=1,                
                               weighting="marginal",   
                               times=c(cutoff_1,cutoff_2,cutoff_3),   
                               iid=TRUE)

df_all2$mstage <- as.numeric(df_all2$mstage)
ROC.mstage <- timeROC(T=df_all2$OS.time,
                               delta=df_all2$OS,   
                               marker=df_all2$mstage,   
                               cause=1,                
                               weighting="marginal",   
                               times=c(cutoff_1,cutoff_2,cutoff_3),   
                               iid=TRUE)

df_all2$stage <- as.numeric(df_all2$stage)
ROC.stage <- timeROC(T=df_all2$OS.time,
                      delta=df_all2$OS,   
                      marker=df_all2$stage,   
                      cause=1,                
                      weighting="marginal",   
                      times=c(cutoff_1,cutoff_2,cutoff_3),   
                      iid=TRUE)
n = 3
timecut = as.numeric(c(cutoff_1,cutoff_2,cutoff_3)[n])
title = paste0("Factor ROC in ",timecut," years")

pdf(file = paste0(timecut,"_factor_train_specific.pdf"),width = 5,height = 5)
plot(ROC.risk, time = timecut, col="#ff9898", lwd=2, title = "")
plot(ROC.age, time = timecut, col="#A65628", lwd=2, add = T)
plot(ROC.nstage, time = timecut, col="#4DAF4A", lwd=2, add = T)
plot(ROC.tstage, time = timecut, col="#377EB8", lwd=2, add = T)
plot(ROC.mstage, time = timecut, col="#dedb8e", lwd=2, add = T)
plot(ROC.stage, time = timecut, col="#007172", lwd=2, add = T)
title(title)
legend("bottomright",
       c(paste0("Risk score: ",round(ROC.risk[["AUC"]][n],3)), 
         paste0("age: ",round(ROC.age[["AUC"]][n],3)),
         paste0("Stage: ",round(ROC.stage[["AUC"]][n],3)),
         paste0("N: ",round(ROC.nstage[["AUC"]][n],3)),
         paste0("T: ",round(ROC.tstage[["AUC"]][1],3)),
         paste0("M: ",round(ROC.mstage[["AUC"]][n],3))
       ),
       col=c("#ff9898", "#A65628", "#4DAF4A","#377EB8","#dedb8e","#007172"),
       lty=1, lwd=2,bty = "n")  
dev.off()
# ------------------TCGA样本中不同预后基因之间的相关性--------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
#           "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
#           "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr', 'ComplexHeatmap')
# for(i in 1:length(depens)){
#   depen<-depens[i]
#   if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
# }
# if (!requireNamespace("IOBR", quietly = TRUE))
#   devtools::install_github("IOBR/IOBR")

gene_expr <- readRDS("../00data/tcga_fpkm.rds")
multi <- read.table("../07CoxLasso/Lasso_Coefficients.xls", header = T) %>% as.data.frame()#cox_result_step2.rds;mul_cox_result.rds
pregene <- read.csv("../07CoxLasso/hub_gene.csv")

multi <- multi[multi$gene %in% pregene$x,]

gene_expr <- t(gene_expr) %>% as.data.frame()
# multi <- t(multi) %>% as.data.frame()
gene_expr <- na.omit(gene_expr)

pre_gene <- multi[,1]
gene_expr <- gene_expr[,pre_gene]
gene_expr <- t(gene_expr) %>% as.data.frame()
cor_r <- cor(gene_expr,multi[,2],method = "spearman") 

d <- corr.test(gene_expr,multi[,2],use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)#转换数据长短 cor
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)#转换数据长短 p

cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
cor_dat$cell <- "different_gene_cor"

write.csv(cor_dat,"04.dif_pregene_cor.csv")

# 预后基因相关性

cor_data <- read.csv("04.dif_pregene_cor.csv")

gene_expr <- t(gene_expr) %>% as.data.frame()
cor_data <- cor(gene_expr,method="spearman")
corp <- cor_pmat(gene_expr)
multi <- t(multi) %>% as.data.frame

write.csv(cor_data,'03.cell_cor_r.csv',quote=F)
write.csv(corp,'03.cell_cor_p.csv',quote=F)
env.cor <- round(cor((gene_expr),method="spearman"), 3)
# env.p <-round(cor_pmat((gene_dat),method = "spearman"),3)
cor_p <- WGCNA::corPvalueStudent(env.cor, nrow(gene_expr))

pdf("dif_pregene_cor_all.pdf", width = 5, height = 5)
corrplot(corr =env.cor,p.mat = cor_p,type="upper",
         col = colorRampPalette(c("blue", "white", "red"))(50),
         tl.pos="lt",tl.col="black", 
         insig = "label_sig", sig.level = c(.001,.01, .05),
         pch.cex=1,pch.col = "black",order = "AOE")
corrplot(corr = env.cor,type="lower",add=TRUE,method="number",
         col = colorRampPalette(c("blue", "white", "red"))(50),
         tl.pos="n",tl.col="black",tl.cex=1.2,
         diag=FALSE, cl.pos="n",pch.col = "black",
         number.cex = 0.7,order = "AOE")

dev.off()

# ------------------------“GOSemSim”对预后基因之间的功能相似性评分进行计算------------
# 数据准备------

# 数据准备

hubgene_id <- mapIds(x = org.Hs.eg.db,#注释包
                     keys = colnames(gene_expr), #需要转换的基因Symbol
                     keytype = "SYMBOL", #需要转换的类型
                     column = "ENTREZID") #需要转换为的类型

print(hubgene_id)

hubgene_id <- data.frame(hubgene_id)
colnames(hubgene_id) <- 'ENTREZID'
hubgene_id$SYMBOL <- rownames(hubgene_id)
hubgene_id$ENTREZID <- as.character(hubgene_id$ENTREZID)
rt <- hubgene_id

# 计算相似性
bp <- godata('org.Hs.eg.db', ont= "BP", computeIC = FALSE)
cc <- godata('org.Hs.eg.db', ont= "CC", computeIC = FALSE)
mf <- godata('org.Hs.eg.db', ont= "MF", computeIC = FALSE)

#完成注释文件准备，计算相似度
simbp<- mgeneSim(rt$ENTREZID, semData= bp, measure= "Wang", drop= NULL, combine= "BMA")
simcc<- mgeneSim(rt$ENTREZID, semData= cc, measure= "Wang", drop= NULL, combine= "BMA")
simmf<- mgeneSim(rt$ENTREZID, semData= mf, measure= "Wang", drop= NULL, combine= "BMA")

#得到评分
fsim<- (simmf * simcc * simbp)^( 1/ 3)
colnames(fsim) <- rt$SYMBOL
rownames(fsim) <- rt$SYMBOL

#去除基因与本身基因之间的相关性，使用melt函数将宽格式数据转化成长格式数据
for(i in 1:ncol(fsim))
{ fsim[i,i] <- NA }
dat <- melt(fsim)
dat <- dat[! is.na(dat$ value),]
dat <- dat[,c(1,3)]
head(dat)

#绘图
dat.mean <- aggregate(value~Var1, dat, mean)
m<- dat.mean$value 
names(m) <- dat.mean$Var1
dat$Var1 <- factor(dat$Var1, levels=names(sort(m))) 
str(dat)
#云雨图
colnames(dat)<-c("genes","value")

mycolor<-c("#85C0D9","#FFEFDB","#CE97B0","#F4A9A8","#9DCD82","#F8B62D")
p<-ggplot(dat,aes(x=genes,y=value,fill=factor(genes)))+    #建立映射
  scale_color_manual(values=rev(mycolor))+
  scale_fill_manual(values=rev(mycolor))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),    #绘制一半的小提琴图
                   side='R',adjust=1.2,trim=F,alpha=0.8)+   #参数调整：
  #position：位置调整，这里将其向右水平移动0.1；
  #side：显示哪一侧， "I"代表左侧，"R"代表右侧，默认"I"；
  #adjust：调整带宽，这里设为1.2使宽带略变平滑；
  #trim：小提琴图尾部的数据修整，默认为"T",表示将尾部修整到数据范围；"F"表示不修剪尾部；
  geom_point(aes(x = as.numeric(genes)-0.1,#散点位置向左平移0.1
                 y = value,
                 color = factor(genes)),
             position = position_jitter(width =0.03),size =0.8, shape = 20)+   #调整散点，使取值相同的原重合散点分散开
  geom_boxplot(outlier.shape = NA, #隐藏离群点；
               width =0.1,         #在散点和二分之一小提琴图中间添加箱线图
               alpha=0.7)+
  coord_flip()+                    #图形翻转
  theme_bw()+                      #去掉灰底
  theme(panel.grid=element_blank())#去掉背景网格线
p
ggsave('01.Friends.png',p,width =8,height = 6)
ggsave('01.Friends.pdf',p,width =8,height = 6)

# 
# #箱线图
# png('02.friends_boxplot.png',w=500,h=400)
# ggplot(dat,aes(x=genes,y=value,fill=factor(genes)))+
#   geom_boxplot(width = 0.4)+
#   coord_flip()+
#   xlab("")+ylab("")+
#   theme_bw()+
#   theme(legend.position = 'none')+
#   theme(axis.title.x =element_text(size=20, face = "bold"),
#         axis.text.x =element_text(size=16, face = "bold",color='black'),
#         axis.title.y =element_text(size=20, face = "bold"),
#         axis.text.y=element_text(size=16, face = "bold",color='black'),
#         strip.text = element_text(size = 14, face = "bold"),
#         legend.position = "none")
# dev.off()
# 
# pdf('02.friends_boxplot.pdf',w=6,h=5)
# ggplot(dat,aes(x=genes,y=value,fill=factor(genes)))+
#   geom_boxplot(width = 0.4)+
#   coord_flip()+
#   xlab("")+ylab("")+
#   theme_bw()+
#   theme(legend.position = 'none')+
#   theme(axis.title.x =element_text(size=20, face = "bold"),
#         axis.text.x =element_text(size=16, face = "bold",color='black'),
#         axis.title.y =element_text(size=20, face = "bold"),
#         axis.text.y=element_text(size=16, face = "bold",color='black'),
#         strip.text = element_text(size = 14, face = "bold"),
#         legend.position = "none")
# dev.off()
# 
