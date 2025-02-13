rm(list = ls())
setwd("D:\\workdir\\10stad/00data/")

# TCGA数据----
# ------------------------------------------TCGA数据----------------------------
library(tinyarray)
library(dplyr) 
library(tidyr)
library(GEOquery)
library(limma)
library(affy)
library(AnnoProbe)
library(Biobase)
library(tidyverse)
library(EnsDb.Hsapiens.v86)

genesymbol <- read.table("gencode.v36.annotation.gtf.gene.probemap",sep = "\t",check.names = F,header = T,row.names = 1)
genesymbol$id <- rownames(genesymbol)
which(duplicated(genesymbol$gene))
# cuont数据
train_count <- read.table("TCGA-STAD.star_counts.tsv.gz",sep = "\t",check.names = F,header = T,row.names = 1)
count1 <- trans_array(train_count,genesymbol,from = "id",to="gene")

# BiocManager::install('EnsDb.Hsapiens.v86')

edb <- EnsDb.Hsapiens.v86
columns(edb)
keytypes(edb)
keys <- keys(edb, keytype="GENEID")
## Get the data
gene2sym<-select(edb, keys=keys, 
                 columns=c("SYMBOL","ENTREZID","GENEBIOTYPE",'GENENAME'),
                 keytype="GENEID")
mrnagene <- gene2sym[gene2sym$GENEBIOTYPE == "protein_coding",]
count1 <- count1[mrnagene$SYMBOL,]
count1 <- na.omit(count1)
#  lncRNA数据
# count1 <- trans_exp(train_count,lncrna_only = T)    

count2 <- 2^count1-1
count2 <- round(count2)   
samples_group <- colnames(train_count) %>% as.data.frame()
colnames(samples_group) <- "sample"
samples_group$group <-ifelse(as.numeric(substr(samples_group$sample,14,15)) < 10,'Tumor','Normal')
samples_group <- samples_group[order(samples_group$group,decreasing = T),]

table(samples_group$group)
# Normal  Tumor 
# 36    412
tcga_count <- count2[,match(samples_group$sample,colnames(count2))]
# x <- intersect(tem,colnames(tcga_count))
# tcga_count <- tcga_count[,match(x,colnames(tcga_count))]
# colnames(tcga_fpkm) <- substr(colnames(tcga_fpkm),1,15)

##fpkm获取
train_fpkm <- read.table("TCGA-STAD.star_fpkm.tsv.gz",sep = "\t",check.names = F,header = T,row.names = 1)

tcga_fpkm <- trans_array(train_fpkm,genesymbol,from = "id",to="gene")
tcga_fpkm <- tcga_fpkm[mrnagene$SYMBOL,]

# tcga_fpkm <- tinyarray::trans_exp(train_fpkm,mrna_only = T)#筛选出蛋白基因
# lncRNA数据
# tcga_fpkm <- tinyarray::trans_exp(train_fpkm,lncrna_only = T)
# tcga_fpkm <- tcga_fpkm[,match(x,colnames(tcga_fpkm))]
# tcga_fpkm <- tcga_fpkm[,match(samples_group$sample,colnames(tcga_fpkm))]

colnames(tcga_fpkm) <- substr(colnames(tcga_fpkm),1,16)
samples_group$sample <- substr(samples_group$sample,1,16)

# tem <- unique(colnames(tcga_count))
# tem <- as.data.frame(tem)
# tem <- substr(tem$tem,1,15)
# tem <- as.character(tem)
# samples_group <- samples_group[match(x,samples_group$sample),]

colnames(tcga_count) <- substr(colnames(tcga_count),1,16)
# tcga_count <- tcga_count[,match(tem,colnames(tcga_count))]
table(samples_group$group)

# 保存
# write.table(samples_group1,file = "00.tcga_group_RR.txt",sep = "\t",quote=F,row.names = F)
write.table(samples_group,file = "00.tcga_group.txt",sep = "\t",quote=F,row.names = F)   
write.table(tcga_count,file = "00.tcga_count.txt",sep = "\t",quote=F,row.names = T)  
write.table(tcga_fpkm,file = "00.tcga_fpkm.txt",sep = "\t",quote=F,row.names = F)   

# saveRDS(samples_group1,"../00_data/tcga_group_RR.rds")
saveRDS(samples_group,"tcga_group.rds")
saveRDS(tcga_count,"tcga_count.rds")

saveRDS(tcga_fpkm,"tcga_fpkm.rds")

# survival 
samples_group <- read.table("00.tcga_group.txt",sep = "\t",header = T,check.names = F)
survival <- read.table("TCGA-STAD.survival.tsv.gz",sep = "\t",header = T,check.names = F)
survival$sample <- substr(survival$sample ,1,16)
# survival <- survival[survival$OS.time > 30,]
usesample <- samples_group$sample[samples_group$group == "Tumor"]
use_sample <- intersect(usesample,survival$sample) 
saveRDS(use_sample,"used_sample.rds")

survival <- survival[match(use_sample,survival$sample),] %>% na.omit()
survival_dat <- survival[,c(1,2,3)]
colnames(survival_dat) <- c("sample","futime","fustat")
survival_dat$futime <- survival_dat$futime/365 %>% as.numeric()
# survival_dat$OS.time <- survival_dat$OS.time/365 %>% as.numeric()
# survival_dat$DSS.time <- survival_dat$DSS.time/365 %>% as.numeric()
saveRDS(survival_dat,"survival_dat.rds")








