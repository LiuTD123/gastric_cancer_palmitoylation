rm(list = ls())

foldpath <- paste0("D:/workdir/12stadb/04venn")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
#韦恩图（VennDiagram 包，适用样本数 2-5）
library(VennDiagram)
library(ggvenn)
library(readr)
library(readxl)
library(tidyverse)

csrgs <- read.csv("../03relatgene/hub_gene.csv")
wgcna <- read.csv("../02wgcna/04.WGCNA_genes.csv")

inter <- Reduce(intersect,list(wgcna$symbol,csrgs$x))
# inter <- Reduce(intersect,list(rownames(DEGs1),DEGs2$X,cluster$genesymbol))

write.csv(inter,"wgcna_relategenes.csv",row.names = T)
vennlist <- list("WGCNA"=wgcna$symbol,"Palmitoylation"=csrgs$x)
pdf('04.venn.pdf',w=6,h=8)
ggvenn(vennlist,
       c("WGCNA","Palmitoylation"),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = T,
       fill_alpha=0.5,
       stroke_alpha=1,
       stroke_size=0.4,
       text_size=5,
       stroke_color=NA,
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       set_name_size=6,
       text_color="black")
dev.off()
