options(timeout = Inf)

rm(list = ls())
foldpath <- "D:/workdir/12stadb/14gsea"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
#  09_GSEA----
library(enrichplot)
library(psych)
# BiocManager::install("GSEABase",force = TRUE)
library(GSEABase)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
require(DOSE)
library(Hmisc)
library(future)
library(future.apply)
library(limma)
library(psych)
library(RColorBrewer)
# BiocManager::install("GSVA")
library(GSVA)
# devtools::install_github("junjunlab/GseaVis",force = TRUE)
library(GseaVis)
# devtools::install_github("satijalab/seurat-data")

KEGG_df = msigdbr(species = "Homo sapiens",category = "C2", subcategory ="KEGG")%>% 
  dplyr::select(gs_description,gene_symbol)
head(KEGG_df)
hub_gene <- read.csv("../08auc/01.lasso_genes.csv")
hub_gene <- hub_gene$symbol
exp <- readRDS("../00data/tcga_count.rds")
genes2 <- hub_gene

# ii <- "A2M"
# ii <- "CFH"
# tar.exp <- t(exp[ii,])
for (ii in genes2){
  tar.exp <- t(exp[ii,])
  cor_df <- cor(t(exp), tar.exp, method="pearson")
  # print(head(cor_df))
  data1 = data.frame(gene=rownames(cor_df),cor=cor_df[, 1])
  # print(head(data1))
  ge = data1$cor
  names(ge) = data1$gene
  ge = sort(ge,decreasing = T)
  print(head(ge))
  em_kegg <- GSEA(ge, TERM2GENE = KEGG_df,pAdjustMethod = "none")
  print(head(em_kegg))
  
  p1 <- GseaVis::gseaNb(object = em_kegg,
                        # rank.gene = ii,
                        # rank.gene.nudgey= 8,#the gene label nudge y on rank plot, defalut is 2.
                        # addGene = T, #whether add gene name on the curve, defalut is FALSE.
                        # geneSize = 6,# gene label text size, defalut is 4.
                        geneSetID = em_kegg@result$ID[1:5], 
                        addPval = T,
                        pvalX = 1,pvalY = 1,
                        ght.facet = T,
                        termWidth = 80, #the width or the term name, defalut is 40.
                        curveCol = brewer.pal(5,'Paired'))+
    labs(title = ii)+
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18,vjust=55),
          plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
          plot.title.position = "plot")
  
  p1
  write.csv(em_kegg ,file=paste0("GSEA",ii, "_KEGG_gsea.csv"))
  ggsave(paste0("GSEA",ii,'.KEGG.pdf'),width = 8,height = 5, p1) #GSEA分析
  ggsave(paste0("GSEA",ii,'.KEGG.png'),width = 8,height = 5, p1)
}
