rm(list = ls())

foldpath <- paste0("D:/workdir/12stadb/03relatgene")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

options(stringsAsFactors = F)
library(TCGAbiolinks)
library(WGCNA)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(corrplot)
library(pheatmap)
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(readxl)

hlmrg <- read.table("D:/workdir/12stadb/palmitoylation.txt",sep = "\t", header = T)

geneexpr <- readRDS("../00data/tcga_count.rds")
diff_express <- readRDS("../01deg/tcga_filter_diff2.rds")

geneset1 <- hlmrg$Gene
geneset2 <- geneexpr[rownames(diff_express),]
geneset2 <- rownames(geneset2)

tpm <- geneexpr

if(length(intersect(rownames(tpm),geneset1))> 0& length(intersect(rownames(tpm),geneset2))> 0 ){
  
  ##正常组织样本ID
  gs1 <- intersect(rownames(tpm),geneset1)
  gs2 <- intersect(rownames(tpm),geneset2)
  gs2 <- gs2[!gs2 %in% gs1]
  g1exp <- data.frame(t(tpm[gs1,]))
  g2exp <- data.frame(t(tpm[gs2,]))
  
  # method = "pearson"
  cor_cut <- 0.5
  p_cut <- 0.05
  
  for(method in c("pearson","spearman")){
    
    cor_coef <- cor(g2exp,g1exp,method = method)
    # %>% as.data.frame()
    cor_coef <- na.omit(cor_coef)
    cor_p <- corPvalueStudent(cor_coef,nrow(g2exp))
    
    cor_coef <- data.frame(cor_coef)
    
    df <- data.frame(lacgene = character(0), 
                     corgene = character(0), 
                     cor = numeric(0),
                     pvalue  = numeric(0)
    )

    for (i in 1:ncol(cor_coef)){
      for (n in 1:nrow(cor_coef)){
        inputgene <- colnames(cor_coef)[i]
        tesgene <- rownames(cor_coef)[n]
        
        cor <- cor_coef[n,i]
        pvalue <- cor_p[n,i]
        if (abs(cor) > cor_cut & pvalue < p_cut){
          inputdata <- data.frame(lacgene = inputgene,
                                  corgene = tesgene,
                                  cor = cor,
                                  pvalue = pvalue)
          print(paste0(inputgene," ",tesgene," ",cor," ",pvalue))
          df <- rbind(df, inputdata)
        }
      }
    }
    
    write.csv(df, file = paste0(method,".csv"))
  }
}

pgenes <- read.csv("pearson.csv")
sgenes <- read.csv("spearman.csv")

lactate_relatedgenes <- intersect(pgenes$corgene,sgenes$corgene)
lactate_relatedgenes <- unique(lactate_relatedgenes)

cor_coef <- sgenes[sgenes$corgene %in% lactate_relatedgenes,]

table(pgenes$lacgene)
table(sgenes$lacgene)
table(cor_coef$lacgene)

length(unique(cor_coef$corgene))
length(unique(sgenes$corgene))
length(unique(pgenes$corgene))

write.csv(cor_coef,"lacgenes_spearmanvalue.csv")

# 交集基因
hub_gene <- Reduce(intersect, list(pgenes$corgene,sgenes$corgene)) # 4个

which(!geneset1 %in% hub_gene)
hub_gene <- c(hub_gene,geneset1)
write.csv(hub_gene,'hub_gene.csv',row.names = T)
library(VennDiagram)
library(ggvenn)
library(grid)
# # 绘图
mydata<-list('Spearman'= sgenes$corgene,"Pearson"= pgenes$corgene)
pdf('02.venn.pdf',w=6,h=6)
ggvenn(mydata,
       c('Spearman',"Pearson"),
       fill_color = c("#ffb2b2","#7570B3"),
       show_percentage = T,
       fill_alpha = 0.5,
       stroke_alpha = 1,
       stroke_size = 0.4,
       text_size = 5,
       stroke_color=NA,
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#7570B3"),
       set_name_size = 8,
       text_color = 'black'
)
dev.off()

# ------------------------------桑基图-----------------------

rt=cor_coef
rt <- rt[,-1]

mycol=rainbow(length(unique(rt[,"lacgene"])), s=0.8, v=0.8)
# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)

colnames(rt)[1] <- "Gene" 
p1<-ggplot(data=rt, aes(axis1 = corgene , axis2 = Gene, y = 1)) +
  geom_alluvium(aes(fill = Gene), knot.pos = 0.1, reverse = F, width = 0.3)+ 
  geom_stratum(fill=NA, color= NA, width = 0.3,alpha = .9)+
  geom_text(stat = 'stratum', size =1.5, color='black', label.strata = T)+
  scale_fill_manual(values = mycol) +
  scale_x_discrete(limits = c('CorGene','Palmitoylation\nGene      '), expand=c(0, 0))+
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(),axis.text.x = element_blank()) + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  coord_flip()+ggtitle("")   

pdf(file="ggalluvial.pdf", width=5, height=4)
print(p1)
dev.off()

# # ------------------------桑基图2--------------------
# library(networkD3)
# 
# p <- sankeyNetwork(Links = env_corr_list, Nodes = env_group,
#                    Source = 'IDsource', Target = 'IDtarget',  
#                    Value = 'weight', LinkGroup = 'direction', 
#                    NodeID = 'env', NodeGroup = 'group', 
#                    nodePadding = 50, nodeWidth = 20, fontSize = 12, 
#                    colourScale = color, height = 300, width = 800)
# # -----------------------------范癌分析---------------------
# rm(list = ls())
# setwd("D:\\workdir\\03-Z004-S002\\02WGCNA/")
# 
# options(stringsAsFactors = F)
# library(TCGAbiolinks)
# library(WGCNA)
# library(ggplot2)
# library(ggpubr)
# library(ggrepel)
# library(corrplot)
# library(pheatmap)
# # FilePath <- dir("H:/MedBioInfoCloud/analysis/TCGA/new/processedTCGAdata/TCGA-STAR_Exp/",
# # "STARdata.Rdata$",full.names = T)
# 
# geneexpr <- readRDS("../00data/tcga_fpkm.rds")
# diff_express <- readRDS("../01DEG/tcga_filter_diff2.rds")
# 
# geneset1 <- read.csv("../00data/lactategene.csv")
# geneset1 <- geneset1$Target
# # geneset2 <- geneexpr[rownames(diff_express),]
# geneset2 <- rownames(geneexpr)
# 
# ###TCGA数据库中33中癌症类型
# project <- getGDCprojects()$project_id
# project <- project[grep("TCGA-",project)]
# # proj = "TCGA-LUAD"
# pandata <- data.frame()
# 
# for(proj in project){
#   message("===============================")
#   message(proj)
#   load(FilePath[grep(proj,FilePath)])#STARdata
#   tpm <- STARdata[["tpm"]]
#   tpm <- filterGeneTypeExpr(expr = tpm,
#                             fil_col = "gene_type",
#                             filter = FALSE)
#   ##过滤不表达的基因
#   tpm <- tpm[apply(tpm,1,var) !=0,]
#   if(length(intersect(rownames(tpm),geneset1))> 0& length(intersect(rownames(tpm),geneset2))> 0 ){
#     
#     ##正常组织样本ID
#     SamN <- TCGAquery_SampleTypes(barcode = colnames(tpm),
#                                   typesample = c("NT","NB","NBC","NEBV","NBM"))
#     
#     ##肿瘤组织样本ID
#     SamT <- setdiff(colnames(tpm),SamN)
#     
#     ###去除重复样本
#     tur_exp <- del_dup_sample(tpm[,SamT],col_rename = T)
#     ###long2转换
#     tur_exp <- log2(tur_exp + 1)
#     
#     gs1 <- intersect(rownames(tpm),geneset1)
#     gs2 <- intersect(rownames(tpm),geneset2)
#     g1exp <- data.frame(t(tur_exp[gs1,]))
#     g2exp <- data.frame(t(tur_exp[gs2,]))
#     
#     # method = "pearson"
#     for(method in c("pearson","spearman")){
#       
#       
#       cor_coef <- cor(g1exp,g2exp,method = method )
#       cor_p <- corPvalueStudent(cor_coef,nrow(g1exp))
#       
#       col2 <- colorRampPalette(c("#3300CC","#3399FF","white","#FF3333","#CC0000"),alpha = TRUE)
#       
#       # 相关性热图，基因数量太少，不绘制
#       if(min(length(gs1),length(gs2)) >= 3){
#         fp_corrplot <- paste0(opt,proj,"/corrplot/")
#         ifelse(dir.exists(fp_corrplot),"The folder already exists.",dir.create(fp_corrplot,recursive = T))
#         
#         display_numbers = matrix(ifelse(cor_p > 0.01 | is.na(cor_p), "×", ""),
#                                  nrow(cor_p))
#         
#         pdf(paste0(fp_corrplot,method,"-corrplot.pdf"),
#             width = ifelse(length(gs1)>=4,length(gs1)*0.9,3 ),
#             height = ifelse(length(gs2)>=4,length(gs1)*0.9,3 ))
#         corrplot(cor_coef, col = col2(100),method = 'square',cl.length=5,
#                  tl.col="black",tl.cex = 1,cl.pos = "r",cl.ratio = 0.2,
#                  cl.align.text = "l",
#                  p.mat = cor_p, sig.level = c(.001, .01, .05),outline="white",
#                  insig = "label_sig",pch.cex = 1.2, pch.col = "black")%>% print()
#         
#         pheatmap(cor_coef,
#                  #annotation_col =annotation_col,
#                  # annotation_row = annotation_row,
#                  color = col2(50),
#                  cluster_cols =F,
#                  fontsize=6,
#                  cluster_rows = F,
#                  fontsize_row=6,
#                  cellwidth =10,
#                  cellheight =10,
#                  fontsize_number = 16,
#                  display_numbers= display_numbers,
#                  #scale="row",
#                  show_colnames=T,
#                  fontsize_col=8,
#                  main = paste0(method," correlation")) %>% print()
#         
#         dev.off()
#       }
#       
#       ###
#       # g1 = gs1[1]
#       # g2 = gs2[1]
#       fp_point <- paste0(opt,proj,"/point/")
#       ifelse(dir.exists(fp_point),"The folder already exists.",dir.create(fp_point,recursive = T))
#       
#       for(g1 in gs1){
#         for(g2 in gs2){
#           coefficient = round(cor_coef[g1,g2],2)
#           pvalue <- round(cor_p[g1,g2],3)
#           pv <- ifelse(pvalue < 0.001,"p value < 0.001",paste0("p value = ",pvalue))
#           
#           txt <- ifelse(method == "pearson",paste0("Pearson correlation coefficient:",coefficient,"\n",pv),
#                         paste0(paste0("Spearman correlation coefficient:",coefficient,"\n",pv)))
#           data <- data.frame(g1exp[,g1],g2exp[,g2])
#           colnames(data) <- c("x","y")
#           p <- ggplot(data = data, aes(x = x, y = y)) + #数据映射
#             geom_point(alpha = 0.6,shape = 19,size=3,color="#DC143C") +#散点图，alpha就是点的透明度
#             #geom_abline()+
#             labs(title = txt)+
#             geom_smooth(method = lm, formula = y ~ x,aes(colour = "lm"), size = 1.2,se = T)+
#             scale_color_manual(values = c("#808080")) + #手动调颜色c("#DC143C","#00008B", "#808080")
#             theme_bw() +#设定主题
#             theme(axis.title=element_text(size=15,face="plain",color="black"),
#                   axis.text = element_text(size=12,face="plain",color="black"),
#                   legend.position =  "none",
#                   panel.background = element_rect(fill = "transparent",colour = "black"),
#                   plot.background = element_blank(),
#                   plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
#                   legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
#             ylab(paste0("The expression of ",g2,"\nlog2(TPM +1)")) + #expression的作用就是让log10的10下标
#             xlab(paste0("The expression of ",g1,"\nlog2(TPM +1)"))
#           ggsave(paste0(fp_point,g1,"-",g2,"-",method,"-point.pdf"),
#                  plot = p,
#                  width = 5,height = 5)
#           
#           onedata <- data.frame(gene1 = g1,
#                                 gene2 = g2,
#                                 coef = coefficient,p = pvalue,method = method,
#                                 cancer = proj)
#           pandata <- rbind(pandata,onedata)
#         }
#       }
#     }
#     
#   }
# }
# 
# ###==========泛癌=======
# fp_pan <- paste0(opt,"pan-cancer/")
# ifelse(dir.exists(fp_pan),"The folder already exists.",dir.create(fp_pan,recursive = T))
# 
# record_files <- dir(fp_pan,".txt$",full.names = T)
# ifelse(length(record_files)>0,file.remove(record_files),
#        "There is no record file that can be removed")
# 
# 
# cut_ceof <- 0.3
# cut_pval <- 0.01
# method = "pearson"
# for(method in c("pearson","spearman")){
#   subdata <- pandata[pandata$method == method,]
#   # g1 <- unique(subdata$gene1)[1]
#   for(g1 in unique(subdata$gene1)){
#     subdata1 <- subdata[subdata$gene1 == g1,]
#     # g2 <- unique(subdata$gene2)[1]
#     for(g2 in unique(subdata$gene2)){
#       subdata2 <- subdata1[subdata1$gene2 == g2,]
#       
#       subdata2$significance <- ifelse(abs(subdata2$coef)> cut_ceof &subdata2$p < cut_pval,
#                                       ifelse(subdata2$coef > 0,"positive","negative"),"ns")
#       
#       if(length(unique(subdata2$significance)) > 1){
#         labslogi <- abs(subdata2$coef) > cut_ceof & subdata2$p < cut_pval
#         
#         ###筛选要在散点图中展现的癌症类型标签：前3个===============
#         if(length(unique(subdata2$significance)) == 3){
#           cl <- c("#00008B", "#808080","#DC143C")
#           ###正相关的最多显示3个
#           posdat <- subdata2[subdata2$significance == "positive",] %>% arrange(desc(coef))
#           if(length(posdat$significance)> 3){
#             poscancer <- posdat$cancer[1:3]
#           }else{poscancer <- posdat$cancer}
#           ###负相关的最多显示3个
#           negdat <- subdata2[subdata2$significance == "negative",] %>% arrange(coef)
#           if(length(negdat$significance)> 3){
#             negcancer <- negdat$cancer[1:3]
#           }else{negcancer <- negdat$cancer}
#           
#           labscancer <- c(poscancer,negcancer)
#         }else if(sum(subdata2$significance == "positive")!=0){
#           cl <- c("#808080","#DC143C")
#           posdat <- subdata2[subdata2$significance == "positive",] %>% arrange(desc(coef))
#           if(length(posdat$significance)> 3){
#             labscancer <- posdat$cancer[1:3]
#           }else{labscancer <- posdat$cancer}
#           
#         }else{
#           cl <- c("#DC143C", "#808080")
#           negdat <- subdata2[subdata2$significance == "negative",] %>% arrange(coef)
#           if(length(negdat$significance)> 3){
#             labscancer <- negdat$cancer[1:3]
#           }else{labscancer <- negdat$cancer}
#         }
#         subdata2$labs <- ""
#         subdata2$labs[match(labscancer,subdata2$cancer)] <- labscancer
#         subdata2$labs <- gsub("TCGA-","",subdata2$labs)
#         ###^^^^^^^^^^^^^^^^^^^
#         ###如果p值为0.随机生成极小p值
#         subdata2$nlog10p <- -log10(subdata2$p)
#         if(sum(subdata2$p == 0)>0){
#           nInf <- length(grep("Inf",subdata2$nlog10p))
#           subdata2$nlog10p[grep("Inf",subdata2$nlog10p)] <- -log10(runif(nInf,0.000000000001,0.0000001) )
#         }
#         
#         p <-  ggplot(data = subdata2, aes(x = coef,
#                                           y = nlog10p,
#                                           colour = significance)) + #数据映射
#           geom_point(alpha = 0.9,shape = 19,size=3) +#散点图，alpha就是点的透明度
#           scale_color_manual(values = cl) +
#           theme_bw() +#设定主题
#           geom_text_repel(label = subdata2$labs,
#                           size = 5,
#                           segment.color = "black", #连接线的颜色，就是名字和点之间的线
#                           show.legend = FALSE) +
#           labs(title = paste0())+
#           theme(axis.title=element_text(size=15,face="plain",color="black"),
#                 axis.text = element_text(size=12,face="plain",color="black"),
#                 legend.position =  "none",
#                 panel.background = element_rect(fill = "transparent",colour = "black"),
#                 plot.background = element_blank(),
#                 plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
#                 legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
#           ylab(expression(-log[10]("p.value"))) +#expression的作用就是让log10的10下标
#           xlab(paste0(method," correlation coefficient of\ngenes ",g1, " and ", g2))+
#           geom_vline(xintercept = c(-cut_ceof,cut_ceof),
#                      lty = 2,
#                      col = "black",
#                      lwd = 0.3) +
#           geom_hline(yintercept = -log10(cut_pval),
#                      lty = 2,
#                      col = "black",
#                      lwd = 0.3)
#         # p
#         ggsave(filename = paste0(fp_pan,g1,"-",g2,"-",method,"_in_pan-cancer.pdf"),
#                plot = p,height=5,width=5)
#         
#         ##记录分析的样本信息,只记录一次
#         pos <- subdata2$cancer[subdata2$significance == "positive"]
#         neg <- subdata2$cancer[subdata2$significance == "negative"]
#         opinfo <- paste0(method,":",g1," and ",g2,":")
#         output <- file(paste0(fp_pan,"00-RecordInformation.txt"), open="a+b")
#         writeLines(c(opinfo,"\tpositive:",paste0("\t\t",pos),"\tnegative:",paste0("\t",neg)),
#                    con=output)
#         close(output)
#       }
#     }
#   }
# }
