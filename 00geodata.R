rm(list = ls())

foldpath <- paste0("D:/workdir/10stad/00geodata")

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(GEOquery)
library(idmap3)
options(timeout = 99999999)
library(readxl)

# ----------------------------------GSE13861------------------
rm(list = ls())
GSEID <- "GSE13861"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

table(cli$characteristics_ch1)
cli <- cli[cli$characteristics_ch1 == "gastric adenocarcinoma",]

GSE13861 <- as.data.frame(t(exp_symbol))
GSE13861 <- GSE13861[rownames(cli),]
GSE13861$ID <- rownames(GSE13861)

idcard <- cli[,c(2,1)]
surv <- read_excel("GSE13861_outcome_deal.xlsx")
idcard$title <- substr(idcard$title,12,10000)
idcard$title <- substr(idcard$title,1,5)
surv <- surv[,c(1,10,9)]
colnames(idcard) <- c("ID","Patients_ID")
surv <- merge(idcard,surv,by = "Patients_ID")
surv <- surv[,-1]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- round(surv$OS.time*30)

GSE13861 <- merge(surv,GSE13861,by = "ID")

save(GSE13861,file = "GSE13861.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE28541------------------------------
rm(list = ls())
GSEID <- "GSE28541"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,31)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

GSE28541 <- as.data.frame(t(exp_symbol))
GSE28541$ID <- rownames(GSE28541)
surv <- read_excel("GSE28541_outcome_deal.xlsx")
surv <- surv[,c(1,8,7)]
colnames(surv) <- c("id","OS.time","OS")
surv$OS.time <- round(surv$OS.time*30)
idcard <- cli[,c(1,2)]
colnames(idcard) <- c("id","ID")
surv <- merge(idcard,surv,by = "id")
surv <- surv[,-1]
GSE28541 <- merge(surv,GSE28541,by="ID")

save(GSE28541,file = "GSE28541.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE15459------------------------------
rm(list = ls())
GSEID <- "GSE15459"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

GSE15459 <- as.data.frame(t(exp_symbol))
GSE15459$ID <- rownames(GSE15459)
surv <- read_excel("GSE15459_outcome.xls")
surv <- surv[,c(1,9,10)]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- round(surv$OS.time*30)
GSE15459 <- merge(surv, GSE15459,by = "ID")

save(GSE15459,file = "GSE15459.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE84437------------------------------
rm(list = ls())
GSEID <- "GSE84437"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol <- na.omit(exp_symbol)
# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

GSE84437 <- exp_symbol
GSE84437 <- as.data.frame(t(GSE84437))
GSE84437$ID <- rownames(GSE84437)

surv <- cli[,c("geo_accession","death:ch1","duration overall survival:ch1")]
colnames(surv) <- c("ID","OS","OS.time")
GSE84437 <- merge(surv,GSE84437, by="ID")
GSE84437$OS.time <- as.numeric(GSE84437$OS.time)*30
GSE84437$OS <- as.numeric(GSE84437$OS)
GSE84437 <- GSE84437[,c(1,3,2,4:ncol(GSE84437))]

save(GSE84437,file = "GSE84437.RData")

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE13911------------------------------
rm(list = ls())
GSEID <- "GSE13911"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE29272------------------------------
rm(list = ls())
GSEID <- "GSE29272"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]
exp_symbol <- exp_symbol[,rownames(cli[cli$source_name_ch1 == "tumor tissue",])]
GSE29272 <- as.data.frame(t(exp_symbol))
GSE29272$ID <- rownames(GSE29272)

surv <- read_excel("GSE29272_outcome_deal.xlsx")

idcard <- cli[,c("geo_accession","title")]
idcard <- idcard[rownames(cli[cli$source_name_ch1 == "tumor tissue",]),]
library("stringr")     
idcard$title <- str_sub(idcard$title, -8,-2)
colnames(idcard) <- c("geo","ID")
surv <- merge(idcard,surv,by ="ID")
surv <- surv[,c(2,14,12)]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- as.numeric(surv$OS.time)
surv$OS <- ifelse(surv$OS == "Y",0,1)
surv <- na.omit(surv)
GSE29272 <- merge(surv,GSE29272,by = "ID")

save(GSE29272,file = "GSE29272.RData")

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE79973------------------------------
rm(list = ls())
GSEID <- "GSE79973"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE118916------------------------------
rm(list = ls())
options(timeout = 999999999)
GSEID <- "GSE118916"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE66229------------------------------
rm(list = ls())
options(timeout = 999999999)
GSEID <- "GSE66229"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

cli <- cli[cli$`tissue:ch1`=="Gastric tumor",]
exp_symbol <- exp_symbol[,rownames(cli)]
surv <- read_excel("GSE66229_outcome_deal.xlsx")
surv <- surv[,c(1,7,8)]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- round(surv$OS.time*30)
GSE66229 <- as.data.frame(t(exp_symbol))
GSE66229$ID <- rownames(GSE66229)
GSE66229 <- merge(surv,GSE66229,by = "ID")

save(GSE66229,file = "GSE66229.RData")

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE62254------------------------------
rm(list = ls())
options(timeout = 999999999)
GSEID <- "GSE62254"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE26942------------------------------
rm(list = ls())
options(timeout = 999999999)
GSEID <- "GSE26942"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

table(cli$`tissue:ch1`)
cli <- cli[cli$`tissue:ch1` %in% c("Gastric cancer tissue","Gastric tumor tissue"),]
GSE26942 <- as.data.frame(t(exp_symbol))
GSE26942 <- GSE26942[rownames(cli),]
GSE26942$ID <- rownames(GSE26942)
idcard <- cli[,c(2,46)]
surv <- read_excel("GSE26942_outcome_deal.xlsx")
surv <- surv[,c(1,10,9)]
colnames(surv) <- c("ID","OS.time","OS")
surv$OS.time <- round(surv$OS.time*30)
colnames(idcard) <- c("geo","ID")
surv <- merge(surv,idcard,by = "ID")
surv <- surv[,-1]
surv <- surv[,c(3,1,2)]
colnames(surv)[1] <- c("ID")

GSE26942 <- merge(surv,GSE26942,by = "ID")

save(GSE26942,file = "GSE26942.RData")
write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE29998------------------------------
rm(list = ls())
options(timeout = 999999999)
GSEID <- "GSE29998"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))

# ----------------------------------GSE38024------------------------------
rm(list = ls())
options(timeout = 999999999)
GSEID <- "GSE38024"
gset <- getGEO(filename = paste0(GSEID,"_series_matrix.txt.gz"),AnnotGPL = T, destdir = "./")
# gset <- getGEO("GSE106090",AnnotGPL = T, destdir = "./",getGPL = T)

gset[1]

# 提取表达量
exp <- exprs(gset)

# 临床信息
cli <- pData(gset)  #获取临床信息
save(cli,file = paste0(GSEID,"_cliall.RData"))

head(cli)

gpl<- fData(gset)

# 执行下面代码后，便通过只取第一个Gene symbol来达到去重效果，这时再打开gpl发现已经没有对应多个Gene symbol的情况了
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

# 接下来对表达矩阵exp进行探针ID转化，首先将其转换为数据框的形式，否则一会会报错
exp<-as.data.frame(exp)

# 此时虽然我们能看到exp的探针ID号（最左侧），但这个ID号并不是真正在列中，所以我们还需要将探针ID放到数据框里面，执行下面命令
exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息

gpl<-gpl[,c(1,3)]
# 接下来通过将gpl合并到exp中，从而在exp中将芯片ID号转化为Gene symbol
exp_symbol<-merge(exp,gpl,by="ID")

# 首先查看有多少重复的Gene symbol
table(duplicated(exp_symbol$`ID`))
# 输入以下命令进行去重操作同时将矩阵转换成标准的表达矩阵形式
# exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`ID`)

exp_symbol <- exp_symbol[,-1]

# 删除重复的
exp_symbol <- exp_symbol[!duplicated(exp_symbol$`Gene symbol`),]
exp_symbol <- na.omit(exp_symbol)
rownames(exp_symbol) <- exp_symbol$`Gene symbol`

exp_symbol <- exp_symbol[,-ncol(exp_symbol)]

write.csv(exp_symbol,paste0("dat_",GSEID,".csv"))
save(exp_symbol, file = paste0("dat_",GSEID,".RData"))