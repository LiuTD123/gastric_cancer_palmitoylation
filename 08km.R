rm(list = ls())

foldpath <- "D:/workdir/12stadb/08km"
if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(survival)
library(survminer)
library(stringr)
load("../08auc/riskscore.RData")

for (i in 1:length(rs)){
  # 首先分组
  title <- names(rs[i])
  i = 1
  cut <- survminer::surv_cutpoint(rs[[i]], #数据集
                                  minprop = 0.25,
                                  time = "OS.time", #生存时间
                                  event = "OS", #生存状态
                                  variables = "RS"  #需要计算的数据列名
  )
  cut
  
  rs[[i]]$group <- ifelse(rs[[i]][,"RS"] > summary(cut)[1,1],'High','Low' )
  
  # km拟合
  kmfit <- survfit(Surv(OS.time, OS) ~ group, data = rs[[i]])
  
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
  
  palette_colors <- c("#f16c23","#2b6a99")
  
  train_km <- ggsurvplot(kmfit,
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
  
  pdf_file_name <- paste0(title, ".pdf")
  pdf(pdf_file_name, family = "Times", height = 6, width = 6, onefile = F)
  print(train_km)
  dev.off()
  
  # ----------------风险评分分布----------------
  res <- survminer::surv_cutpoint(rs[[i]], #数据集
                                  minprop = 0.25,
                                  time = "OS.time", #生存时间
                                  event = "OS", #生存状态
                                  variables = "RS"  #需要计算的数据列名
  )
  risk_dis <- ggplot(rs[[i]], aes(x=reorder(rs[[i]]$ID, rs[[i]]$RS), y=rs[[i]]$RS, color = rs[[i]]$group)) +
    geom_point() +
    scale_color_manual(values = palette_colors) + 
    scale_x_discrete(breaks = rownames(rs[[i]]$ID)[order(rs[[i]]$RS)][c(1,50,100,150,200,250,300,350,400,450)],
                     labels = c(1,50,100,150,200,250,300,350,400,450),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(rs[[i]][which(rs[[i]]$group=="Low"),]) + 0.5, lty = 2) +
    # geom_hline(yintercept = summary(res)[1,1], lty =2) +
    labs(x = "Patients(increasing risk score)",
         y = "Risk Score",
         title = paste0(title," Risk Score Distribution")) + 
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          # legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15),
          text = element_text(family = "Times", face = "bold"))
  ggsave(filename = paste0(title,"_riskScore_dis.pdf"), height = 3, width = 5, risk_dis)
  
  # --------------------
  res <- survminer::surv_cutpoint(rs[[i]], #数据集
                                  minprop = 0.25,
                                  time = "OS.time", #生存时间
                                  event = "OS", #生存状态
                                  variables = "RS"  #需要计算的数据列名
  )
  
  surv_stat <- ggplot(rs[[i]], aes(x=reorder(rs[[i]]$ID, rs[[i]]$RS),
                                   y=OS.time,
                                   color = factor(OS,
                                                  levels = c(1,0),
                                                  labels = c("Dead","Alive")))) +
    geom_point() +
    scale_color_manual(values = palette_colors) +
    scale_x_discrete(breaks = rs[[i]]$ID[order(rs[[i]]$RS)][c(1,50,100,150,200,250,300,350,400,450)],
                     labels = c(1,50,100,150,200,250,300,350,400,450),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(rs[[i]][which(rs[[i]]$group=="Low"),]) + 0.5, lty = 2) +
    labs(x = "Patients(increasing risk score)",
         y = paste0("Overall Survival (days)"),
         title = paste0(title," Distribution")) + 
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          # legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15),
          text = element_text(family = "Times", face = "bold"))
  ggsave(filename = paste0(title,"_dis.pdf"), height = 3, width = 5, surv_stat)
}
