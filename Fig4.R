#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.4
# Author: Chao Cheng; Chenyang Li
# Note: CPRIT-MIRA is the MDAMPLC cohort
# 12/15/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx ##################################################################

# [1.1] infer imm infiltration =================================================
rm(list=ls())

outdir <- "./Fig4/"

myinf1 <- "~/Mydata/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- "./BASE/ImmGen_representative_TILs_Profile_4BASE.txt"

myfun1 <- "./BASE/base5.R"

myoutf1 <- paste0(outdir,"/RNAseq_Genes_representativeImm_TRACERx.txt")

data <- read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)

mywt <- read.table(myinf2, sep="\t", header=T, row.names=1)

source(myfun1)
#a) Remove low exp and log transform -------------------------------------------
xx <- apply(data>0, 1, sum)
se <- which(xx>0)
data <- data[se,]
dim(data)
range(data)

data <- log10(data+1)

# BASE  ------------------------------------------------------------------------
reg <- mywt
xx <- base5(data, reg, perm=1000, myoutf1, median.norm=T)


# [1.2] Avg/Max/Min immune cell score ==========================================
rm(list=ls())

outdir <- "./Fig4/"

myinf1 <- paste0(outdir,"/RNAseq_Genes_representativeImm_TRACERx.txt")

myoutf1 <- paste0(outdir,"/Imm_scores_TRACERx.rda")


# a) Immune score --------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1)
cnum <- ncol(data)/4
data <- data[, 1:cnum] - data[, (cnum+1):(2*cnum)]
tmp <- colnames(data)
tmp <- gsub("_up\\.ES", "", tmp)
colnames(data) <- tmp

# b) avg/max/min ---------------------------------------------------------------
p.tag <- substr(row.names(data),17, 24)
pid <- unique(p.tag)

res <- matrix(0, length(pid), ncol(data))
row.names(res) <- pid
colnames(res) <- colnames(data)
avg.mat <- max.mat <- min.mat <- res

for(k in 1:nrow(res))
{
  cat("\r", k)
  se <- which(p.tag==pid[k])
  avg.mat[k,] <- apply(data[se, ], 2, mean)
  max.mat[k,] <- apply(data[se, ], 2, max)
  min.mat[k,] <- apply(data[se, ], 2, min)
}

save(avg.mat, max.mat, min.mat, file= myoutf1)
# [1.3] Survival Analysis ======================================================
rm(list=ls())
library(dplyr)
library(ggplot2)
library(survival)
library(forestplot)
outdir <- "./Fig4/"

myinf1 <- paste0(outdir,"/Imm_scores_TRACERx.rda")
myinf2 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myFig1 <- paste0(outdir,"/Imm_scores_forest_TRACERx.pdf")
myFig2 <- paste0(outdir,"/Imm_scores_CIbar_TRACERx.pdf")
# a) expression and survival data ----------------------------------------------
load(myinf1)
ls()

info <- read.csv(myinf2, header=T, row.names=1)

comxx <- intersect(row.names(avg.mat), row.names(info))
info <- info[comxx,]
avg.mat <- avg.mat[comxx,]
max.mat <- max.mat[comxx,]
min.mat <- min.mat[comxx,]

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months. 
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)



# b) immune cell survival analysis result --------------------------------------
res <- data.frame(matrix(ncol = 0, nrow = 6))
for (i in c("avg","max","min")){
  nam1 <- paste0(i,".mat")
  data <- get(nam1)
  upper <- lower <- ci <- hr <- pv <- rep(0, ncol(data))
  for(k in 1:ncol(data))
  {
    cat("\r", k)
    mytf <- as.numeric(data[,k])
    xx <- cbind(mytf, info)
    xx <- xx[xx[, "t.surv"]>0,]
    mycox <- coxph(Surv(t.surv, e.surv)~mytf, xx) 
    mycox <- summary(mycox)
    pv[k] <- round(mycox$coefficients[5],4)
    tmp <- mycox$conf.int
    hr[k] <- round(tmp[1],3)
    upper[k] <- round(tmp[4],3)
    lower[k] <- round(tmp[3],3)
    ci[k] <- round(mycox$concordance[1],3)
  }
  
  tmp <- data.frame("HR"=hr,
                    "P-Value"=pv,
                    "C-index"=ci, 
                    "lower.95"=lower, 
                    "upper.95"=upper)
  colnames(tmp) <- paste0(colnames(tmp),".",i)
  row.names(tmp) <- colnames(data)
  res <- cbind(res,tmp)
}
res

# c) plot data  ----------------------------------------------------------------
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

p.data <- data.frame(matrix(ncol = 0, nrow = 6))
for (i in c("avg","max","min")){
  nam1 <- paste0(i,".mat")
  data <- get(nam1)
  upper <- lower <- ci <-  hr <- pv <- rep(0, ncol(data))
  for(k in 1:ncol(data))
  {
    cat("\r", k)
    mytf <- as.numeric(data[,k])
    xx <- cbind(mytf, info)
    xx <- xx[xx[, "t.surv"]>0,]
    mycox <- coxph(Surv(t.surv, e.surv)~mytf, xx) 
    mycox <- summary(mycox)
    pv[k] <- round(mycox$coefficients[5],4)
    tmp <- mycox$conf.int
    hr[k] <- round(tmp[1],3)
    upper[k] <- round(tmp[4],3)
    lower[k] <- round(tmp[3],3)
    ci[k] <- round(mycox$concordance[1],3)
  }
  
  tmp <- data.frame(Method =rep(paste0(CapStr(i),"."),length(hr)),
                   Cell = c( "Naive B", "Memory B",
                             "CD8+ T", "CD4+ T", "NK cell", "Monocyte"),
                   "HR"=hr,
                   "P-Value"=pv,
                   "C-index"=ci, 
                   "lower.95"=lower, 
                   "upper.95"=upper)
  p.data <- rbind(p.data,tmp)
}
p.data

# d) forest plot ---------------------------------------------------------------
p.data <- p.data %>% 
  mutate(Cell = factor(Cell, levels = c( "Naive B", "Memory B",
                                         "CD8+ T", "CD4+ T", 
                                         "NK cell", "Monocyte"))) %>%
  arrange(Cell) %>% 
  mutate(Cell = as.character(Cell))
tabletext <- cbind(
                   c("Cell Type",p.data$Cell),
                   c("Method",p.data$Method),
                   c("Hazard Ratio",p.data$HR),
                   c("P-value",p.data$P.Value),
                   c("C-index",p.data$C.index))
pdf(myFig1, width = 6,height = 5)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="Immune Cell Infiltration",
  xlab="Hazard Ratio",
  grid = T,
  new_page = getOption("forestplot_new_page", F),
  txt_gp = fpTxtGp(label = gpar(fontfamily = "serif", cex = 1),
                   ticks = gpar(fontfamily = "serif", cex = 1),
                   xlab  = gpar(fontfamily = "serif", cex = 1),
                   title = gpar(fontfamily = "serif", cex = 1)),
  col=fpColors(box="#d7302790" , lines="#d7302790" , zero = "gray"),
  zero=1, cex=0.9, lineheight = "auto", boxsize=0.3, colgap=unit(3,"mm"),
  lwd.ci=1, lty.ci = 2, ci.vertices=TRUE, ci.vertices.height = 0.2,
  hrzl_lines = list("1" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044"),
                    "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                    "5" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "8" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "11" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "14" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "17" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "20" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()
# e) bar plot ---------------------------------------------------------------
p1 <-  p.data %>% 
  mutate(Cell = factor(Cell, levels = c( "Naive B", "Memory B",
                                         "CD8+ T", "CD4+ T", 
                                         "NK cell", "Monocyte"))) %>%
  mutate(Method = factor(Method, levels = c("Avg.", "Max.", "Min." ))) %>%
  arrange(desc(Cell),desc(Method)) %>%
  mutate(label = paste0(Cell,"_",Method)) %>%
  mutate(label = factor(label,levels = label)) %>% 
  ggplot(mapping = aes(x = label,y=C.index-0.5)) +
  geom_hline(lty=2,size = 0.9,yintercept = 0.1,color="gray") +
  geom_hline(lty=2,size = 0.9,yintercept = 0.2,color="gray") +
  geom_bar(mapping = aes(fill=Method),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  coord_flip() +
  labs(x="Cell Type",y="C-index",
       title=paste0("TRACERx"),
       base_size = 10, base_family = "serif",face="bold") +
  theme(
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",
                                family = "serif", size=5),
    legend.text= element_text(face="bold", color="black",
                              family = "serif", size=5),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=8),
    axis.text.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=10),
    axis.title.y = element_text(face="bold",color="black", size=10)) +
  scale_y_continuous(limits=c(0, 0.22), breaks=seq(0, 0.5, 0.1),
                     labels = seq(0.5, 1, 0.1), expand = c(0,0)) 
pdf(myFig2, width = 3,height = 5)
print(p1)
dev.off()


# [2] CPRIT-MIRA ##################################################################

# [2.1] infer imm infiltration =================================================
rm(list=ls())

outdir <- "./Fig4/"

myinf1 <- "./data/MDAMPLC.rda"
myinf2 <- "./BASE/ImmGen_representative_TILs_Profile_4BASE.txt"

myfun1 <- "./BASE/base5.R"

myoutf1 <- paste0(outdir,"/RNAseq_Genes_representativeImm_CPRIT-MIRA.txt")

data <- read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)

mywt <- read.table(myinf2, sep="\t", header=T, row.names=1)

source(myfun1)
#a) Remove low exp and log transform -------------------------------------------
load(myinf1)
data <- RNAseq

# remove low exp
xx <- apply(data>0, 1, sum)
se <- which(xx>0)
data <- data[se,]
dim(data)
range(data)

# log tranform
data <- log10(data+1)

# BASE  ------------------------------------------------------------------------
reg <- mywt
xx <- base5(data, reg, perm=1000, myoutf1, median.norm=T)


# [2.2] Avg/Max/Min immune cell score ==========================================
rm(list=ls())

outdir <- "./Fig4/"

myinf1 <- paste0(outdir,"/RNAseq_Genes_representativeImm_CPRIT-MIRA.txt")

myoutf1 <- paste0(outdir,"/Imm_scores_CPRIT-MIRA.rda")


# a) Immune score --------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1)
cnum <- ncol(data)/4
data <- data[, 1:cnum] - data[, (cnum+1):(2*cnum)]
tmp <- colnames(data)
tmp <- gsub("_up\\.ES", "", tmp)
colnames(data) <- tmp

# b) avg/max/min ---------------------------------------------------------------
p.tag <- substr(row.names(data),1, 6)
pid <- unique(p.tag)

res <- matrix(0, length(pid), ncol(data))
row.names(res) <- pid
colnames(res) <- colnames(data)
avg.mat <- max.mat <- min.mat <- res

for(k in 1:nrow(res))
{
  cat("\r", k)
  se <- which(p.tag==pid[k])
  avg.mat[k,] <- apply(data[se, ], 2, mean)
  max.mat[k,] <- apply(data[se, ], 2, max)
  min.mat[k,] <- apply(data[se, ], 2, min)
}

save(avg.mat, max.mat, min.mat, file= myoutf1)
# [2.3] Survival Analysis ======================================================
rm(list=ls())
library(dplyr)
library(ggplot2)
library(survival)
library(forestplot)
outdir <- "./Fig4/"

myinf1 <- paste0(outdir,"/Imm_scores_CPRIT-MIRA.rda")
myinf2 <-  "./data/MDAMPLC.rda"

myFig1 <- paste0(outdir,"/Imm_scores_forest_CPRIT-MIRA.pdf")
myFig2 <- paste0(outdir,"/Imm_scores_CIbar_CPRIT-MIRA.pdf")
myFig3 <- paste0(outdir,"/Imm_scores_CIbarFacet_CPRIT-MIRA.pdf")
# a) expression and survival data ----------------------------------------------
load(myinf2)

info <- Clinical.info
colnames(info) <-  c("t.surv","e.surv", 
                     "Smoker.current.never.former",
                     "Gender.Male.Female" ,         
                     "Diagnosis.Age"  ,
                     "Stage" ,
                     "Subtype"         )


comxx <- intersect(row.names(avg.mat), row.names(info))
info <- info[comxx,]
avg.mat <- avg.mat[comxx,]
max.mat <- max.mat[comxx,]
min.mat <- min.mat[comxx,]


# b) immune cell survival analysis result --------------------------------------
res <- data.frame(matrix(ncol = 0, nrow = 6))
for (i in c("avg","max","min")){
  nam1 <- paste0(i,".mat")
  data <- get(nam1)
  upper <- lower <- ci <- hr <- pv <- rep(0, ncol(data))
  for(k in 1:ncol(data))
  {
    cat("\r", k)
    mytf <- as.numeric(data[,k])
    xx <- cbind(mytf, info)
    xx <- xx[xx[, "t.surv"]>0,]
    mycox <- coxph(Surv(t.surv, e.surv)~mytf, xx) 
    mycox <- summary(mycox)
    pv[k] <- round(mycox$coefficients[5],4)
    tmp <- mycox$conf.int
    hr[k] <- round(tmp[1],3)
    upper[k] <- round(tmp[4],3)
    lower[k] <- round(tmp[3],3)
    ci[k] <- round(mycox$concordance[1],3)
  }
  
  tmp <- data.frame("HR"=hr,
                   "P-Value"=pv,
                   "C-index"=ci, 
                   "lower.95"=lower, 
                   "upper.95"=upper)
  colnames(tmp) <- paste0(colnames(tmp),".",i)
  row.names(tmp) <- colnames(data)
  res <- cbind(res,tmp)
}
res

# c) plot data  ----------------------------------------------------------------
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

p.data <- data.frame(matrix(ncol = 0, nrow = 6))
for (i in c("avg","max","min")){
  nam1 <- paste0(i,".mat")
  data <- get(nam1)
  upper <- lower <- ci <- hr <- pv <- rep(0, ncol(data))
  for(k in 1:ncol(data))
  {
    cat("\r", k)
    mytf <- as.numeric(data[,k])
    xx <- cbind(mytf, info)
    xx <- xx[xx[, "t.surv"]>0,]
    mycox <- coxph(Surv(t.surv, e.surv)~mytf, xx) 
    mycox <- summary(mycox)
    pv[k] <- round(mycox$coefficients[5],4)
    tmp <- mycox$conf.int
    hr[k] <- round(tmp[1],3)
    upper[k] <- round(tmp[4],3)
    lower[k] <- round(tmp[3],3)
    ci[k] <- round(mycox$concordance[1],3)
  }
  
  tmp <- data.frame(Method =rep(paste0(CapStr(i),"."),length(hr)),
                   Cell = c( "Naive B", "Memory B",
                             "CD8+ T", "CD4+ T", "NK cell", "Monocyte"),
                   "HR"=hr,
                   "P-Value"=pv,
                   "C-index"=ci, 
                   "lower.95"=lower, 
                   "upper.95"=upper)
  p.data <- rbind(p.data,tmp)
}
p.data

# d) forest plot ---------------------------------------------------------------
p.data <- p.data %>% 
  mutate(Cell = factor(Cell, levels = c( "Naive B", "Memory B",
                                         "CD8+ T", "CD4+ T", 
                                         "NK cell", "Monocyte"))) %>%
  arrange(Cell) %>% 
  mutate(Cell = as.character(Cell))
tabletext <- cbind(
  c("Cell Type",p.data$Cell),
  c("Method",p.data$Method),
  c("Hazard Ratio",p.data$HR),
  c("P-value",p.data$P.Value),
  c("C-index",p.data$C.index))
pdf(myFig1, width = 6,height = 5)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="Immune Cell Infiltration",
  xlab="Hazard Ratio",
  grid = T,
  new_page = getOption("forestplot_new_page", F),
  txt_gp = fpTxtGp(label = gpar(fontfamily = "serif", cex = 1),
                   ticks = gpar(fontfamily = "serif", cex = 1),
                   xlab  = gpar(fontfamily = "serif", cex = 1),
                   title = gpar(fontfamily = "serif", cex = 1)),
  col=fpColors(box="#3288bd90" , lines="#3288bd90" , zero = "gray"),
  zero=1, cex=0.9, lineheight = "auto", boxsize=0.3, colgap=unit(3,"mm"),
  lwd.ci=1, lty.ci = 2, ci.vertices=TRUE, ci.vertices.height = 0.2,
  hrzl_lines = list("1" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044"),
                    "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                    "5" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "8" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "11" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "14" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "17" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "20" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()
# e) bar plot ---------------------------------------------------------------
p1 <-  p.data %>% 
  mutate(Cell = factor(Cell, levels = c( "Naive B", "Memory B",
                                         "CD8+ T", "CD4+ T", 
                                         "NK cell", "Monocyte"))) %>%
  mutate(Method = factor(Method, levels = c("Avg.", "Max.", "Min." ))) %>%
  arrange(desc(Cell),desc(Method)) %>%
  mutate(label = paste0(Cell,"_",Method)) %>%
  mutate(label = factor(label,levels = label)) %>% 
  ggplot(mapping = aes(x = label,y=C.index-0.5)) +
  geom_hline(lty=2,size = 0.9,yintercept = 0.1,color="gray") +
  geom_hline(lty=2,size = 0.9,yintercept = 0.2,color="gray") +
  geom_bar(mapping = aes(fill=Method),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  coord_flip() +
  labs(x="Cell Type",y="C-index",
       title=paste0("CPRIT-MIRA"),
       base_size = 10, base_family = "serif",face="bold") +
  theme(
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",
                                family = "serif", size=5),
    legend.text= element_text(face="bold", color="black",
                              family = "serif", size=5),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=8),
    axis.text.y = element_text(face="bold", color="black", size=8),
    axis.title.x = element_text(face="bold", color="black", size=10),
    axis.title.y = element_text(face="bold",color="black", size=10)) +
  scale_y_continuous(limits=c(0, 0.22), breaks=seq(0, 0.5, 0.1),
                     labels = seq(0.5, 1, 0.1), expand = c(0,0)) 
pdf(myFig2, width = 3,height = 5)
print(p1)
dev.off()


# f) bar facet  ----------------------------------------------------------------
p2 <-  p.data %>% 
  mutate(Cell = factor(Cell, levels = c( "Naive B", "Memory B",
                                         "CD8+ T", "CD4+ T", 
                                         "NK cell", "Monocyte"))) %>%
  mutate(Method = factor(Method, levels = c("Avg.", "Max.", "Min." ))) %>%
  
  ggplot(mapping = aes(x = Method,y=C.index-0.5)) +
  geom_bar(mapping = aes(fill=Method),
           width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  facet_wrap( ~ Cell, nrow=1) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Cell Type",y="C-index",
       title=paste0("CPRIT-MIRA"),
       base_size = 18, base_family = "serif",face="bold") +
  theme(
    legend.position="top",
    legend.box = "horizontal",
    legend.direction= "horizontal",
    panel.grid=element_blank(),
    legend.key.width = unit(0.5,"cm"),
    legend.title = element_text(face="bold", color="black",family = "serif", size=10),
    legend.text= element_text(face="bold", color="black",family = "serif", size=10),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face="bold", color="black", size=15, 
                               angle = 60,  vjust = 1, hjust=1),
    axis.text.y = element_text(face="bold", color="black", size=15),
    axis.title.x = element_text(face="bold", color="black", size=18),
    axis.title.y = element_text(face="bold",color="black", size=18)) +
  scale_y_continuous(limits=c(0, 0.18), breaks=seq(0, 0.5, 0.05),
                     labels = seq(0.5, 1, 0.05), expand = c(0,0)) 
pdf(myFig3, width = 8,height = 6)
print(p2)
dev.off()
