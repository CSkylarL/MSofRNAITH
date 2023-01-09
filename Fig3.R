#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.3
# Author: Chao Cheng; Chenyang Li
# CPRIT-MIRA  Result is not significant due to the limit sample size,
#             But the cox regression and C-index have the same trend 
#             --> show C-index
# 12/14/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx ##################################################################

# [1.1] ORACLE Score ===========================================================
rm(list=ls())
library(survival)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- paste0(outdir,"/ORACLE.txt")
myinf3 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myoutf1 <- paste0(outdir,"/ORACLE_scores_TRACERx.txt")
myoutf2 <- paste0(outdir,"/ORACLE_M2_Region_TRACERx.txt")

# a) expression and signature  -------------------------------------------------

data <- read.table(myinf1, sep="\t", header=T, check.names=F)
oracle <- read.table(myinf2,sep="\t",header = T)

# confirm the same order of gene
row.names(oracle) <- oracle$Gene.Symbol
comxx <- intersect(row.names(data), row.names(oracle))
data <- data[comxx,]
oracle <- oracle[comxx,]
all(rownames(data) == oracle$Gene.Symbol) # T


oracle_coefficient <- oracle$Model.Coefficient
names(oracle_coefficient) <- oracle$Gene.Symbol
# b) log transform -------------------------------------------------------------

data <- log2(data+1)
# c) survival info -------------------------------------------------------------
info <- read.csv(myinf3, header=T, row.names=1)

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 17, 25)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23   161

# d) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
all(row.names(tmp)  == names(oracle_coefficient)) # T
colnames(tmp) <- mypat
dat.rev <- dat.adj <- dat.avg <- dat.max <- dat.min <- tmp
# patient gene expression 
for(k in 1:length(mypat)){
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  tmp <- data[, se, drop=F]
  dat.max[,k] <- apply(tmp, 1, max)			
  dat.min[,k] <- apply(tmp, 1, min)			
  dat.avg[,k] <- apply(tmp, 1, mean)		
  dat.adj[,k] <- ifelse(test = oracle_coefficient > 0, 
                       yes = dat.max[,k],  no = dat.min[,k])
  dat.rev[,k] <- ifelse(test = oracle_coefficient < 0, 
                       yes = dat.max[,k],  no = dat.min[,k])
}

# patient score 
res.m1 <- data.frame(m1.s.avg=t(dat.avg) %*% oracle_coefficient, 
                    m1.s.max=t(dat.max) %*% oracle_coefficient, 
                    m1.s.min=t(dat.min) %*% oracle_coefficient,
                    m1.s.adj=t(dat.adj) %*% oracle_coefficient,
                    m1.s.rev=t(dat.rev) %*% oracle_coefficient )

# e) M2 Region Specific Score ---------------------------------------
# region score
all(row.names(data)  == names(oracle_coefficient)) # T
signature <- t(data) %*% oracle$Model.Coefficient
signature <- as.data.frame(signature)

# patient score
mys <- rep(0, length(mypat))
names(mys) <-  mypat
mys.max <- mys.min <- mys.avg <- mys

for(k in 1:length(mypat)) {
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  mys.max[k]  <- max(signature[se, ])			
  mys.min[k]  <- min(signature[se, ])			 
  mys.avg[k]  <- mean(signature[se, ])			
}

res.m2 <- data.frame(m2.s.avg=mys.avg, 
                  m2.s.max=mys.max, 
                  m2.s.min=mys.min)

all(rownames(res.m1) == rownames(res.m2)) #T
res <- cbind(res.m1,res.m2)
write.table(res, myoutf1, sep="\t", quote=F)
write.table(signature, myoutf2, sep="\t", quote=F)
# [1.2] ORACLE Survival ========================================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- paste0(outdir,"/ORACLE_scores_TRACERx.txt")
myinf2 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myFig1 <- paste0(outdir,"/ORACLE_scores_forest_TRACERx.pdf")
myFig2 <- paste0(outdir,"/ORACLE_scores_CIbar_TRACERx.pdf")

# a) load data -----------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info <- read.csv(myinf2, header=T, row.names=1)


comxx <- intersect(row.names(data), row.names(info))
info <- info[comxx,]
data <- data[comxx,]

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

mydata <- cbind(data, info)
mydata <- mydata[mydata[, "t.surv"]>0,]

# b) survival analysis ---------------------------------------------------------

p.data <- matrix(0, ncol(data), 5)
rownames(p.data) <- colnames(data)
colnames(p.data) <- c("HR","P-Value","C-index", "lower.95", "upper.95")
p.data <- as.data.frame(p.data)


for(k in 1:nrow(p.data)){
  cat("\r", k)
  # max
  formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
  mycox <- coxph(formula,mydata) 
  mycox <- summary(mycox)
  tmp <- mycox$conf.int
  p.data$HR[k] <- round(tmp[1],3)
  p.data$lower.95[k] <- tmp[3]
  p.data$upper.95[k] <- tmp[4]
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],3)
  p.data$`C-index`[k] <- round(mycox$concordance[1],3)
  
}

# c) forest plot ---------------------------------------------------------------
tabletext <- cbind(c("Method",c(rep("M1",5),rep("M2",3))),
                     c("Function","Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min."),
                     c("Hazard Ratio",p.data$HR),
                     c("P-value",p.data$`P-Value`),
                     c("C-index",p.data$`C-index`))
pdf(myFig1, width = 6,height = 4)
forestplot(
    mean=c(NA,p.data$HR), 
    lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
    labeltext=tabletext, graph.pos=3, 
    title="ORACLE",
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
                      "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                      "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                      "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()

# d) bar plot ---------------------------------------------------------------

p1 <-  p.data %>% 
  mutate( Method = c(rep("M1",5),rep("M2",3))) %>%
  mutate(label = factor(rownames(.),levels = rownames(.))) %>%
  arrange(desc(label)) %>%
  mutate(label = factor(label,levels = label)) %>%
  ggplot(mapping = aes(x = label,y=`C-index`-0.5,fill = Method)) +
  geom_hline(lty=2,size = 0.9,yintercept = 0.1,color="gray") +
  geom_hline(lty=2,size = 0.9,yintercept = 0.2,color="gray") +
  geom_bar(width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  scale_fill_manual(guide = "none",
                    values = c("#d73027" , "#c51b7d" )) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  coord_flip() +
  labs(x="Function",y="C-index",
       title=paste0("ORACLE"),
       base_size = 8, base_family = "serif",face="bold") +
  theme(aspect.ratio = 4/2,
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
pdf(myFig2, width = 3, height = 4)
print(p1)
dev.off()

# [1.3] WTGS Score ===========================================================
rm(list=ls())
library(survival)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- "~/ChaoCheng/Rprogram/BASE/TCGA-Lung-Prog-Profile_4base.txt"
myinf3 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myfun1 <- "~/ChaoCheng/Rprogram/BASE/base5.R"

myoutf1 <- paste0(outdir,"/WTGS_scores_TRACERx.txt")
myoutf2 <- paste0(outdir,"/WTGS_M2_Region_TRACERx.txt")

# a) expression data  -------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, check.names=F)

# remove low expression 
xx = apply(data>0,1,sum)
se = which(xx>=20)
data = data[se,]

# log2
data <- log2(data+1)

# b) signature -----------------------------------------------------------------
mywt <- read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

comxx = intersect(row.names(data), row.names(mywt))
data = data[comxx,]
mywt = mywt[comxx,]
mywt = mywt[, c(1:2, 5:6)]
dim(data) # 13469   161

all(row.names(data)  == rownames(mywt)) # T
# c) survival info -------------------------------------------------------------
info <- read.csv(myinf3, header=T, row.names=1)

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 17, 25)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23   161

# d) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
all(row.names(tmp)  == names(WTGS_coefficient)) # T
colnames(tmp) <- mypat
dat.rev <- dat.adj <- dat.avg <- dat.max <- dat.min <- tmp
# patient gene expression 
for(k in 1:length(mypat)){
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  tmp <- data[, se, drop=F]
  dat.avg[,k] <- apply(tmp,1,mean)
  xx1 <- apply(tmp,1, max)
  xx2 <- apply(tmp,1, min)
  dat.adj[,k] <- ifelse(mywt[,1] > 0, xx2, xx1)
  dat.rev[,k] <- ifelse(mywt[,1] > 0, xx1, xx2)
  dat.max[,k] <- xx1
  dat.min[,k] <- xx2
}

# patient score 
source(myfun1)
reg <- mywt
xx1 <- base5(dat.avg, reg, perm=1000, myoutf1, median.norm=T)
xx2 <- base5(dat.adj, reg, perm=1000, myoutf1, median.norm=T)
xx3 <- base5(dat.rev, reg, perm=1000, myoutf1, median.norm=T)
xx4 <- base5(dat.max, reg, perm=1000, myoutf1, median.norm=T)
xx5 <- base5(dat.min, reg, perm=1000, myoutf1, median.norm=T)
tmp <- xx1[[1]]
m1.s.avg <- tmp[,1] - tmp[,3]
tmp <- xx2[[1]]
m1.s.adj <- tmp[,1] - tmp[,3]
tmp <- xx3[[1]]
m1.s.rev <- tmp[,1] - tmp[,3]
tmp <- xx4[[1]]
m1.s.max <- tmp[,1] - tmp[,3]
tmp <- xx5[[1]]
m1.s.min <- tmp[,1] - tmp[,3]
res.m1 <- data.frame(m1.s.avg,  m1.s.max, m1.s.min, m1.s.adj, m1.s.rev)

row.names(res.m1) <- colnames(dat.avg)


# e) M2 Region Specific Score ---------------------------------------
# region score
all(row.names(data)  == rownames(mywt)) # T

reg <-mywt
xx <-base5(data, reg, perm=1000, myoutf1, median.norm=T)
tmp <-xx[[1]]
score <-tmp[,1] - tmp[,3]
signature <-data.frame(row.names =colnames(data), score)


# patient score
mys <- rep(0, length(mypat))
names(mys) <-  mypat
mys.max <- mys.min <- mys.avg <- mys

for(k in 1:length(mypat)) {
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  mys.max[k]  <- max(signature[se, ])			
  mys.min[k]  <- min(signature[se, ])			 
  mys.avg[k]  <- mean(signature[se, ])			
}

res.m2 <- data.frame(m2.s.avg=mys.avg, 
                     m2.s.max=mys.max, 
                     m2.s.min=mys.min)

all(rownames(res.m1) == rownames(res.m2)) #T
res <- cbind(res.m1,res.m2)
write.table(res, myoutf1, sep="\t", quote=F)
write.table(signature, myoutf2, sep="\t", quote=F)

# [1.4] WTGS Survival ========================================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- paste0(outdir,"/WTGS_scores_TRACERx.txt")
myinf2 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myFig1 <- paste0(outdir,"/WTGS_scores_forest_TRACERx.pdf")
myFig2 <- paste0(outdir,"/WTGS_scores_CIbar_TRACERx.pdf")

# a) load data -----------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info <- read.csv(myinf2, header=T, row.names=1)


comxx <- intersect(row.names(data), row.names(info))
info <- info[comxx,]
data <- data[comxx,]

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

mydata <- cbind(data, info)
mydata <- mydata[mydata[, "t.surv"]>0,]

# b) survival analysis ---------------------------------------------------------

p.data <- matrix(0, ncol(data), 5)
rownames(p.data) <- colnames(data)
colnames(p.data) <- c("HR","P-Value","C-index", "lower.95", "upper.95")
p.data <- as.data.frame(p.data)


for(k in 1:nrow(p.data)){
  cat("\r", k)
  # max
  formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
  mycox <- coxph(formula,mydata) 
  mycox <- summary(mycox)
  tmp <- mycox$conf.int
  p.data$HR[k] <- round(tmp[1],3)
  p.data$lower.95[k] <- tmp[3]
  p.data$upper.95[k] <- tmp[4]
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],3)
  p.data$`C-index`[k] <- round(mycox$concordance[1],3)
  
}

# c) forest plot ---------------------------------------------------------------
tabletext <- cbind(c("Method",c(rep("M1",5),rep("M2",3))),
                   c("Function","Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min."),
                   c("Hazard Ratio",p.data$HR),
                   c("P-value",p.data$`P-Value`),
                   c("C-index",p.data$`C-index`))
pdf(myFig1, width = 6,height = 4)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="WTGS",
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()

# d) bar plot ---------------------------------------------------------------

p1 <-  p.data %>% 
  mutate( Method = c(rep("M1",5),rep("M2",3))) %>%
  mutate(label = factor(rownames(.),levels = rownames(.))) %>%
  arrange(desc(label)) %>%
  mutate(label = factor(label,levels = label)) %>%
  ggplot(mapping = aes(x = label,y=`C-index`-0.5,fill = Method)) +
  geom_hline(lty=2,size = 0.9,yintercept = 0.1,color="gray") +
  geom_hline(lty=2,size = 0.9,yintercept = 0.2,color="gray") +
  geom_bar(width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  scale_fill_manual(guide = "none",
                    values = c("#d73027" , "#c51b7d" )) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  coord_flip() +
  labs(x="Function",y="C-index",
       title=paste0("WTGS"),
       base_size = 8, base_family = "serif",face="bold") +
  theme(aspect.ratio = 4/2,
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
pdf(myFig2, width = 3, height = 4)
print(p1)
dev.off()

# [1.5] PACEG Score ===========================================================
rm(list=ls())
library(survival)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- paste0(outdir,"/PACEG.txt")
myinf3 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myoutf1 <- paste0(outdir,"/PACEG_scores_TRACERx.txt")
myoutf2 <- paste0(outdir,"/PACEG_M2_Region_TRACERx.txt")
# a) expression and signature  -------------------------------------------------

data <- read.table(myinf1, sep="\t", header=T, check.names=F)
PACEG <- read.table(myinf2,sep="\t",header = T)

# confirm the same order of gene
row.names(PACEG) <- PACEG$gene
comxx <- intersect(row.names(data), row.names(PACEG))
data <- data[comxx,]
PACEG <- PACEG[comxx,]
all(rownames(data) == PACEG$gene) # T


PACEG_coefficient <- PACEG$coeff
names(PACEG_coefficient) <- PACEG$gene
# b) log transform  ------------------------------------------------------------

data <- log2(data+1)
# c) survival info -------------------------------------------------------------
info <- read.csv(myinf3, header=T, row.names=1)

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 17, 25)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23   161

# d) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
all(row.names(tmp)  == names(PACEG_coefficient)) # T
colnames(tmp) <- mypat
dat.rev <- dat.adj <- dat.avg <- dat.max <- dat.min <- tmp
# patient gene expression 
for(k in 1:length(mypat)){
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  tmp <- data[, se, drop=F]
  dat.max[,k] <- apply(tmp, 1, max)			
  dat.min[,k] <- apply(tmp, 1, min)			
  dat.avg[,k] <- apply(tmp, 1, mean)		
  dat.adj[,k] <- ifelse(test = PACEG_coefficient > 0, 
                        yes = dat.max[,k],  no = dat.min[,k])
  dat.rev[,k] <- ifelse(test = PACEG_coefficient < 0, 
                        yes = dat.max[,k],  no = dat.min[,k])
}

# patient score 
res.m1 <- data.frame(m1.s.avg=t(dat.avg) %*% PACEG_coefficient, 
                     m1.s.max=t(dat.max) %*% PACEG_coefficient, 
                     m1.s.min=t(dat.min) %*% PACEG_coefficient,
                     m1.s.adj=t(dat.adj) %*% PACEG_coefficient,
                     m1.s.rev=t(dat.rev) %*% PACEG_coefficient )

# e) M2 Region Specific Score ---------------------------------------
# region score
all(row.names(data)  == names(PACEG_coefficient)) # T
signature <- t(data) %*%  PACEG$coeff
signature <- as.data.frame(signature)

# patient score
mys <- rep(0, length(mypat))
names(mys) <-  mypat
mys.max <- mys.min <- mys.avg <- mys

for(k in 1:length(mypat)) {
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  mys.max[k]  <- max(signature[se, ])			
  mys.min[k]  <- min(signature[se, ])			 
  mys.avg[k]  <- mean(signature[se, ])			
}

res.m2 <- data.frame(m2.s.avg=mys.avg, 
                     m2.s.max=mys.max, 
                     m2.s.min=mys.min)

all(rownames(res.m1) == rownames(res.m2)) #T
res <- cbind(res.m1,res.m2)
write.table(res, myoutf1, sep="\t", quote=F)
write.table(signature, myoutf2, sep="\t", quote=F)
# [1.6] PACEG Survival ========================================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- paste0(outdir,"/PACEG_scores_TRACERx.txt")
myinf2 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

myFig1 <- paste0(outdir,"/PACEG_scores_forest_TRACERx.pdf")
myFig2 <- paste0(outdir,"/PACEG_scores_CIbar_TRACERx.pdf")

# a) load data -----------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info <- read.csv(myinf2, header=T, row.names=1)


comxx <- intersect(row.names(data), row.names(info))
info <- info[comxx,]
data <- data[comxx,]

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

mydata <- cbind(data, info)
mydata <- mydata[mydata[, "t.surv"]>0,]

# b) survival analysis ---------------------------------------------------------

p.data <- matrix(0, ncol(data), 5)
rownames(p.data) <- colnames(data)
colnames(p.data) <- c("HR","P-Value","C-index", "lower.95", "upper.95")
p.data <- as.data.frame(p.data)


for(k in 1:nrow(p.data)){
  cat("\r", k)
  # max
  formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
  mycox <- coxph(formula,mydata) 
  mycox <- summary(mycox)
  tmp <- mycox$conf.int
  p.data$HR[k] <- round(tmp[1],3)
  p.data$lower.95[k] <- tmp[3]
  p.data$upper.95[k] <- tmp[4]
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],3)
  p.data$`C-index`[k] <- round(mycox$concordance[1],3)
  
}

# c) forest plot ---------------------------------------------------------------
tabletext <- cbind(c("Method",c(rep("M1",5),rep("M2",3))),
                   c("Function","Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min."),
                   c("Hazard Ratio",p.data$HR),
                   c("P-value",p.data$`P-Value`),
                   c("C-index",p.data$`C-index`))
pdf(myFig1, width = 6,height = 4)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="PACEG",
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()

# d) bar plot ---------------------------------------------------------------

p1 <-  p.data %>% 
  mutate( Method = c(rep("M1",5),rep("M2",3))) %>%
  mutate(label = factor(rownames(.),levels = rownames(.))) %>%
  arrange(desc(label)) %>%
  mutate(label = factor(label,levels = label)) %>%
  ggplot(mapping = aes(x = label,y=`C-index`-0.5,fill = Method)) +
  geom_hline(lty=2,size = 0.9,yintercept = 0.1,color="gray") +
  geom_hline(lty=2,size = 0.9,yintercept = 0.2,color="gray") +
  geom_bar(width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  scale_fill_manual(guide = "none",
                    values = c("#d73027" , "#c51b7d" )) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  coord_flip() +
  labs(x="Function",y="C-index",
       title=paste0("PACEG"),
       base_size = 8, base_family = "serif",face="bold") +
  theme(aspect.ratio = 4/2,
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
pdf(myFig2, width = 3, height = 4)
print(p1)
dev.off()
# [2] CPRIT-MIRA ##################################################################

# [2.1] ORACLE Score ===========================================================
rm(list=ls())
library(survival)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- paste0(outdir,"/ORACLE.txt")
myinf3 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

myoutf1 <- paste0(outdir,"/ORACLE_scores_CPRIT-MIRA.txt")
myoutf2 <- paste0(outdir,"/ORACLE_M2_Region_CPRIT-MIRA.txt")
# a) expression and signature  -------------------------------------------------

data <- read.table(myinf1, sep="\t", header=T, check.names=F)
oracle <- read.table(myinf2,sep="\t",header = T)

# remove normal
se <- which(substr(colnames(data), 11, 11) =="T"  )
data <- data[,se]
dim(data)		

# confirm the same order of gene
row.names(oracle) <- oracle$Gene.Symbol
comxx <- intersect(row.names(data), row.names(oracle))
data <- data[comxx,]
oracle <- oracle[comxx,]
all(rownames(data) == oracle$Gene.Symbol) # T


oracle_coefficient <- oracle$Model.Coefficient
names(oracle_coefficient) <- oracle$Gene.Symbol
# b) log transform -------------------------------------------------------------

data <- log2(data+1)
# c) survival info -------------------------------------------------------------
info <- read.csv(myinf3, header=T, row.names=1)

e.surv <- replace(info$Relapse, 17, "Y") 
# 17 is continuous but other column show it has relapse
t.surv <- info$PFS..Date.of.the.biopsy.to.confirm.relapse.Date.of.surgery..months

# To be consistent with cohortD, the event is Recurrence.or.death
se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "Y")
info[se,]
e.surv[se] <- info$Deceased..Y.N.[se]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "N")
info[se,]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

e.surv <- ifelse(e.surv == "Y",1,0)
info <- cbind(t.surv, e.surv, info)

pid <-  substr(colnames(data), 1, 6)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23 64

# d) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
all(row.names(tmp)  == names(oracle_coefficient)) # T
colnames(tmp) <- mypat
dat.rev <- dat.adj <- dat.avg <- dat.max <- dat.min <- tmp
# patient gene expression 
for(k in 1:length(mypat)){
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  tmp <- data[, se, drop=F]
  dat.max[,k] <- apply(tmp, 1, max)			
  dat.min[,k] <- apply(tmp, 1, min)			
  dat.avg[,k] <- apply(tmp, 1, mean)		
  dat.adj[,k] <- ifelse(test = oracle_coefficient > 0, 
                        yes = dat.max[,k],  no = dat.min[,k])
  dat.rev[,k] <- ifelse(test = oracle_coefficient < 0, 
                        yes = dat.max[,k],  no = dat.min[,k])
}

# patient score 
res.m1 <- data.frame(m1.s.avg=t(dat.avg) %*% oracle_coefficient, 
                     m1.s.max=t(dat.max) %*% oracle_coefficient, 
                     m1.s.min=t(dat.min) %*% oracle_coefficient,
                     m1.s.adj=t(dat.adj) %*% oracle_coefficient,
                     m1.s.rev=t(dat.rev) %*% oracle_coefficient )

# e) M2 Region Specific Score ---------------------------------------
# region score
all(row.names(data)  == names(oracle_coefficient)) # T
signature <- t(data) %*% oracle$Model.Coefficient
signature <- as.data.frame(signature)

# patient score
mys <- rep(0, length(mypat))
names(mys) <-  mypat
mys.max <- mys.min <- mys.avg <- mys

for(k in 1:length(mypat)) {
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  mys.max[k]  <- max(signature[se, ])			
  mys.min[k]  <- min(signature[se, ])			 
  mys.avg[k]  <- mean(signature[se, ])			
}

res.m2 <- data.frame(m2.s.avg=mys.avg, 
                     m2.s.max=mys.max, 
                     m2.s.min=mys.min)

all(rownames(res.m1) == rownames(res.m2)) #T
res <- cbind(res.m1,res.m2)
write.table(res, myoutf1, sep="\t", quote=F)
write.table(signature, myoutf2, sep="\t", quote=F)
# [2.2] ORACLE Survival ========================================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- paste0(outdir,"/ORACLE_scores_CPRIT-MIRA.txt")
myinf2 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

myFig1 <- paste0(outdir,"/ORACLE_scores_forest_CPRIT-MIRA.pdf")
myFig2 <- paste0(outdir,"/ORACLE_scores_CIbar_CPRIT-MIRA.pdf")
# a) load data -----------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info <- read.csv(myinf2, header=T, row.names=1)


comxx <- intersect(row.names(data), row.names(info))
info <- info[comxx,]
data <- data[comxx,]

# survival info
e.surv <- replace(info$Relapse, 17, "Y") 
# 17 is continuous but other column show it has relapse
t.surv <- info$PFS..Date.of.the.biopsy.to.confirm.relapse.Date.of.surgery..months

# To be consistent with cohortD, the event is Recurrence.or.death
se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "Y")
info[se,]
e.surv[se] <- info$Deceased..Y.N.[se]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "N")
info[se,]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

e.surv <- ifelse(e.surv == "Y",1,0)
info <- cbind(t.surv, e.surv, info)

mydata <- cbind(data, info)
mydata <- mydata[mydata[, "t.surv"]>0,]

# b) survival analysis ---------------------------------------------------------

p.data <- matrix(0, ncol(data), 5)
rownames(p.data) <- colnames(data)
colnames(p.data) <- c("HR","P-Value","C-index", "lower.95", "upper.95")
p.data <- as.data.frame(p.data)


for(k in 1:nrow(p.data)){
  cat("\r", k)
  # max
  formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
  mycox <- coxph(formula,mydata) 
  mycox <- summary(mycox)
  tmp <- mycox$conf.int
  p.data$HR[k] <- round(tmp[1],3)
  p.data$lower.95[k] <- tmp[3]
  p.data$upper.95[k] <- tmp[4]
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],3)
  p.data$`C-index`[k] <- round(mycox$concordance[1],3)
  
}

# The p-value is not significant due to limited cohort size,
# So we forcused on C-index as external validation of our findings.

# HR(m1.s.min, m1.s.rev and m2.s.min) < 1, 
# so it is C-index is calculated opposite to the m1.s.avg and m2.s.avg
# that is, m1.s.avg and m2.s.avg(higher score, higher risk)
# that is, m1.s.min, m1.s.rev and m2.s.min(higher score, lower risk)
# To calculate the consistent C-index for HR < 1
p.data$`C-index` <- ifelse(p.data$HR < 1, 1-p.data$`C-index`, p.data$`C-index`)


# c) forest plot ---------------------------------------------------------------
tabletext <- cbind(c("Method",c(rep("M1",5),rep("M2",3))),
                   c("Function","Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min."),
                   c("Hazard Ratio",p.data$HR),
                   c("P-value",p.data$`P-Value`),
                   c("C-index",p.data$`C-index`))
pdf(myFig1, width = 6,height = 4)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="ORACLE",
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()
# d) bar plot ---------------------------------------------------------------
p1 <-  p.data %>% 
  mutate( Method = c(rep("M1",5),rep("M2",3)),
  Function = c("Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min.")) %>%
  mutate(label = factor(rownames(.),levels = rownames(.))) %>%
  ggplot(mapping = aes(x = label,y=`C-index`-0.4, fill = Method)) +
  geom_hline(lty=4,size = 0.9,yintercept = 0.1,color="gray") +
  geom_vline(lty=4,size = 0.9,xintercept = 5.5,color="black") +
  annotate("text", x = 3, y= 0.175, size = 7, 
           label =  "Method 1", family = "serif") +
  annotate("text", x = 7, y= 0.175, size = 7, 
           label =  "Method 2", family = "serif") +
  geom_bar(width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  scale_fill_manual(guide = "none",
                      values = c("#3288bd" , "#5e4fa2" )) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Function",y="C-index",
       title=paste0("ORACLE"),
       base_size = 8, base_family = "serif",face="bold") +
  theme(aspect.ratio = 1/1,
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
    axis.text.x = element_text(face="bold", color="black", size=15),
    axis.text.y = element_text(face="bold", color="black", size=15),
    axis.title.x = element_text(face="bold", color="black", size=18),
    axis.title.y = element_text(face="bold",color="black", size=18)) +
    scale_x_discrete(labels = c("Avg.","Max.","Min.","Adj.","Rev.",
                                "Avg.","Max.","Min.")) +
  scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.6, 0.05),
                     labels = seq(0.4, 1, 0.05), expand = c(0,0)) 
pdf(myFig2, width = 8, height = 8)
print(p1)
dev.off()
# [2.3] WTGS Score ===========================================================
rm(list=ls())
library(survival)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- "~/ChaoCheng/Rprogram/BASE/TCGA-Lung-Prog-Profile_4base.txt"
myinf3 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

myfun1 <- "~/ChaoCheng/Rprogram/BASE/base5.R"

myoutf1 <- paste0(outdir,"/WTGS_scores_CPRIT-MIRA.txt")
myoutf2 <- paste0(outdir,"/WTGS_M2_Region_CPRIT-MIRA.txt")

# a) expression data  -------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, check.names=F)

# remove normal
se <- which(substr(colnames(data), 11, 11) =="T"  )
data <- data[,se]
dim(data)		

# remove low expression 
xx = apply(data>0,1,sum)
se = which(xx>=20)
data = data[se,]

# log2
data <- log2(data+1)

# b) signature -----------------------------------------------------------------
mywt <- read.table(myinf2, sep="\t", header=T, row.names=1, quote="")

comxx = intersect(row.names(data), row.names(mywt))
data = data[comxx,]
mywt = mywt[comxx,]
mywt = mywt[, c(1:2, 5:6)]
dim(data) # 13469   161

all(row.names(data)  == rownames(mywt)) # T
# c) survival info -------------------------------------------------------------
info <- read.csv(myinf3, header=T, row.names=1)

e.surv <- replace(info$Relapse, 17, "Y") 
# 17 is continuous but other column show it has relapse
t.surv <- info$PFS..Date.of.the.biopsy.to.confirm.relapse.Date.of.surgery..months

# To be consistent with cohortD, the event is Recurrence.or.death
se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "Y")
info[se,]
e.surv[se] <- info$Deceased..Y.N.[se]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "N")
info[se,]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

e.surv <- ifelse(e.surv == "Y",1,0)
info <- cbind(t.surv, e.surv, info)

pid <-  substr(colnames(data), 1, 6)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23   161

# d) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
all(row.names(tmp)  == names(WTGS_coefficient)) # T
colnames(tmp) <- mypat
dat.rev <- dat.adj <- dat.avg <- dat.max <- dat.min <- tmp
# patient gene expression 
for(k in 1:length(mypat)){
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  tmp <- data[, se, drop=F]
  dat.avg[,k] = apply(tmp,1,mean)
  xx1 = apply(tmp,1, max)
  xx2 = apply(tmp,1, min)
  dat.adj[,k] = ifelse(mywt[,1] > 0, xx2, xx1)
  dat.rev[,k] = ifelse(mywt[,1] > 0, xx1, xx2)
  dat.max[,k] = xx1
  dat.min[,k] = xx2
}

# patient score 
source(myfun1)
reg = mywt
xx1 = base5(dat.avg, reg, perm=1000, myoutf1, median.norm=T)
xx2 = base5(dat.adj, reg, perm=1000, myoutf1, median.norm=T)
xx3 = base5(dat.rev, reg, perm=1000, myoutf1, median.norm=T)
xx4 = base5(dat.max, reg, perm=1000, myoutf1, median.norm=T)
xx5 = base5(dat.min, reg, perm=1000, myoutf1, median.norm=T)
tmp = xx1[[1]]
m1.s.avg = tmp[,1] - tmp[,3]
tmp = xx2[[1]]
m1.s.adj = tmp[,1] - tmp[,3]
tmp = xx3[[1]]
m1.s.rev = tmp[,1] - tmp[,3]
tmp = xx4[[1]]
m1.s.max = tmp[,1] - tmp[,3]
tmp = xx5[[1]]
m1.s.min = tmp[,1] - tmp[,3]
res.m1 = data.frame(m1.s.avg,  m1.s.max, m1.s.min, m1.s.adj, m1.s.rev)

row.names(res.m1) = colnames(dat.avg)


# e) M2 Region Specific Score ---------------------------------------
# region score
all(row.names(data)  == rownames(mywt)) # T

reg = mywt
xx = base5(data, reg, perm=1000, myoutf1, median.norm=T)
tmp = xx[[1]]
score = tmp[,1] - tmp[,3]
signature = data.frame(row.names =colnames(data), score)


# patient score
mys <- rep(0, length(mypat))
names(mys) <-  mypat
mys.max <- mys.min <- mys.avg <- mys

for(k in 1:length(mypat)) {
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  mys.max[k]  <- max(signature[se, ])			
  mys.min[k]  <- min(signature[se, ])			 
  mys.avg[k]  <- mean(signature[se, ])			
}

res.m2 <- data.frame(m2.s.avg=mys.avg, 
                     m2.s.max=mys.max, 
                     m2.s.min=mys.min)

all(rownames(res.m1) == rownames(res.m2)) #T
res <- cbind(res.m1,res.m2)
write.table(res, myoutf1, sep="\t", quote=F)
write.table(signature, myoutf2, sep="\t", quote=F)
# [2.4] WTGS Survival ========================================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- paste0(outdir,"/WTGS_scores_CPRIT-MIRA.txt")
myinf2 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

myFig1 <- paste0(outdir,"/WTGS_scores_forest_CPRIT-MIRA.pdf")
myFig2 <- paste0(outdir,"/WTGS_scores_CIbar_CPRIT-MIRA.pdf")

# a) load data -----------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info <- read.csv(myinf2, header=T, row.names=1)

comxx <- intersect(row.names(data), row.names(info))
info <- info[comxx,]
data <- data[comxx,]

# survival info
e.surv <- replace(info$Relapse, 17, "Y") 
# 17 is continuous but other column show it has relapse
t.surv <- info$PFS..Date.of.the.biopsy.to.confirm.relapse.Date.of.surgery..months

# To be consistent with cohortD, the event is Recurrence.or.death
se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "Y")
info[se,]
e.surv[se] <- info$Deceased..Y.N.[se]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "N")
info[se,]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

e.surv <- ifelse(e.surv == "Y",1,0)
info <- cbind(t.surv, e.surv, info)

mydata <- cbind(data, info)
mydata <- mydata[mydata[, "t.surv"]>0,]

# b) survival analysis ---------------------------------------------------------

p.data <- matrix(0, ncol(data), 5)
rownames(p.data) <- colnames(data)
colnames(p.data) <- c("HR","P-Value","C-index", "lower.95", "upper.95")
p.data <- as.data.frame(p.data)


for(k in 1:nrow(p.data)){
  cat("\r", k)
  # max
  formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
  mycox <- coxph(formula,mydata) 
  mycox <- summary(mycox)
  tmp <- mycox$conf.int
  p.data$HR[k] <- round(tmp[1],3)
  p.data$lower.95[k] <- tmp[3]
  p.data$upper.95[k] <- tmp[4]
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],3)
  p.data$`C-index`[k] <- round(mycox$concordance[1],3)
  
}

# The p-value is not significant due to limited cohort size,
# So we forcused on C-index as external validation of our findings.

# HR(m1.s.min, m1.s.rev and m2.s.min) > 1, 
# so it is C-index is calculated opposite to the m1.s.avg and m2.s.avg
# that is, m1.s.avg and m2.s.avg(higher score, lower risk)
# that is, m1.s.min, m1.s.rev and m2.s.min(higher score, higher risk)
# To calculate the consistent C-index for HR > 1
k=7 # HR = 1
formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
fit <- survreg(formula,mydata) 
concordance(fit)
p.data$`C-index` <- ifelse(p.data$HR >= 1, 1-p.data$`C-index`, p.data$`C-index`)



# c) forest plot ---------------------------------------------------------------
tabletext <- cbind(c("Method",c(rep("M1",5),rep("M2",3))),
                   c("Function","Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min."),
                   c("Hazard Ratio",p.data$HR),
                   c("P-value",p.data$`P-Value`),
                   c("C-index",p.data$`C-index`))
pdf(myFig1, width = 6,height = 4)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="WTGS",
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()

# d) bar plot ---------------------------------------------------------------
p1 <-  p.data %>% 
  mutate( Method = c(rep("M1",5),rep("M2",3)),
          Function = c("Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min.")) %>%
  mutate(label = factor(rownames(.),levels = rownames(.))) %>%
  ggplot(mapping = aes(x = label,y=`C-index`-0.4, fill = Method)) +
  geom_hline(lty=4,size = 0.9,yintercept = 0.1,color="gray") +
  geom_vline(lty=4,size = 0.9,xintercept = 5.5,color="black") +
  annotate("text", x = 3, y= 0.175, size = 7, 
           label =  "Method 1", family = "serif") +
  annotate("text", x = 7, y= 0.175, size = 7, 
           label =  "Method 2", family = "serif") +
  geom_bar(width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  scale_fill_manual(guide = "none",
                    values = c("#3288bd" , "#5e4fa2" )) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Function",y="C-index",
       title=paste0("WTGS"),
       base_size = 8, base_family = "serif",face="bold") +
  theme(aspect.ratio = 1/1,
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
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=18),
        axis.title.y = element_text(face="bold",color="black", size=18)) +
  scale_x_discrete(labels = c("Avg.","Max.","Min.","Adj.","Rev.",
                              "Avg.","Max.","Min.")) +
  scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.6, 0.05),
                     labels = seq(0.4, 1, 0.05), expand = c(0,0)) 
pdf(myFig2, width = 8, height = 8)
print(p1)
dev.off()
# [2.5] PACEG Score ===========================================================
rm(list=ls())
library(survival)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- paste0(outdir,"/PACEG.txt")
myinf3 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

myoutf1 <- paste0(outdir,"/PACEG_scores_CPRIT-MIRA.txt")
myoutf2 <- paste0(outdir,"/PACEG_M2_Region_CPRIT-MIRA.txt")
# a) expression and signature  -------------------------------------------------

data <- read.table(myinf1, sep="\t", header=T, check.names=F)
PACEG <- read.table(myinf2,sep="\t",header = T)

# remove normal
se <- which(substr(colnames(data), 11, 11) =="T"  )
data <- data[,se]
dim(data)	
# confirm the same order of gene
row.names(PACEG) <- PACEG$gene
comxx <- intersect(row.names(data), row.names(PACEG))
data <- data[comxx,]
PACEG <- PACEG[comxx,]
all(rownames(data) == PACEG$gene) # T


PACEG_coefficient <- PACEG$coeff
names(PACEG_coefficient) <- PACEG$gene
# b) log transform  ------------------------------------------------------------

data <- log2(data+1)
# c) survival info -------------------------------------------------------------
info <- read.csv(myinf3, header=T, row.names=1)

e.surv <- replace(info$Relapse, 17, "Y") 
# 17 is continuous but other column show it has relapse
t.surv <- info$PFS..Date.of.the.biopsy.to.confirm.relapse.Date.of.surgery..months

# To be consistent with cohortD, the event is Recurrence.or.death
se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "Y")
info[se,]
e.surv[se] <- info$Deceased..Y.N.[se]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "N")
info[se,]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

e.surv <- ifelse(e.surv == "Y",1,0)
info <- cbind(t.surv, e.surv, info)

pid <-  substr(colnames(data), 1, 6)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23   161

# d) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
all(row.names(tmp)  == names(PACEG_coefficient)) # T
colnames(tmp) <- mypat
dat.rev <- dat.adj <- dat.avg <- dat.max <- dat.min <- tmp
# patient gene expression 
for(k in 1:length(mypat)){
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  tmp <- data[, se, drop=F]
  dat.max[,k] <- apply(tmp, 1, max)			
  dat.min[,k] <- apply(tmp, 1, min)			
  dat.avg[,k] <- apply(tmp, 1, mean)		
  dat.adj[,k] <- ifelse(test = PACEG_coefficient > 0, 
                        yes = dat.max[,k],  no = dat.min[,k])
  dat.rev[,k] <- ifelse(test = PACEG_coefficient < 0, 
                        yes = dat.max[,k],  no = dat.min[,k])
}

# patient score 
res.m1 <- data.frame(m1.s.avg=t(dat.avg) %*% PACEG_coefficient, 
                     m1.s.max=t(dat.max) %*% PACEG_coefficient, 
                     m1.s.min=t(dat.min) %*% PACEG_coefficient,
                     m1.s.adj=t(dat.adj) %*% PACEG_coefficient,
                     m1.s.rev=t(dat.rev) %*% PACEG_coefficient )

# e) M2 Region Specific Score ---------------------------------------
# region score
all(row.names(data)  == names(PACEG_coefficient)) # T
signature <- t(data) %*%  PACEG$coeff
signature <- as.data.frame(signature)

# patient score
mys <- rep(0, length(mypat))
names(mys) <-  mypat
mys.max <- mys.min <- mys.avg <- mys

for(k in 1:length(mypat)) {
  cat("\r", k)
  se <- which(pid==mypat[k])
  if(length(se)==0){
    mys[k,] <- NA
    next
  }
  mys.max[k]  <- max(signature[se, ])			
  mys.min[k]  <- min(signature[se, ])			 
  mys.avg[k]  <- mean(signature[se, ])			
}

res.m2 <- data.frame(m2.s.avg=mys.avg, 
                     m2.s.max=mys.max, 
                     m2.s.min=mys.min)

all(rownames(res.m1) == rownames(res.m2)) #T
res <- cbind(res.m1,res.m2)
write.table(res, myoutf1, sep="\t", quote=F)
write.table(signature, myoutf2, sep="\t", quote=F)
# [2.6] PACEG Survival ========================================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"

myinf1 <- paste0(outdir,"/PACEG_scores_CPRIT-MIRA.txt")
myinf2 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

myFig1 <- paste0(outdir,"/PACEG_scores_forest_CPRIT-MIRA.pdf")
myFig2 <- paste0(outdir,"/PACEG_scores_CIbar_CPRIT-MIRA.pdf")

# a) load data -----------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
info <- read.csv(myinf2, header=T, row.names=1)


comxx <- intersect(row.names(data), row.names(info))
info <- info[comxx,]
data <- data[comxx,]

# survival info
e.surv <- replace(info$Relapse, 17, "Y") 
# 17 is continuous but other column show it has relapse
t.surv <- info$PFS..Date.of.the.biopsy.to.confirm.relapse.Date.of.surgery..months

# To be consistent with cohortD, the event is Recurrence.or.death
se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "Y")
info[se,]
e.surv[se] <- info$Deceased..Y.N.[se]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

se <- which(info$Relapse == "N" & info$Deceased..Y.N. == "N")
info[se,]
t.surv[se] <- info$OS..Date.of.Death.Date.of.Biopsy..months[se]

e.surv <- ifelse(e.surv == "Y",1,0)
info <- cbind(t.surv, e.surv, info)

mydata <- cbind(data, info)
mydata <- mydata[mydata[, "t.surv"]>0,]

# b) survival analysis ---------------------------------------------------------

p.data <- matrix(0, ncol(data), 5)
rownames(p.data) <- colnames(data)
colnames(p.data) <- c("HR","P-Value","C-index", "lower.95", "upper.95")
p.data <- as.data.frame(p.data)


for(k in 1:nrow(p.data)){
  cat("\r", k)
  # max
  formula <- as.formula(paste0( "Surv(t.surv, e.surv) ~ ", rownames(p.data)[k]))
  mycox <- coxph(formula,mydata) 
  mycox <- summary(mycox)
  tmp <- mycox$conf.int
  p.data$HR[k] <- round(tmp[1],3)
  p.data$lower.95[k] <- tmp[3]
  p.data$upper.95[k] <- tmp[4]
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],3)
  p.data$`C-index`[k] <- round(mycox$concordance[1],3)
  
}

# The p-value is not significant due to limited cohort size,
# So we forcused on C-index as external validation of our findings.

# HR(m1.s.min, m1.s.rev and m2.s.min) < 1, 
# so it is C-index is calculated opposite to the m1.s.avg and m2.s.avg
# that is, m1.s.avg and m2.s.avg(higher score, higher risk)
# that is, m1.s.min, m1.s.rev and m2.s.min(higher score, lower risk)
# To calculate the consistent C-index for HR < 1
p.data$`C-index` <- ifelse(p.data$HR < 1, 1-p.data$`C-index`, p.data$`C-index`)

# c) forest plot ---------------------------------------------------------------
tabletext <- cbind(c("Method",c(rep("M1",5),rep("M2",3))),
                   c("Function","Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min."),
                   c("Hazard Ratio",p.data$HR),
                   c("P-value",p.data$`P-Value`),
                   c("C-index",p.data$`C-index`))
pdf(myFig1, width = 6,height = 4)
forestplot(
  mean=c(NA,p.data$HR), 
  lower=c(NA,p.data$lower.95), upper=c(NA,p.data$upper.95),
  labeltext=tabletext, graph.pos=3, 
  title="PACEG",
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()

# d) bar plot ---------------------------------------------------------------
p1 <-  p.data %>% 
  mutate( Method = c(rep("M1",5),rep("M2",3)),
          Function = c("Avg.","Max.","Min.","Adj.","Rev.","Avg.","Max.","Min.")) %>%
  mutate(label = factor(rownames(.),levels = rownames(.))) %>%
  ggplot(mapping = aes(x = label,y=`C-index`-0.4, fill = Method)) +
  geom_hline(lty=4,size = 0.9,yintercept = 0.1,color="gray") +
  geom_vline(lty=4,size = 0.9,xintercept = 5.5,color="black") +
  annotate("text", x = 3, y= 0.2, size = 7, 
           label =  "Method 1", family = "serif") +
  annotate("text", x = 7, y= 0.2, size = 7, 
           label =  "Method 2", family = "serif") +
  geom_bar(width=0.8, position = position_dodge(width=0.01),
           stat="identity",alpha=0.5) +
  scale_fill_manual(guide = "none",
                    values = c("#3288bd" , "#5e4fa2" )) +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Function",y="C-index",
       title=paste0("PACEG"),
       base_size = 8, base_family = "serif",face="bold") +
  theme(aspect.ratio = 1/1,
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
        axis.text.x = element_text(face="bold", color="black", size=15),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=18),
        axis.title.y = element_text(face="bold",color="black", size=18)) +
  scale_x_discrete(labels = c("Avg.","Max.","Min.","Adj.","Rev.",
                              "Avg.","Max.","Min.")) +
  scale_y_continuous(limits=c(0, 0.22), breaks=seq(0, 0.6, 0.05),
                     labels = seq(0.4, 1, 0.05), expand = c(0,0)) 
pdf(myFig2, width = 8, height = 8)
print(p1)
dev.off()
