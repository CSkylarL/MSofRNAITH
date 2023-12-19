#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Supp. Fig 2
# Author: Chenyang Li
# 09/20/2023
# TRACERx 421
# LUAD subtype: TRACERx421_LUAD
# Additional validation: TRACERx 421 excluding TRACERx 100 (TRACERx421-100)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx421_LUAD #########################################################

# [1.1] ORACLE Score ===========================================================
rm(list=ls())
library(survival)
library(dplyr)
library(fst)
library(forestplot)
library(ggplot2)

outdir <- "./Supp_Fig2/"
dir.create(outdir)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"
myinf5 <- "./Fig3/ORACLE.txt"

myoutf1 <- paste0(outdir,"/ORACLE_scores_TRACERx421_LUAD.txt")
myoutf2 <- paste0(outdir,"/ORACLE_M2_Region_TRACERx421_LUAD.txt")

myFig1 <- paste0(outdir,"/ORACLE_scores_forest_TRACERx421_LUAD.pdf")
# a) expression and clinical  --------------------------------------------------
patient <- readRDS(myinf1)
tumor <- readRDS(myinf2)
region <- read_fst(myinf3) %>% 
  mutate(tumor_id_new = gsub(pattern = "Cluster",
                             replacement = "Tumour",
                             tumour_id_muttable_cruk))
rna <- read_fst(myinf4) 
rownames(rna) <- rna$gene_id
rna <- rna[-1]

table(patient$histology_multi_full_genomically.confirmed)
# LUAD  LUAD&LUSC LUAD&Other     LUADx2     LUADx3       LUSC      Other 
# 231          4          1          5          1        134         45 

# LUAD + all have LUAD subtypes: LUAD&LUSC LUAD&Other     LUADx2     LUADx3  

se <- grep(pattern = "LUAD", 
           x = patient$histology_multi_full_genomically.confirmed)
info <- patient[se,]
dim(info) # 242
table(info$histology_multi_full_genomically.confirmed)
# LUAD  LUAD&LUSC LUAD&Other     LUADx2     LUADx3 
# 231          4          1          5          1 
PID <- info$cruk_id
rownames(info) <- PID
tumor_LUAD <- tumor %>% filter (cruk_id %in% PID)
region_LUAD <- region %>% filter (cruk_id %in% PID)

se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
data <- rna[se]

# b) signature -----------------------------------------------------------------
oracle <- read.table(myinf5,sep="\t",header = T)

# confirm the same order of gene
row.names(oracle) <- oracle$Gene.Symbol
comxx <- intersect(row.names(data), row.names(oracle))
data <- data[comxx,]
oracle <- oracle[comxx,]
all(rownames(data) == oracle$Gene.Symbol) # T


oracle_coefficient <- oracle$Model.Coefficient
names(oracle_coefficient) <- oracle$Gene.Symbol
# c) log transform -------------------------------------------------------------

data <- log2(data+1)
# d) survival info -------------------------------------------------------------

e.surv <- info$cens_dfs_any_event 
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
dim(info) # 187 54

se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23 472


# e) M1 Trasformed gene expression score ---------------------------------------
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

# f) M2 Region Specific Score ---------------------------------------
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
# a) load data -----------------------------------------------------------------
data <- res
all(rownames(data) == rownames(info))

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
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],5)
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


# [1.3] WTGS Score ===========================================================
rm(list=ls())
library(survival)
library(dplyr)
library(fst)
library(forestplot)
library(ggplot2)

outdir <- "./Supp_Fig2/"
dir.create(outdir)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"
myinf5 <- "./BASE/TCGA-Lung-Prog-Profile_4base.txt"

myfun1 <- "./BASE/base5.R"

myoutf1 <- paste0(outdir,"/WTGS_scores_TRACERx421_LUAD.txt")
myoutf2 <- paste0(outdir,"/WTGS_M2_Region_TRACERx421_LUAD.txt")

myFig1 <- paste0(outdir,"/WTGS_scores_forest_TRACERx421_LUAD.pdf")
# a) expression and clinical  --------------------------------------------------
patient <- readRDS(myinf1)
tumor <- readRDS(myinf2)
region <- read_fst(myinf3) %>% 
  mutate(tumor_id_new = gsub(pattern = "Cluster",
                             replacement = "Tumour",
                             tumour_id_muttable_cruk))
rna <- read_fst(myinf4) 
rownames(rna) <- rna$gene_id
rna <- rna[-1]

table(patient$histology_multi_full_genomically.confirmed)
# LUAD  LUAD&LUSC LUAD&Other     LUADx2     LUADx3       LUSC      Other 
# 231          4          1          5          1        134         45 

# LUAD + all have LUAD subtypes: LUAD&LUSC LUAD&Other     LUADx2     LUADx3  

se <- grep(pattern = "LUAD", 
           x = patient$histology_multi_full_genomically.confirmed)
info <- patient[se,]
dim(info) # 242
table(info$histology_multi_full_genomically.confirmed)
# LUAD  LUAD&LUSC LUAD&Other     LUADx2     LUADx3 
# 231          4          1          5          1 
PID <- info$cruk_id
rownames(info) <- PID
tumor_LUAD <- tumor %>% filter (cruk_id %in% PID)
region_LUAD <- region %>% filter (cruk_id %in% PID)

se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
data <- rna[se]

# b) preprocess ----------------------------------------------------------------
# remove low expression 
xx <- apply(data>0,1,sum)
se <- which(xx>=20)
data <- data[se,]

# log2
data <- log2(data+1)


# b) signature -----------------------------------------------------------------
mywt <- read.table(myinf5, sep="\t", header=T, row.names=1, quote="")

comxx <- intersect(row.names(data), row.names(mywt))
data <- data[comxx,]
mywt <- mywt[comxx,]
mywt <- mywt[, c(1:2, 5:6)]
dim(data) # 12845   472

all(row.names(data)  == rownames(mywt)) # T

# d) survival info -------------------------------------------------------------

e.surv <- info$cens_dfs_any_event
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #   12845   472


# e) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
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


# f) M2 Region Specific Score ---------------------------------------
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
# a) load data -----------------------------------------------------------------
data <- res
all(rownames(data) == rownames(info))

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
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],7)
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()


# [1.5] PACEG Score ===========================================================
rm(list=ls())
library(survival)
library(dplyr)
library(fst)
library(forestplot)
library(ggplot2)

outdir <- "./Supp_Fig2/"
dir.create(outdir)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"
myinf5 <- "./Fig3/PACEG.txt"

myoutf1 <- paste0(outdir,"/PACEG_scores_TRACERx421_LUAD.txt")
myoutf2 <- paste0(outdir,"/PACEG_M2_Region_TRACERx421_LUAD.txt")

myFig1 <- paste0(outdir,"/PACEG_scores_forest_TRACERx421_LUAD.pdf")
# a) expression and clinical  --------------------------------------------------
patient <- readRDS(myinf1)
tumor <- readRDS(myinf2)
region <- read_fst(myinf3) %>% 
  mutate(tumor_id_new = gsub(pattern = "Cluster",
                             replacement = "Tumour",
                             tumour_id_muttable_cruk))
rna <- read_fst(myinf4) 
rownames(rna) <- rna$gene_id
rna <- rna[-1]

table(patient$histology_multi_full_genomically.confirmed)
# LUAD  LUAD&LUSC LUAD&Other     LUADx2     LUADx3       LUSC      Other 
# 231          4          1          5          1        134         45 

# LUAD + all have LUAD subtypes: LUAD&LUSC LUAD&Other     LUADx2     LUADx3  

se <- grep(pattern = "LUAD", 
           x = patient$histology_multi_full_genomically.confirmed)
info <- patient[se,]
dim(info) # 242
table(info$histology_multi_full_genomically.confirmed)
# LUAD  LUAD&LUSC LUAD&Other     LUADx2     LUADx3 
# 231          4          1          5          1 
PID <- info$cruk_id
rownames(info) <- PID
tumor_LUAD <- tumor %>% filter (cruk_id %in% PID)
region_LUAD <- region %>% filter (cruk_id %in% PID)

se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
data <- rna[se]

# b) signature -----------------------------------------------------------------
PACEG <- read.table(myinf5,sep="\t",header = T)

# confirm the same order of gene
row.names(PACEG) <- PACEG$gene
comxx <- intersect(row.names(data), row.names(PACEG))

# check the name problem
PACEG[which(! PACEG$gene %in% comxx),]
# gene     coeff
# C1orf172 C1orf172 -6.51e-02 --> KDF1
# C9orf69   C9orf69  2.60e-07 --> TMEM250
"KDF1" %in% rownames(data) # T
"TMEM250" %in% rownames(data) # T

PACEG["C1orf172",1] <- "KDF1"
PACEG["C9orf69",1] <- "TMEM250"

# Redo: confirm the same order of gene
row.names(PACEG) <- PACEG$gene
comxx <- intersect(row.names(data), row.names(PACEG))


data <- data[comxx,]
PACEG <- PACEG[comxx,]
all(rownames(data) == PACEG$gene) # T


PACEG_coefficient <- PACEG$coeff
names(PACEG_coefficient) <- PACEG$gene
# c) log transform -------------------------------------------------------------

data <- log2(data+1)
# d) survival info -------------------------------------------------------------

e.surv <- info$cens_dfs_any_event
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #   26 472


# e) M1 Trasformed gene expression score ---------------------------------------
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

# f) M2 Region Specific Score ---------------------------------------
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
# a) load data -----------------------------------------------------------------
data <- res
all(rownames(data) == rownames(info))

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
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],6)
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()



# [2] TRACERx421_nonTRACERx100 #########################################################

# [2.1] ORACLE Score ===========================================================
rm(list=ls())
library(survival)
library(dplyr)
library(fst)
library(forestplot)
library(ggplot2)

outdir <- "./Supp_Fig2/"
dir.create(outdir)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"
myinf5 <- "./Fig3/ORACLE.txt"

myoutf1 <- paste0(outdir,"/ORACLE_scores_TRACERx421_nonTRACERx100.txt")
myoutf2 <- paste0(outdir,"/ORACLE_M2_Region_TRACERx421_nonTRACERx100.txt")

myFig1 <- paste0(outdir,"/ORACLE_scores_forest_TRACERx421_nonTRACERx100.pdf")
# a) expression and clinical  --------------------------------------------------
patient <- readRDS(myinf1)
tumor <- readRDS(myinf2)
region <- read_fst(myinf3) %>% 
  mutate(tumor_id_new = gsub(pattern = "Cluster",
                             replacement = "Tumour",
                             tumour_id_muttable_cruk))
rna <- read_fst(myinf4) 
rownames(rna) <- rna$gene_id
rna <- rna[-1]

table(patient$tx100)
# FALSE  TRUE 
# 323    98 

se <- which(patient$tx100 == FALSE)
info <- patient[se,]
dim(info) # 323  50

PID <- info$cruk_id
rownames(info) <- PID
tumor_LUAD <- tumor %>% filter (cruk_id %in% PID)
region_LUAD <- region %>% filter (cruk_id %in% PID)

se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
data <- rna[se]

# b) signature -----------------------------------------------------------------
oracle <- read.table(myinf5,sep="\t",header = T)

# confirm the same order of gene
row.names(oracle) <- oracle$Gene.Symbol
comxx <- intersect(row.names(data), row.names(oracle))
data <- data[comxx,]
oracle <- oracle[comxx,]
all(rownames(data) == oracle$Gene.Symbol) # T


oracle_coefficient <- oracle$Model.Coefficient
names(oracle_coefficient) <- oracle$Gene.Symbol
# c) log transform -------------------------------------------------------------

data <- log2(data+1)
# d) survival info -------------------------------------------------------------

e.surv <- info$cens_dfs_any_event
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #  23 472


# e) M1 Trasformed gene expression score ---------------------------------------
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

# f) M2 Region Specific Score ---------------------------------------
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
# a) load data -----------------------------------------------------------------
data <- res
all(rownames(data) == rownames(info))

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
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],5)
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


# [2.3] WTGS Score ===========================================================
rm(list=ls())
library(survival)
library(dplyr)
library(fst)
library(forestplot)
library(ggplot2)

outdir <- "./Supp_Fig2/"
dir.create(outdir)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"
myinf5 <- "./BASE/TCGA-Lung-Prog-Profile_4base.txt"

myfun1 <- "./BASE/base5.R"

myoutf1 <- paste0(outdir,"/WTGS_scores_TRACERx421_nonTRACERx100.txt")
myoutf2 <- paste0(outdir,"/WTGS_M2_Region_TRACERx421_nonTRACERx100.txt")

myFig1 <- paste0(outdir,"/WTGS_scores_forest_TRACERx421_nonTRACERx100.pdf")
# a) expression and clinical  --------------------------------------------------
patient <- readRDS(myinf1)
tumor <- readRDS(myinf2)
region <- read_fst(myinf3) %>% 
  mutate(tumor_id_new = gsub(pattern = "Cluster",
                             replacement = "Tumour",
                             tumour_id_muttable_cruk))
rna <- read_fst(myinf4) 
rownames(rna) <- rna$gene_id
rna <- rna[-1]

table(patient$tx100)
# FALSE  TRUE 
# 323    98 

se <- which(patient$tx100 == FALSE)
info <- patient[se,]
dim(info) # 323  50

PID <- info$cruk_id
rownames(info) <- PID
tumor_LUAD <- tumor %>% filter (cruk_id %in% PID)
region_LUAD <- region %>% filter (cruk_id %in% PID)

se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
data <- rna[se]


# b) preprocess ----------------------------------------------------------------
# remove low expression 
xx <- apply(data>0,1,sum)
se <- which(xx>=20)
data <- data[se,]

# log2
data <- log2(data+1)


# b) signature -----------------------------------------------------------------
mywt <- read.table(myinf5, sep="\t", header=T, row.names=1, quote="")

comxx <- intersect(row.names(data), row.names(mywt))
data <- data[comxx,]
mywt <- mywt[comxx,]
mywt <- mywt[, c(1:2, 5:6)]
dim(data) # 12845   472

all(row.names(data)  == rownames(mywt)) # T

# d) survival info -------------------------------------------------------------

e.surv <- info$cens_dfs_any_event
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #   12845   472


# e) M1 Trasformed gene expression score ---------------------------------------
tmp <- matrix(0, nrow(data), length(mypat))
row.names(tmp) <- row.names(data)
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


# f) M2 Region Specific Score ---------------------------------------
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
# [2.4] WTGS Survival ========================================================
# a) load data -----------------------------------------------------------------
data <- res
all(rownames(data) == rownames(info))

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
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],7)
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()


# [2.5] PACEG Score ===========================================================
rm(list=ls())
library(survival)
library(dplyr)
library(fst)
library(forestplot)
library(ggplot2)

outdir <- "./Supp_Fig2/"
dir.create(outdir)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"
myinf5 <- "./Fig3/PACEG.txt"

myoutf1 <- paste0(outdir,"/PACEG_scores_TRACERx421_nonTRACERx100.txt")
myoutf2 <- paste0(outdir,"/PACEG_M2_Region_TRACERx421_nonTRACERx100.txt")

myFig1 <- paste0(outdir,"/PACEG_scores_forest_TRACERx421_nonTRACERx100.pdf")
# a) expression and clinical  --------------------------------------------------
patient <- readRDS(myinf1)
tumor <- readRDS(myinf2)
region <- read_fst(myinf3) %>% 
  mutate(tumor_id_new = gsub(pattern = "Cluster",
                             replacement = "Tumour",
                             tumour_id_muttable_cruk))
rna <- read_fst(myinf4) 
rownames(rna) <- rna$gene_id
rna <- rna[-1]

table(patient$tx100)
# FALSE  TRUE 
# 323    98 

se <- which(patient$tx100 == FALSE)
info <- patient[se,]
dim(info) # 323  50

PID <- info$cruk_id
rownames(info) <- PID
tumor_LUAD <- tumor %>% filter (cruk_id %in% PID)
region_LUAD <- region %>% filter (cruk_id %in% PID)

se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
data <- rna[se]
# b) signature -----------------------------------------------------------------
PACEG <- read.table(myinf5,sep="\t",header = T)

# confirm the same order of gene
row.names(PACEG) <- PACEG$gene
comxx <- intersect(row.names(data), row.names(PACEG))

# check the name problem
PACEG[which(! PACEG$gene %in% comxx),]
# gene     coeff
# C1orf172 C1orf172 -6.51e-02 --> KDF1
# C9orf69   C9orf69  2.60e-07 --> TMEM250
"KDF1" %in% rownames(data) # T
"TMEM250" %in% rownames(data) # T

PACEG["C1orf172",1] <- "KDF1"
PACEG["C9orf69",1] <- "TMEM250"

# Redo: confirm the same order of gene
row.names(PACEG) <- PACEG$gene
comxx <- intersect(row.names(data), row.names(PACEG))


data <- data[comxx,]
PACEG <- PACEG[comxx,]
all(rownames(data) == PACEG$gene) # T


PACEG_coefficient <- PACEG$coeff
names(PACEG_coefficient) <- PACEG$gene
# c) log transform -------------------------------------------------------------

data <- log2(data+1)
# d) survival info -------------------------------------------------------------

e.surv <- info$cens_dfs_any_event
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(data), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
data <- data[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(data) #   26 472


# e) M1 Trasformed gene expression score ---------------------------------------
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

# f) M2 Region Specific Score ---------------------------------------
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
# a) load data -----------------------------------------------------------------
data <- res
all(rownames(data) == rownames(info))

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
  p.data$`P-Value`[k] <- round(mycox$coefficients[5],5)
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
                    "7" = gpar(lty = 2, lwd = 2, columns = 1:6, col = "#737373"),
                    "2" = gpar(lwd = 1,lwd = 2,  columns = 1:6, col = "#000044"),
                    "10" = gpar(lwd = 1,lwd = 3,  columns = 1:6, col = "#000044")))
dev.off()