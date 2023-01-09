#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.5
# Author: Chao Cheng; Chenyang Li
# TRACERx Color: "#67000d",[ "#d73027" ], "#c51b7d"
# CPRIT-MIRA Color: "#08306b", [ "#3288bd" ], "#5e4fa2" 
# CPRIT-MIRA  Result is not significant --> remove it
# 12/14/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx ##################################################################
# [1.1] Integrate score & Predict Region Specific Risk =========================
rm(list=ls())

library(survival)


indir1 <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"
myinf1 <- paste0(indir1,"/ORACLE_M2_Region_TRACERx.txt")
myinf2 <- paste0(indir1,"/WTGS_M2_Region_TRACERx.txt")
myinf3 <- paste0(indir1,"/PACEG_M2_Region_TRACERx.txt")

indir2 <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig4/"
myinf4 <- paste0(indir2,"/RNAseq_Genes_representativeImm_TRACERx.txt")

myinf5 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/RNAseq_Sample_FullID.txt"
myinf6 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig5/"

myoutf1 <- paste0(outdir,"/Region_specific_scores_TRACERx.txt")

# a) ORACLE --------------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1)
oracle <- data
# c) WTGS ----------------------------------------------------------------------
data <- read.table(myinf2, sep="\t", header=T, row.names=1)
wtgs <- data
# b) PACEG ---------------------------------------------------------------------
data <- read.table(myinf3, sep="\t", header=T, row.names=1)
paceg <- data

# d) immune cell ---------------------------------------------------------------
data <- read.table(myinf4, sep="\t", header=T, row.names=1)
cnum <- ncol(data)/4
data <- data[, 1:cnum] - data[, (cnum+1):(2*cnum)]
tmp <- colnames(data)
tmp <- gsub("_up\\.ES", "", tmp)
colnames(data) <- tmp
# since NK is not significant -- remove it
xx <- data[, -5]
imm <- xx

# e) integrate -----------------------------------------------------------------
comxx <- intersect(row.names(imm), row.names(oracle))
comxx <- intersect(row.names(wtgs), comxx)
comxx <- intersect(row.names(paceg), comxx)

oracle <- oracle[comxx,1]
paceg <- paceg[comxx,1]
wtgs <- wtgs[comxx,1]
imm <- imm[comxx,]

data <- cbind(imm, oracle, paceg, wtgs)

# f) map ID --------------------------------------------------------------------
info <- read.table(myinf5, sep="\t", header=T, 
                  stringsAsFactors=F, quote="", comment.char="")
map <- info[,1]
names(map) <- info[,2]

xx <- substr(row.names(data), 1, 15)
row.names(data) <- map[xx]

# g) Use max/min to construct model at the patient level -----------------------
colnames(data)
# [1] "NavB"  "MemB"  "CD8T"     "CD4T"  "Monocyte" "oracle"   "paceg"   "wtgs"    
flag <- c(-1,-1,-1,-1, 1, 1, 1, -1)

pid <- substr(row.names(data), 1, 8)
mypat <- unique(pid)

mydat <- matrix(0, length(mypat), ncol(data))
row.names(mydat) <- mypat
colnames(mydat) <- colnames(data)
mydat <- as.data.frame(mydat)

for(k in 1:length(mypat)) {
  se <- which(pid==mypat[k])
  tmp <- data[se, , drop=F]
  for(i in 1:ncol(data)) {
    if(flag[i]==-1) {
      mydat[k, i] <- min(tmp[,i])
    }
    if(flag[i]==1) {
      mydat[k, i] <- max(tmp[,i])
    }
    if(flag[i]==0) {
      mydat[k, i] <- mean(tmp[,i])
    } 
  }
}


# h) survival info -------------------------------------------------------------
info <- read.csv(myinf6, header=T, row.names=1)

comxx <- intersect(row.names(mydat), row.names(info))
info <- info[comxx,]
mydat <- mydat[comxx,]

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)



# i) predict risk score for each regions ---------------------------------------
score <- data
all(colnames(score) == colnames(mydat)) # T
for(k in 1:ncol(data)) {
  mytf <- mydat[,k]
  xx <- cbind(info, mytf)
  mycox <- coxph(Surv(t.surv, e.surv)~mytf, xx)
  tmp.dat <- data[,k,drop=F]
  colnames(tmp.dat)="mytf"
  xx <- predict(mycox, newdata= tmp.dat, type = "risk")
  score[,k] <- xx
}
score[1:3,]

# j) export --------------------------------------------------------------------
write.table(score, myoutf1, sep="\t", quote=F)

# [1.2] Correlation of Region Specific Risk ====================================
rm(list=ls())

library(ComplexHeatmap)
library(circlize)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig5/"

myinf1 <- paste0(outdir,"/Region_specific_scores_TRACERx.txt")
myFig1 <- paste0(outdir,"/Region_risk_scores_TRACERx.pdf")
myFig2 <- paste0(outdir,"/Correlation_TRACERx.pdf")



# a) choose samples with multiple regions (>=5) --------------------------------
data <- read.table(myinf1, sep="\t", header=T, stringsAsFactors=F, 
                  quote="", comment.char="")


# /Region_risk_scores_TRACERx_median-11.pdf: Normalize to [-1,1]
# data <- apply(data, 2, function(x) (x-median(x))/(max(x)-min(x)))
# Region_risk_scores_TRACERx_scale.pdf: [column-mean]/sd
# data <- scale(data)
# select regions >=3 use the below codes
# Heatmap(as.matrix(data),  cluster_rows = F, 
# heatmap_legend_param = list(title = "Risk Score"), 
# cluster_columns = F)

tag <- substr(row.names(data), 1, 8)
xx <- table(tag)
mul.sam  <- names(xx)[xx>=5]

se <- which(tag %in% mul.sam)

data <- data[se, ]
data <- round(data, 4)
colnames(data) <- c("Naive B", "Memory B", "CD8+ T", "CD4+ T", 
                    "Monocyte","ORACLE", "WTGS", "PACEG")
# b) heatmap of risk score  ----------------------------------------------------
col_fun <- colorRamp2(c(0.5, 2.5), c("white", "#2171b5"))
p1 <- Heatmap(as.matrix(data),  cluster_rows = F, 
                heatmap_legend_param = list(title = "Risk Score", 
                                            at = c(0,2,4,6,8)),
                cluster_columns = F,col = col_fun,
              column_title_gp = gpar(fontsize = 15),
              width = ncol(data)*unit(10, "mm"), 
              height = nrow(data)*unit(5, "mm"))

pdf(myFig1,
    width = 8,  
    height = 12)
print(p1)
dev.off()
# b) heatmap of correlation  ----------------------------------------------------

cor.all <- cor(data, method ="pearson")
idx <- c( "WTGS", "PACEG","ORACLE","Naive B", "Memory B", "CD8+ T", "CD4+ T", 
          "Monocyte")
cor.all <- cor.all[idx,idx]

col_fun <- colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027"))
p2 <- Heatmap(cor.all, name = "PCC", show_column_names = T,
              show_row_names = T, cluster_rows = F,cluster_columns = F,
              column_title_gp = gpar(fontsize = 15),
              width = ncol(cor.all)*unit(10, "mm"), 
              height = nrow(cor.all)*unit(10, "mm"),
              col = col_fun)

pdf(myFig2,
    width = 8,  
    height = 8)
print(p2)
dev.off()

# [2] CPRIT-MIRA ##################################################################
# [2.1] Integrate score & Predict Region Specific Risk =========================
rm(list=ls())

library(survival)


indir1 <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig3/"
myinf1 <- paste0(indir1,"/ORACLE_M2_Region_CPRIT-MIRA.txt")
myinf2 <- paste0(indir1,"/WTGS_M2_Region_CPRIT-MIRA.txt")
myinf3 <- paste0(indir1,"/PACEG_M2_Region_CPRIT-MIRA.txt")

indir2 <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig4/"
myinf4 <- paste0(indir2,"/RNAseq_Genes_representativeImm_CPRIT-MIRA.txt")

myinf5 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/Clinical_Sample_info/CohortC_CPRIT-MIRA_clinical_info.csv"

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig5/"

myoutf1 <- paste0(outdir,"/Region_specific_scores_CPRIT-MIRA.txt")

# a) ORACLE --------------------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1)
oracle <- data
# c) WTGS ----------------------------------------------------------------------
data <- read.table(myinf2, sep="\t", header=T, row.names=1)
wtgs <- data
# b) PACEG ---------------------------------------------------------------------
data <- read.table(myinf3, sep="\t", header=T, row.names=1)
paceg <- data

# d) immune cell ---------------------------------------------------------------
data <- read.table(myinf4, sep="\t", header=T, row.names=1)
cnum <- ncol(data)/4
data <- data[, 1:cnum] - data[, (cnum+1):(2*cnum)]
tmp <- colnames(data)
tmp <- gsub("_up\\.ES", "", tmp)
colnames(data) <- tmp
# since NK is not significant -- remove it
xx <- data[, -5]
imm <- xx

# e) integrate -----------------------------------------------------------------
comxx <- intersect(row.names(imm), row.names(oracle))
comxx <- intersect(row.names(wtgs), comxx)
comxx <- intersect(row.names(paceg), comxx)

oracle <- oracle[comxx,1]
paceg <- paceg[comxx,1]
wtgs <- wtgs[comxx,1]
imm <- imm[comxx,]

data <- cbind(imm, oracle, paceg, wtgs)


# f) Use max/min to construct model at the patient level -----------------------
colnames(data)
# [2] "NavB"  "MemB"  "CD8T"     "CD4T"  "Monocyte" "oracle"   "paceg"   "wtgs"    
flag <- c(-1,-1,-1,-1, 1, 1, 1, -1)

pid <- substr(row.names(data), 1, 6)
mypat <- unique(pid)

mydat <- matrix(0, length(mypat), ncol(data))
row.names(mydat) <- mypat
colnames(mydat) <- colnames(data)
mydat <- as.data.frame(mydat)

for(k in 1:length(mypat)) {
  se <- which(pid==mypat[k])
  tmp <- data[se, , drop=F]
  for(i in 1:ncol(data)) {
    if(flag[i]==-1) {
      mydat[k, i] <- min(tmp[,i])
    }
    if(flag[i]==1) {
      mydat[k, i] <- max(tmp[,i])
    }
    if(flag[i]==0) {
      mydat[k, i] <- mean(tmp[,i])
    } 
  }
}


# g) survival info -------------------------------------------------------------
info <- read.csv(myinf5, header=T, row.names=1)

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

comxx <- intersect(row.names(mydat), row.names(info))
info <- info[comxx,]
mydat <- mydat[comxx,]




# i) predict risk score for each regions ---------------------------------------
score <- data
all(colnames(score) == colnames(mydat)) # T
for(k in 1:ncol(data)) {
  mytf <- mydat[,k]
  xx <- cbind(info, mytf)
  mycox <- coxph(Surv(t.surv, e.surv)~mytf, xx)
  tmp.dat <- data[,k,drop=F]
  colnames(tmp.dat)="mytf"
  xx <- predict(mycox, newdata= tmp.dat, type = "risk")
  score[,k] <- xx
}
score[1:3,]

# j) export --------------------------------------------------------------------
write.table(score, myoutf1, sep="\t", quote=F)

# [2.2] Correlation of Region Specific Risk ====================================
rm(list=ls())

library(ComplexHeatmap)
library(circlize)

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig5/"

myinf1 <- paste0(outdir,"/Region_specific_scores_CPRIT-MIRA.txt")
myFig1 <- paste0(outdir,"/Region_risk_scores_CPRIT-MIRA.pdf")
myFig2 <- paste0(outdir,"/Correlation_CPRIT-MIRA.pdf")



# a) choose samples with multiple regions (>=3) --------------------------------
data <- read.table(myinf1, sep="\t", header=T, stringsAsFactors=F, 
                   quote="", comment.char="")

tag <- substr(row.names(data), 1, 6)
xx <- table(tag)
mul.sam  <- names(xx)[xx>=3]

se <- which(tag %in% mul.sam)

data <- data[se, ]
data <- round(data, 4)
colnames(data) <- c("Naive B", "Memory B", "CD8+ T", "CD4+ T", 
                    "Monocyte","ORACLE", "WTGS", "PACEG")
# b) heatmap of risk score  ----------------------------------------------------
col_fun <- colorRamp2(c(0.5, 2.5), c("white", "#2171b5"))
p1 <- Heatmap(as.matrix(data),  cluster_rows = F, 
              heatmap_legend_param = list(title = "Risk Score", 
                                          at = c(0,2,4,6,8)),
              cluster_columns = F,col = col_fun,
              column_title_gp = gpar(fontsize = 15),
              width = ncol(data)*unit(10, "mm"), 
              height = nrow(data)*unit(5, "mm"))

pdf(myFig1,
    width = 8,  
    height = 12)
print(p1)
dev.off()
# b) heatmap of correlation  ----------------------------------------------------

cor.all <- cor(data, method ="pearson")
idx <- c( "WTGS", "PACEG","ORACLE","Naive B", "Memory B", "CD8+ T", "CD4+ T", 
          "Monocyte")
cor.all <- cor.all[idx,idx]

col_fun <- colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027"))
p2 <- Heatmap(cor.all, name = "PCC", show_column_names = T,
              show_row_names = T, cluster_rows = F,cluster_columns = F,
              column_title_gp = gpar(fontsize = 15),
              width = ncol(cor.all)*unit(10, "mm"), 
              height = nrow(cor.all)*unit(10, "mm"),
              col = col_fun)

pdf(myFig2,
    width = 8,  
    height = 8)
print(p2)
dev.off()