#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.5
# Author: Chao Cheng; Chenyang Li
# Note: CPRIT-MIRA is the MDAMPLC cohort
# 12/14/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx ##################################################################
# [1.1] Integrate score & Predict Region Specific Risk =========================
rm(list=ls())

library(survival)


indir1 <- "./Fig3/"
myinf1 <- paste0(indir1,"/ORACLE_M2_Region_TRACERx.txt")
myinf2 <- paste0(indir1,"/WTGS_M2_Region_TRACERx.txt")
myinf3 <- paste0(indir1,"/PACEG_M2_Region_TRACERx.txt")

indir2 <- "./Fig4/"
myinf4 <- paste0(indir2,"/RNAseq_Genes_representativeImm_TRACERx.txt")

myinf5 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/RNAseq_Sample_FullID.txt"
myinf6 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

outdir <- "./Fig5/"

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

data <- cbind(imm, oracle,  wtgs, paceg)

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
library(dplyr)

library(ComplexHeatmap)
library(circlize)

outdir <- "./Fig5/"

myinf1 <- paste0(outdir,"/Region_specific_scores_TRACERx.txt")
myFig1 <- paste0(outdir,"/Region_risk_scores_TRACERx.pdf")
myFig2 <- paste0(outdir,"/Correlation_TRACERx.pdf")



# a) choose samples with multiple regions (>=3) --------------------------------
data <- read.table(myinf1, sep="\t", header=T, stringsAsFactors=F, 
                   quote="", comment.char="")
data <- data[order(rownames(data),decreasing = F),]
data <- apply(data, 2, function(x) (x-median(x))/(max(x)-min(x)))
tag <- substr(row.names(data), 1, 8)
xx <- table(tag)
mul.sam  <- names(xx)[xx>=3]

se <- which(tag %in% mul.sam)

data <- data[se, ]
data <- round(data, 4)

# b) heatmap of risk score  ----------------------------------------------------
rowAnn <- data.frame(
  row.names = rownames(data),
  patient = rownames(data) %>% 
    strsplit("_") %>% 
    lapply("[[", 1) %>% 
    unlist() ,
  region = rownames(data) %>% 
    strsplit("_") %>% 
    lapply("[[", 2) %>% 
    unlist() 
)

rowAnn$pa_col <- ifelse (test = as.numeric(as.factor(rowAnn$patient)) %% 2 == 1,
                         yes = "#9970ab",no = "#5aae61")
cellAnn <- cbind(rowAnn,data) %>% 
  group_by(patient) %>% 
  mutate(NavB_rank =rank(-NavB ,ties.method = "min"),
         MemB_rank =rank(-MemB ,ties.method = "min"),
         CD8T_rank =rank(-CD8T ,ties.method = "min"),
         CD4T_rank =rank(-CD4T ,ties.method = "min"),
         Monocyte_rank =rank(-Monocyte ,ties.method = "min"),
         oracle_rank =rank(-oracle ,ties.method = "min"),
         wtgs_rank =rank(-wtgs ,ties.method = "min"),
         paceg_rank =rank(-paceg ,ties.method = "min"),
         
  ) %>%
  ungroup() %>%
  select(NavB_rank, MemB_rank, CD8T_rank, CD4T_rank,
         Monocyte_rank, oracle_rank, wtgs_rank, paceg_rank) %>%
  data.frame()

colnames(data)
colnames(data) <- c("Naive B", "Memory B", "CD8+ T", "CD4+ T", 
                    "Monocyte","ORACLE", "WTGS", "PACEG")

p1 <-  Heatmap(as.matrix(data),  cluster_rows = F, 
               heatmap_legend_param = list(title = "Risk Score"), 
               cluster_columns = F,
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 
                 if(cellAnn[i, j] == 1) {
                   gb = textGrob("*")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("*", x, y - gb_h*0.6 + gb_w*0.4,
                             gp = gpar(fontsize = 15))
                 } else if(cellAnn[i, j] == 2) {
                   gb = textGrob("+")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("+", x, y - gb_h*0.2 + gb_w*0.4,
                             gp = gpar(fontsize = 12))
                 }
               },
               border_gp = gpar(col = "black", lty = 1),
               rect_gp = gpar(col = "white", lwd = 0.5),
               row_split = rowAnn$patient,
               left_annotation = 
                 rowAnnotation(
                   Patient = anno_block(gp = 
                                          gpar(fill = rep(c("#1f78b480", 
                                                            "#33a02c80"),4)),
                                        labels = unique(rowAnn$patient), 
                                        labels_gp = gpar(col = "white", 
                                                         fontsize = 6))),
               row_title=NULL,
               right_annotation = 
                 rowAnnotation(
                   Region= anno_text(rowAnn$region,
                                     just = "center", 
                                     location = unit(0.5, "npc"), 
                                     show_name = F)),
               show_row_names = F,
               
               row_gap = unit(1, "mm"),
               column_title_gp = gpar(fontsize = 15),
               width = ncol(data)*unit(10, "mm"),
               height = nrow(data)*unit(4, "mm"))

pdf(myFig1,
    width = 6,  
    height = 20)
print(p1)
dev.off()
# b) heatmap of correlation  ----------------------------------------------------

cor.all <- cor(data, method ="pearson")
idx <- c( "WTGS", "PACEG","ORACLE","Naive B", "Memory B", "CD8+ T", "CD4+ T", 
          "Monocyte")
cor.all <- cor.all[idx,idx]

p2 <- Heatmap(cor.all, name = "PCC", show_column_names = T,
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              border_gp = gpar(col = "black", lty = 1),
              show_row_names = T, cluster_rows = T,cluster_columns = T,
              column_names_rot = 60,
              column_title_gp = gpar(fontsize = 15),
              width = ncol(cor.all)*unit(10, "mm"), 
              height = nrow(cor.all)*unit(10, "mm"))

pdf(myFig2,
    width = 8,  
    height = 8)
print(p2)
dev.off()

# [1.3] Select patients to show in main figures ================================

rm(list=ls())
library(dplyr)

library(ComplexHeatmap)
library(circlize)

outdir <- "./Fig5/"

myinf1 <- paste0(outdir,"/Region_specific_scores_TRACERx.txt")
myFig1 <- paste0(outdir,"/Region_risk_scores_TRACERx_select.pdf")



# a) select patients  ----------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, stringsAsFactors=F, 
                   quote="", comment.char="")
data <- data[order(rownames(data),decreasing = F),]
data <- apply(data, 2, function(x) (x-median(x))/(max(x)-min(x)))
tag <- substr(row.names(data), 1, 8)
xx <- table(tag)

# the patient are selected 
# because their highest risk scores are 
# from the same region using different signatures 
sel.sam  <- c("CRUK0005", "CRUK0017", "CRUK0018", "CRUK0035","CRUK0050", 
              "CRUK0052", "CRUK0084", "CRUK0086", "CRUK0094", "CRUK0098")

se <- which(tag %in% sel.sam)

data <- data[se, ]
data <- round(data, 4)

# b) heatmap of risk score  ----------------------------------------------------
rowAnn <- data.frame(
  row.names = rownames(data),
  patient = rownames(data) %>% 
    strsplit("_") %>% 
    lapply("[[", 1) %>% 
    unlist() ,
  region = rownames(data) %>% 
    strsplit("_") %>% 
    lapply("[[", 2) %>% 
    unlist() 
)

rowAnn$pa_col <- ifelse (test = as.numeric(as.factor(rowAnn$patient)) %% 2 == 1,
                         yes = "#9970ab",no = "#5aae61")
cellAnn <- cbind(rowAnn,data) %>% 
  group_by(patient) %>% 
  mutate(NavB_rank =rank(-NavB ,ties.method = "min"),
         MemB_rank =rank(-MemB ,ties.method = "min"),
         CD8T_rank =rank(-CD8T ,ties.method = "min"),
         CD4T_rank =rank(-CD4T ,ties.method = "min"),
         Monocyte_rank =rank(-Monocyte ,ties.method = "min"),
         oracle_rank =rank(-oracle ,ties.method = "min"),
         wtgs_rank =rank(-wtgs ,ties.method = "min"),
         paceg_rank =rank(-paceg ,ties.method = "min"),
         
  ) %>%
  ungroup() %>%
  select(NavB_rank, MemB_rank, CD8T_rank, CD4T_rank,
         Monocyte_rank, oracle_rank, wtgs_rank, paceg_rank) %>%
  data.frame()

colnames(data)
colnames(data) <- c("Naive B", "Memory B", "CD8+ T", "CD4+ T", 
                    "Monocyte","ORACLE", "WTGS", "PACEG")
p1 <-  Heatmap(as.matrix(data),  cluster_rows = F, 
               heatmap_legend_param = list(title = "Risk Score"), 
               cluster_columns = F,
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 
                 if(cellAnn[i, j] == 1) {
                   gb = textGrob("*")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("*", x, y - gb_h*0.6 + gb_w*0.4,
                             gp = gpar(fontsize = 15))
                 } else if(cellAnn[i, j] == 2) {
                   gb = textGrob("+")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("+", x, y - gb_h*0.2 + gb_w*0.4,
                             gp = gpar(fontsize = 12))
                 }
               },
               border_gp = gpar(col = "black", lty = 1),
               rect_gp = gpar(col = "white", lwd = 0.5),
               row_split = rowAnn$patient,
               left_annotation = 
                 rowAnnotation(
                   Patient = anno_block(gp = 
                                          gpar(fill = rep(c("#1f78b480", 
                                                            "#33a02c80"),4)),
                                        labels = unique(rowAnn$patient), 
                                        labels_gp = gpar(col = "white", 
                                                         fontsize = 7))),
               row_title=NULL,
               right_annotation = 
                 rowAnnotation(
                   Region= anno_text(rowAnn$region,
                                     just = "center", 
                                     location = unit(0.5, "npc"), 
                                     show_name = F)),
               show_row_names = F,
               
               row_gap = unit(1, "mm"),
               column_title_gp = gpar(fontsize = 15),
               width = ncol(data)*unit(10, "mm"),
               height = nrow(data)*unit(5, "mm"))

pdf(myFig1,
    width = 6,  
    height = 10)
print(p1)
dev.off()

# [2] CPRIT-MIRA ##################################################################
# [2.1] Integrate score & Predict Region Specific Risk =========================
rm(list=ls())

library(survival)


indir1 <- "./Fig3/"
myinf1 <- paste0(indir1,"/ORACLE_M2_Region_CPRIT-MIRA.txt")
myinf2 <- paste0(indir1,"/WTGS_M2_Region_CPRIT-MIRA.txt")
myinf3 <- paste0(indir1,"/PACEG_M2_Region_CPRIT-MIRA.txt")

indir2 <- "./Fig4/"
myinf4 <- paste0(indir2,"/RNAseq_Genes_representativeImm_CPRIT-MIRA.txt")

myinf5 <-"./data/MDAMPLC.rda"

outdir <- "./Fig5/"

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

data <- cbind(imm, oracle,  wtgs, paceg)


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
load(myinf5)

info <- Clinical.info
colnames(info) <-  c("t.surv","e.surv", 
                     "Smoker.current.never.former",
                     "Gender.Male.Female" ,         
                     "Diagnosis.Age"  ,
                     "Stage" ,
                     "Subtype"         )

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



# [1.2] Display Region Specific Risk ===========================================
rm(list=ls())
library(dplyr)

library(ComplexHeatmap)
library(circlize)

outdir <- "./Fig5/"

myinf1 <- paste0(outdir,"/Region_specific_scores_CPRIT-MIRA.txt")
myFig1 <- paste0(outdir,"/Region_risk_scores_CPRIT-MIRA.pdf")



# a) choose samples with multiple regions (>=5) --------------------------------
data <- read.table(myinf1, sep="\t", header=T, stringsAsFactors=F, 
                   quote="", comment.char="")
data <- data[order(rownames(data),decreasing = F),]
data <- apply(data, 2, function(x) (x-median(x))/(max(x)-min(x)))
tag <- substr(row.names(data), 1, 6)
xx <- table(tag)
mul.sam  <- names(xx)[xx>=3]

se <- which(tag %in% mul.sam)

data <- data[se, ]
data <- round(data, 4)

# b) heatmap of risk score  ----------------------------------------------------
rowAnn <- data.frame(
  row.names = rownames(data),
  patient = rownames(data) %>% 
    strsplit("-") %>% 
    lapply("[[", 1) %>% 
    unlist() ,
  region = rownames(data) %>% 
    strsplit("-") %>% 
    lapply("[[", 2) %>% 
    unlist() %>%
    gsub(pattern = "T",replacement = "R")
)

rowAnn$pa_col <- ifelse (test = as.numeric(as.factor(rowAnn$patient)) %% 2 == 1,
                         yes = "#9970ab",no = "#5aae61")
cellAnn <- cbind(rowAnn,data) %>% 
  group_by(patient) %>% 
  mutate(NavB_rank =rank(-NavB ,ties.method = "min"),
         MemB_rank =rank(-MemB ,ties.method = "min"),
         CD8T_rank =rank(-CD8T ,ties.method = "min"),
         CD4T_rank =rank(-CD4T ,ties.method = "min"),
         Monocyte_rank =rank(-Monocyte ,ties.method = "min"),
         oracle_rank =rank(-oracle ,ties.method = "min"),
         wtgs_rank =rank(-wtgs ,ties.method = "min"),
         paceg_rank =rank(-paceg ,ties.method = "min"),
         
  ) %>%
  ungroup() %>%
  select(NavB_rank, MemB_rank, CD8T_rank, CD4T_rank,
         Monocyte_rank, oracle_rank, wtgs_rank, paceg_rank) %>%
  data.frame()

colnames(data)
colnames(data) <- c("Naive B", "Memory B", "CD8+ T", "CD4+ T", 
                    "Monocyte","ORACLE", "WTGS", "PACEG")

p1 <-  Heatmap(as.matrix(data),  cluster_rows = F, 
               heatmap_legend_param = list(title = "Risk Score"), 
               cluster_columns = F,
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 
                 if(cellAnn[i, j] == 1) {
                   gb = textGrob("*")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("*", x, y - gb_h*0.6 + gb_w*0.4,
                             gp = gpar(fontsize = 15))
                 } else if(cellAnn[i, j] == 2) {
                   gb = textGrob("+")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("+", x, y - gb_h*0.2 + gb_w*0.4,
                             gp = gpar(fontsize = 12))
                 }
               },
               border_gp = gpar(col = "black", lty = 1),
               rect_gp = gpar(col = "white", lwd = 0.5),
               row_split = rowAnn$patient,
               left_annotation = 
                 rowAnnotation(
                   Patient = anno_block(gp = 
                                          gpar(fill = rep(c("#1f78b480", 
                                                            "#33a02c80"),4)),
                                        labels = unique(rowAnn$patient), 
                                        labels_gp = gpar(col = "white", 
                                                         fontsize = 8))),
               row_title=NULL,
               right_annotation = 
                 rowAnnotation(
                   Region= anno_text(rowAnn$region,
                                     just = "center", 
                                     location = unit(0.5, "npc"), 
                                     show_name = F)),
               show_row_names = F,
               
               row_gap = unit(1, "mm"),
               column_title_gp = gpar(fontsize = 15),
               width = ncol(data)*unit(10, "mm"),
               height = nrow(data)*unit(4, "mm"))

pdf(myFig1,
    width = 6,  
    height = 10)
print(p1)
dev.off()


# [1.3] Select patients to show in main figures ================================
rm(list=ls())
library(dplyr)

library(ComplexHeatmap)
library(circlize)

outdir <- "./Fig5/"

myinf1 <- paste0(outdir,"/Region_specific_scores_CPRIT-MIRA.txt")
myFig1 <- paste0(outdir,"/Region_risk_scores_CPRIT-MIRA_select.pdf")



# a) choose samples with multiple regions (>=5) --------------------------------
data <- read.table(myinf1, sep="\t", header=T, stringsAsFactors=F, 
                   quote="", comment.char="")
data <- data[order(rownames(data),decreasing = F),]
data <- apply(data, 2, function(x) (x-median(x))/(max(x)-min(x)))
tag <- substr(row.names(data), 1, 6)
xx <- table(tag)

# the patient are selected 
# because their highest risk scores are 
# from the same region using different signatures 
sel.sam  <- c("213109", "244136", "250854", "841808", "886403")

se <- which(tag %in% sel.sam)

data <- data[se, ]
data <- round(data, 4)

# b) heatmap of risk score  ----------------------------------------------------
rowAnn <- data.frame(
  row.names = rownames(data),
  patient = rownames(data) %>% 
    strsplit("-") %>% 
    lapply("[[", 1) %>% 
    unlist() ,
  region = rownames(data) %>% 
    strsplit("-") %>% 
    lapply("[[", 2) %>% 
    unlist() %>%
    gsub(pattern = "T",replacement = "R")
)

rowAnn$pa_col <- ifelse (test = as.numeric(as.factor(rowAnn$patient)) %% 2 == 1,
                         yes = "#9970ab",no = "#5aae61")
cellAnn <- cbind(rowAnn,data) %>% 
  group_by(patient) %>% 
  mutate(NavB_rank =rank(-NavB ,ties.method = "min"),
         MemB_rank =rank(-MemB ,ties.method = "min"),
         CD8T_rank =rank(-CD8T ,ties.method = "min"),
         CD4T_rank =rank(-CD4T ,ties.method = "min"),
         Monocyte_rank =rank(-Monocyte ,ties.method = "min"),
         oracle_rank =rank(-oracle ,ties.method = "min"),
         wtgs_rank =rank(-wtgs ,ties.method = "min"),
         paceg_rank =rank(-paceg ,ties.method = "min"),
         
  ) %>%
  ungroup() %>%
  select(NavB_rank, MemB_rank, CD8T_rank, CD4T_rank,
         Monocyte_rank, oracle_rank, wtgs_rank, paceg_rank) %>%
  data.frame()

colnames(data)
colnames(data) <- c("Naive B", "Memory B", "CD8+ T", "CD4+ T", 
                    "Monocyte","ORACLE", "WTGS", "PACEG")

p1 <-  Heatmap(as.matrix(data),  cluster_rows = F, 
               heatmap_legend_param = list(title = "Risk Score"), 
               cluster_columns = F,
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 
                 if(cellAnn[i, j] == 1) {
                   gb = textGrob("*")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("*", x, y - gb_h*0.6 + gb_w*0.4,
                             gp = gpar(fontsize = 15))
                 } else if(cellAnn[i, j] == 2) {
                   gb = textGrob("+")
                   gb_w = convertWidth(grobWidth(gb), "mm")
                   gb_h = convertHeight(grobHeight(gb), "mm")
                   grid.text("+", x, y - gb_h*0.2 + gb_w*0.4,
                             gp = gpar(fontsize = 12))
                 }
               },
               border_gp = gpar(col = "black", lty = 1),
               rect_gp = gpar(col = "white", lwd = 0.5),
               row_split = rowAnn$patient,
               left_annotation = 
                 rowAnnotation(
                   Patient = anno_block(gp = 
                                          gpar(fill = rep(c("#1f78b480", 
                                                            "#33a02c80"),4)),
                                        labels = unique(rowAnn$patient), 
                                        labels_gp = gpar(col = "white", 
                                                         fontsize = 8))),
               row_title=NULL,
               right_annotation = 
                 rowAnnotation(
                   Region= anno_text(rowAnn$region,
                                     just = "center", 
                                     location = unit(0.5, "npc"), 
                                     show_name = F)),
               show_row_names = F,
               
               row_gap = unit(1, "mm"),
               column_title_gp = gpar(fontsize = 15),
               width = ncol(data)*unit(10, "mm"),
               height = nrow(data)*unit(4, "mm"))

pdf(myFig1,
    width = 6,  
    height = 6)
print(p1)
dev.off()