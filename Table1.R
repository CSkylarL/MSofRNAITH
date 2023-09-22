#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Table 1 : Try Other signatures
# Author: Chenyang Li
# TRACERx Color: "#67000d",[ "#d73027" ], "#c51b7d"
# CPRIT-MIRA Color: "#08306b", [ "#3288bd" ], "#5e4fa2" 
# CPRIT-MIRA  Result is not significant --> show C-index
# 01/14/2023
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Collected_Gene_Signatures   #################################################

# [I] LUAD Training ============================================================
library(survival)
#1. Load data LUAD -------------------------------------------------------------
rm(list = ls())
outdir <- "./Table1/"
myinf0 <- paste0(outdir,"/Collected_Gene_Signatures.csv")
myinf1 <- "~/Mydata/TCGA_LUAD/LUAD_RNAseqv2_ALL_Symbol.rda"
myinf2 <- "~/Mydata/TCGA_LUAD/LUAD_Clincial_info.txt"

myoutf1 <- paste0(outdir,"/Collected_Gene_Signatures_LUAD_training_summary.csv")
myoutf2 <- paste0(outdir,"/Collected_Gene_Signatures_coefficients.csv")


sig.sum <- read.csv(myinf0)
sig.name <- unique(sig.sum$signature)


load(myinf1)
rna <- mydata
rm(mydata)
# keep primary samples
xx <- as.numeric(substr(colnames(rna), 14, 15))
table(xx)
# xx
# 1   2  11 
# 515   2  59 
# 01 Primary Solid Tumor
# 02 Recurrent Solid Tumor
# 11 Solid Tissue Normal

se <- which(xx==1)		## Primary Solid Tumor
length(se) # 515 samples
rna <- rna[,se]
dim(rna) # 20501   515
colnames(rna) <- substr(colnames(rna), 1, 12)

# info 
info <- read.table(myinf2, sep="\t", header=T, row.names=1, quote="",fill = TRUE)
se <- c("vital_status", "days_to_death", "days_to_last_followup")
info <- info[,se]
event <- ifelse(info$vital_status=="dead", 1, 0)
xx <- info[, "days_to_death"] 
time <-  ifelse(!is.na(xx), xx, info[, "days_to_last_followup"])

info <- data.frame(row.names = row.names(info),
                   e.surv= event,
                   t.surv = as.numeric(time))

comxx <- intersect(rownames(info) ,colnames(rna))
rna <- rna[,comxx]
info <- info[comxx,]


rna <- log2(rna+1)

all(rownames(info) == colnames(rna))
all(rownames(info) == rownames(t(rna))) # T
rownames(rna) <- gsub(rownames(rna), pattern = "-",replacement = ".")

# traing in LUAD
Training_sum <- data.frame(row.names = sig.name, 
                           "Signature" = sig.name)
sig_coeff <- data.frame()

for ( iii in sig.name) {
  se <- which(sig.sum$signature == iii)
  mysig <- sig.sum[se,]
  mygene <- gsub(mysig$Gene, pattern = "-",replacement = ".",fixed = T)
  
  se <- which(rownames(rna) %in% mygene)
  cat("\n", iii, "\t",length(se), "/",  length(mygene))
  Training_sum[iii,"No.Gene"] <- paste0(length(se), "/",  length(mygene))
  
  if (length(se) > 2) {
    coxData <- cbind(info,t(rna[se,]))
    mygene <- rownames(rna)[se]
    fmla <- as.formula(paste("Surv(t.surv, e.surv) ~", 
                             paste(mygene,collapse= "+")))
    
    coxModel <- coxph(formula = fmla,
                      data    = coxData)
    mycox <- summary(coxModel)
    Training_sum[iii, "P-value_logtest" ] <- mycox$logtest[3]
    Training_sum[iii, "P-value_sctest"  ] <- mycox$sctest[3]
    Training_sum[iii, "P-value_waldtest"] <- mycox$waldtest[3]
    Training_sum[iii,"C-index" ] <- round(mycox$concordance[1],3)
    mytmp <- data.frame(Gene = gsub(rownames(mycox$coefficients), 
                                    pattern = ".",replacement = "-",fixed = T),
                        beta = mycox$coefficients[,1],
                        signature = rep(iii,nrow(mycox$coefficients)))
    sig_coeff <- rbind(sig_coeff,mytmp)
  } else next
}

write.csv(file = myoutf1,x = Training_sum,quote = F,row.names = F)
write.csv(file = myoutf2,x = sig_coeff,quote = F,row.names = F)


# [II] Applied to Multiregional RNA-seq ========================================
rm(list=ls())
library(survival)
library(forestplot)
library(ggplot2)
outdir <- "./Table1/"

myinf0 <- paste0(outdir,"/Collected_Gene_Signatures_coefficients.csv")
myoutf0 <- paste0(outdir,"/Collected_Gene_Signatures_summary.csv")
sig.sum <- read.csv(myinf0)

sig.name <- unique(sig.sum$signature)

myres <- data.frame("HR","P-Value","C-index","cohort")
for ( iii in sig.name) {
  se <- which(sig.sum$signature == iii)
  mysig <- sig.sum[se,]
  
  # [1] TRACERx ================================================================
  # [1.1] mysig Score ----------------------------------------------------------
  myinf1 <- "~/Mydata/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
  myinf2 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"
  
  myoutf1 <- paste0(outdir,"/Collected_Gene_Signatures_", iii ,"_scores_TRACERx.txt")
  # a) expression and signature  -------------------------------------------------
  data <- read.table(myinf1, sep="\t", header=T, check.names=F)
  
  # confirm the same order of gene
  row.names(mysig) <- mysig$Gene
  comxx <- intersect(row.names(data), row.names(mysig))
  data <- data[comxx,]
  mysig <- mysig[comxx,]
  my_gene_number <- length(comxx)
  all(rownames(data) == mysig$Gene) # T
  
  
  mysig_coefficient <- mysig$beta
  names(mysig_coefficient) <- mysig$Gene
  # b) log transform -------------------------------------------------------------
  
  data <- log2(data+1)
  # c) survival info -------------------------------------------------------------
  info <- read.csv(myinf2, header=T, row.names=1)
  
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
  
  dim(data) 
  
  # d) M1 Trasformed gene expression score ---------------------------------------
  tmp <- matrix(0, nrow(data), length(mypat))
  row.names(tmp) <- row.names(data)
  all(row.names(tmp)  == names(mysig_coefficient)) # T
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
    dat.adj[,k] <- ifelse(test = mysig_coefficient > 0, 
                          yes = dat.max[,k],  no = dat.min[,k])
    dat.rev[,k] <- ifelse(test = mysig_coefficient < 0, 
                          yes = dat.max[,k],  no = dat.min[,k])
  }
  
  # patient score 
  res.m1 <- data.frame(m1.s.avg=t(dat.avg) %*% mysig_coefficient, 
                       m1.s.max=t(dat.max) %*% mysig_coefficient, 
                       m1.s.min=t(dat.min) %*% mysig_coefficient,
                       m1.s.adj=t(dat.adj) %*% mysig_coefficient,
                       m1.s.rev=t(dat.rev) %*% mysig_coefficient )
  
  # e) M2 Region Specific Score ---------------------------------------
  # region score
  all(row.names(data)  == names(mysig_coefficient)) # T
  signature <- t(data) %*% mysig$beta
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
  # [1.2] mysig Survival ---------------------------------------------------------
  myinf1 <- paste0(outdir,"/Collected_Gene_Signatures_", iii ,"_scores_TRACERx.txt")
  
  # a) load data -----------------------------------------------------------------
  rm(data)
  data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
  
  
  comxx <- intersect(row.names(data), row.names(info))
  info <- info[comxx,]
  data <- data[comxx,]
  
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
  
  mytmp <- p.data[1:3]
  # The last row show info: signature name, cohort, signature gene number
  mytmp[nrow(mytmp)+1,] <- c(paste0("Collected_Gene_Signatures_", iii ),
                             "TRACERx",
                             paste0(my_gene_number,"/", nrow(mysig) )
  )
  myres <- cbind(mytmp,myres)
  
  
  
  # [2] CPRIT-MIRA ===============================================================
  
  # [2.1] mysig Score ------------------------------------------------------------
  rm(data,info,p.data,mytmp,mydata,myinf1,myinf2, myoutf1)
  myinf1 <-  "./data/MDAMPLC.rda"
  
  myoutf1 <- paste0(outdir,"/Collected_Gene_Signatures_", iii ,"_scores_CPRIT-MIRA.txt")
  # a) expression and signature  -------------------------------------------------
  load(myinf1)
  data <- RNAseq
  
  # confirm the same order of gene
  row.names(mysig) <- mysig$Gene
  comxx <- intersect(row.names(data), row.names(mysig))
  data <- data[comxx,]
  mysig <- mysig[comxx,]
  my_gene_number <- length(comxx)
  all(rownames(data) == mysig$Gene) # T
  
  
  mysig_coefficient <- mysig$beta
  names(mysig_coefficient) <- mysig$Gene
  # b) log transform -------------------------------------------------------------
  
  data <- log2(data+1)
  # c) survival info -------------------------------------------------------------
  info <- Clinical.info
  colnames(info) <-  c("t.surv","e.surv", 
                      "Smoker.current.never.former",
                      "Gender.Male.Female" ,         
                      "Diagnosis.Age"  ,
                      "Stage" ,
                      "Subtype"         )
  
  pid <-  substr(colnames(data), 1, 6)
  comxx <- intersect(row.names(info), pid)
  info <- info[comxx,]
  se <- which(pid %in% comxx)
  # patient ID of expression data
  pid <- pid[se]
  data <- data[,se]
  # unique patient ID of survival data
  mypat <- row.names(info)
  
  dim(data) 
  
  # d) M1 Trasformed gene expression score ---------------------------------------
  tmp <- matrix(0, nrow(data), length(mypat))
  row.names(tmp) <- row.names(data)
  all(row.names(tmp)  == names(mysig_coefficient)) # T
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
    dat.adj[,k] <- ifelse(test = mysig_coefficient > 0, 
                          yes = dat.max[,k],  no = dat.min[,k])
    dat.rev[,k] <- ifelse(test = mysig_coefficient < 0, 
                          yes = dat.max[,k],  no = dat.min[,k])
  }
  
  # patient score 
  res.m1 <- data.frame(m1.s.avg=t(dat.avg) %*% mysig_coefficient, 
                       m1.s.max=t(dat.max) %*% mysig_coefficient, 
                       m1.s.min=t(dat.min) %*% mysig_coefficient,
                       m1.s.adj=t(dat.adj) %*% mysig_coefficient,
                       m1.s.rev=t(dat.rev) %*% mysig_coefficient )
  
  # e) M2 Region Specific Score ---------------------------------------
  # region score
  all(row.names(data)  == names(mysig_coefficient)) # T
  signature <- t(data) %*% mysig$beta
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
  # [2.2] mysig Survival ---------------------------------------------------------
  
  myinf1 <- paste0(outdir,"/Collected_Gene_Signatures_", iii ,"_scores_CPRIT-MIRA.txt")
  
  # a) load data -----------------------------------------------------------------
  rm(data)
  data <- read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
  
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
  
  # The C-index will be manually corrected 
  # since there are both hazardous and protective signatures
  # p.data$`C-index` <- ifelse(p.data$HR < 1, 1-p.data$`C-index`, p.data$`C-index`)
  
  mytmp <- p.data[1:3]
  # The last row show info: signature name, cohort, signature gene number
  mytmp[nrow(mytmp)+1,] <- c(paste0("Collected_Gene_Signatures_", iii ),
                             "CPRIT-MIRA",
                             paste0(my_gene_number,"/", nrow(mysig) )
  )
  myres <- cbind(myres,mytmp)
  
}

write.csv(file = myoutf0,x = myres,quote = F,row.names = T)
