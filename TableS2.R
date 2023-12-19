#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Table S2 : Try Other signatures (Table 1 supp)
# Author: Chenyang Li
# 09/20/2023
# TRACERx 421
# Additional validation: TRACERx 421 excluding TRACERx 100 (TRACERx421-100)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#  Collected_Gene_Signatures   #################################################
# [II] TRACERx421_nonTRACERx100 ========================================================
rm(list=ls())
library(fst)
library(survival)
library(dplyr)
outdir <- "./TableS2/"
dir.create(outdir)

myinf0 <- "./Table1/Collected_Gene_Signatures_coefficients.csv"
myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_tumour_df.rds"
myinf3 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-14_clinicohistopathological_data.fst"
myinf4 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"


myoutf0 <- paste0(outdir,"/Collected_Gene_Signatures_summary_TRACERx421_nonTRACERx100.csv")

# data -------------------------------------------------------------------------

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

# rna 
se <- which(colnames(rna) %in% region_LUAD$sample_name_cruk)
rna <- rna[se]
rna <- log2(rna+1)
# survival info

e.surv <- info$cens_dfs_any_event 
t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
info <- cbind(t.surv, e.surv, info)

pid <- substr(colnames(rna), 1, 8)
comxx <- intersect(row.names(info), pid)
info <- info[comxx,]
dim(info) # 261  52

se <- which(pid %in% comxx)
# patient ID of expression data
pid <- pid[se]
rna <- rna[,se]
# unique patient ID of survival data
mypat <- row.names(info)

dim(rna) #   28073   652

# Signature --------------------------------------------------------------------
sig.sum <- read.csv(myinf0)

# order 
sig.name <- c("Boutros2008","Krzystanek2016" ,"Bianchi2007","Kratz2012","Zhu2010" ,
              "Garber2001", "Wistuba2013","Raz2008", "Beer2002" )

# multiregional survival analysis ----------------------------------------------
myres <- data.frame()
for ( iii in sig.name) {
  se <- which(sig.sum$signature == iii)
  mysig <- sig.sum[se,]
  
  #  mysig Score ----------------------------------------------------------
  
  myoutf1 <- paste0(outdir,"/Collected_Gene_Signatures_", iii ,"_scores_TRACERx421_nonTRACERx100.txt")
  # a) expression and signature  -------------------------------------------------
  data <- rna
  
  # confirm the same order of gene
  row.names(mysig) <- mysig$Gene
  comxx <- intersect(row.names(data), row.names(mysig))
  data <- data[comxx,]
  mysig <- mysig[comxx,]
  my_gene_number <- length(comxx)
  all(rownames(data) == mysig$Gene) # T
  
  
  mysig_coefficient <- mysig$beta
  names(mysig_coefficient) <- mysig$Gene
  
  # b) M1 Trasformed gene expression score ---------------------------------------
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
  
  # a) load data -----------------------------------------------------------------
  data <- res
  
  
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
    p.data$`P-Value`[k] <- round(mycox$coefficients[5],8)
    p.data$`C-index`[k] <- round(mycox$concordance[1],3)
    
  }
  
  mytmp <- data.frame("Cohort" = rep("TRACERx421_nonTRACERx100",8),
                      "Signature" = rep(iii,8),
                      "N" = rep(paste0(my_gene_number,"/", nrow(mysig) ),8))
  mytmp <- cbind(mytmp,p.data[1:3])
  
  myres <- rbind(myres,mytmp)
  
}

write.csv(myres,myoutf0, quote = F, row.names = T)  
