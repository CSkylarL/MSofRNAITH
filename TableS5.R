#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Table S5: compare PACEG+Adj. to key clinical prognostic variables 
# Author: Chenyang Li
# 09/20/2023
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TRACERx421_LUAD #############################################################

rm(list=ls())
library(survival)
library(dplyr)

myinf1 <- "~/Mydata/TRACERx421_zenodo.7819449/20221014_transcriptomic_DATA/20221109_TRACERx421_all_patient_df.rds"
myinf2 <-  "./Supp_Fig2/PACEG_scores_TRACERx421_LUAD.txt"

# survival info -------------------------------------------------------------
info <- readRDS(myinf1)
rownames(info) <- info$cruk_id
# our model risk score----------------------------------------------------------
mymodel <- read.table(myinf2, sep="\t", header=T, check.names=F)

# filter info (with RNA-seq) ------------------------------------------------
comxx <- intersect(row.names(info), rownames(mymodel))
info <- info[comxx,]
mymodel <- mymodel[comxx,]
# preprocess info ------------------------------------------------

t.surv <- info$dfs_time_any_event * 0.0329 # Change day to month
e.surv <- info$cens_dfs_any_event
PACEG_Adj <- as.numeric(mymodel$m1.s.adj)
Smoke_Status <- info$smoking_status_merged %>% as.factor()
Stage <- info$pathologyTNM%>% as.factor()
Age <- info$age %>% as.numeric()
Gender <- info$sex %>% as.factor()

info_new <- data.frame(t.surv, e.surv,PACEG_Adj, Smoke_Status,Stage,Age,Gender)
rownames(info_new) <- rownames(info)
# Multivariate Cox regression  --------------------------------------------------

myformular <- as.formula(paste0("Surv(t.surv,e.surv) ~ ",
                                paste0(colnames(info_new[3:ncol(info_new)]), collapse = "+")))
cox_model <- coxph(myformular, data = info_new)
print(summary(cox_model))

# Call:
#   coxph(formula = myformular, data = info_new)
# 
# n= 187, number of events= 95 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# PACEG_Adj             0.52462   1.68982  0.15728 3.336    0.000851 ***
#   Smoke_StatusEx-Smoker 0.96549   2.62609  0.47247 2.044  0.041001 *  
#   Smoke_StatusSmoker    0.49960   1.64806  0.50091 0.997  0.318585    
# StageIB               0.36398   1.43905  0.36968 0.985    0.324830    
# StageIIA              0.94190   2.56486  0.38929 2.420    0.015539 *  
#   StageIIB              0.71191   2.03788  0.43350 1.642  0.100541    
# StageIIIA             1.10261   3.01203  0.36233 3.043    0.002342 ** 
#   StageIIIB             1.58581   4.88322  1.05666 1.501  0.133415    
# Age                   0.02598   1.02632  0.01277 2.034    0.041944 *  
#   GenderMale            0.06264   1.06464  0.22549 0.278  0.781166     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# PACEG_Adj                 1.690     0.5918    1.2415     2.300
# Smoke_StatusEx-Smoker     2.626     0.3808    1.0403     6.629
# Smoke_StatusSmoker        1.648     0.6068    0.6174     4.399
# StageIB                   1.439     0.6949    0.6973     2.970
# StageIIA                  2.565     0.3899    1.1959     5.501
# StageIIB                  2.038     0.4907    0.8713     4.766
# StageIIIA                 3.012     0.3320    1.4806     6.127
# StageIIIB                 4.883     0.2048    0.6156    38.738
# Age                       1.026     0.9744    1.0009     1.052
# GenderMale                1.065     0.9393    0.6843     1.656
# 
# Concordance= 0.728  (se = 0.028 )
# Likelihood ratio test= 56.84  on 10 df,   p=1e-08
# Wald test            = 52.45  on 10 df,   p=9e-08
# Score (logrank) test = 56.84  on 10 df,   p=1e-08

# This will be used for Table S5