#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.1
# Author: Chao Cheng; Chenyang Li
# 12/12/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx ##################################################################
# [1.1] Pairwise Similarity of Expression data =================================
rm(list=ls())
myinf0 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/Clinical_Sample_info/RNAseq_Sample_FullID.txt"
myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig1/"
myoutf1 <- paste0(outdir,"/Fig1_pairwise_distance_TRACERx.rda")

info <- read.table(myinf0, sep="\t", header=T, stringsAsFactors=F, quote="", comment.char="")

# a) calculate similarity ------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)

xx <- apply(data>0, 1, sum)
se <- which(xx>0)
data <- data[se,]
dim(data)
range(data)

data <- log10(data+1)

mysam <- colnames(data)

nn <- length(mysam)*(length(mysam)-1)/2
res <- matrix(0, nn, 3)
colnames(res) <- c("ID1", "ID2","E.Dist")
res <- as.data.frame(res)

# Euclidean Distance of gene expression 
cnum <- ncol(data)
count <- 0
for(i in 1:(cnum-1)) {
  tmpx <- as.numeric(data[,i])
  for(j in (i+1):cnum)
  {
    cat("\r\r\r", i, "->", j)
    count <- count + 1
    tmpy <- as.numeric(data[,j])
    res[count,1] <- mysam[i]
    res[count,2] <- mysam[j]
    xx <- sqrt(sum((tmpx-tmpy)^2))
    res[count,3] <- xx
  }
}
# b) map ID  -------------------------------------------------------------------
data <- res
sam1 <- substr(data$ID1, 17, 24)
sam2 <- substr(data$ID2, 17, 24)

map <- info[,1]
names(map) <- info[,2]

xx <- substr(data[,1], 1, 15)
tmp <- map[xx]
data[,1] <- tmp
xx <- substr(data[,2], 1, 15)
tmp <- map[xx]
data[,2] <- tmp

exp.mat <- data


# c) Export Pairwise Similarity data for Fig.1 ---------------------------------
save(exp.mat, file=myoutf1)
# [1.2] Compare RNA inter & intra heterogeneity ================================
rm(list=ls())
library(ggplot2)
library(dplyr)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig1/"
myinf1 <- paste0(outdir,"/Fig1_pairwise_distance_TRACERx.rda")
myFig1 <- paste0(outdir,"/Fig1_RNA-ITH_distribution_TRACERx.pdf")


load(myinf1)
ls()
# a) Normalize -----------------------------------------------------------------
data <- exp.mat
data[1:3,]
xx <- data$E.Dist
tmp <- (xx-min(xx))/(max(xx)-min(xx))
range(tmp)

data$myv <- tmp
tag1 <- substr(data$ID1, 1, 8)
tag2 <- substr(data$ID2, 1, 8)
sum(tag1==tag2)		## 194
sum(tag1!=tag2)		## 13172

se <- which(tag1==tag2)
myx1 <- data$myv[se] 
se <- which(tag1!=tag2)
myx2 <- data$myv[se] 

wilcox.test(myx1, myx2, alternative="l")$p.value			## 5.425435e-116

# b) Plot data -----------------------------------------------------------------

p.data <- rbind(
  data.frame(Value = myx1, Feature = rep("Intratumor", length(myx1))),
  data.frame(Value = myx2, Feature = rep("Intertumor", length(myx2)))
) %>%
  mutate(Feature = factor(Feature, 
                          levels = c("Intratumor", "Intertumor")),
         Value = as.numeric(Value))
# c) Violin plot ---------------------------------------------------------------

p <- ggplot(p.data, aes(x = Feature, y = Value, fill = Feature)) + 
  geom_violin(alpha = 0.6) +
  geom_boxplot(alpha = 0.8, width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("#d73027", "#c51b7d"))  + 
  annotate("text", x = 2.0, y = 0.02, size = 7, 
           label = "p < 0.0001", family = "serif") +
  guides(fill = "none") +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="",y="RNA Heterogeneity",
       title="TRACERx",
       base_size = 20, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 4.5/3,
        axis.text.x = element_text(face="bold", color="black", size=18),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=18),
        axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1, 0.25),expand = c(0,0)) 

# d) Export --------------------------------------------------------------------

pdf(myFig1, width = 4, height = 6)
print(p)
dev.off()

# [1.3] Patient Specific RNA-ITH ===============================================
rm(list=ls())
library(ggplot2)
library(dplyr)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig1/"
myinf1 <- paste0(outdir,"/Fig1_pairwise_distance_TRACERx.rda")
myFig1 <- paste0(outdir,"/Fig1_RNA-ITH_patient_TRACERx.pdf")

# a) Normalize -----------------------------------------------------------------

load(myinf1)
data <- exp.mat
xx <- data$E.Dist
xx <- (xx-min(xx))/(max(xx)-min(xx))
data$myv <- xx
# b) Only Keep multi-region samples --------------------------------------------
xx <- unique(c(data$ID1, data$ID2))
length(xx) #164 regions
yy <- substr(xx, 1, 8)
length(unique(yy)) # 64 patients
sort(table(yy))
# yy
# CRUK0002 CRUK0009 CRUK0010 CRUK0012 CRUK0014 CRUK0019 CRUK0026 CRUK0028 
# 1        1        1        1        1        1        1        1 
# CRUK0031 CRUK0032 CRUK0042 CRUK0044 CRUK0048 CRUK0053 CRUK0055 CRUK0060 
# 1        1        1        1        1        1        1        1 
# CRUK0075 CRUK0081 CRUK0097 CRUK0020 CRUK0021 CRUK0023 CRUK0025 CRUK0027 
# 1        1        1        2        2        2        2        2 
# CRUK0030 CRUK0033 CRUK0047 CRUK0057 CRUK0061 CRUK0067 CRUK0076 CRUK0078 
# 2        2        2        2        2        2        2        2 
# CRUK0082 CRUK0090 CRUK0095 CRUK0099 CRUK0004 CRUK0005 CRUK0017 CRUK0024 
# 2        2        2        2        3        3        3        3 
# CRUK0039 CRUK0052 CRUK0084 CRUK0086 CRUK0094 CRUK0098 CRUK0100 CRUK0018 
# 3        3        3        3        3        3        3        4 
# CRUK0035 CRUK0036 CRUK0041 CRUK0046 CRUK0062 CRUK0077 CRUK0079 CRUK0083 
# 4        4        4        4        4        4        4        4 
# CRUK0013 CRUK0037 CRUK0050 CRUK0069 CRUK0070 CRUK0096 CRUK0029 CRUK0065 
# 5        5        5        5        5        5        6        6 
# only keep the tumor has more than 1 region

tag1 <- substr(data$ID1, 1, 8)
tag2 <- substr(data$ID2, 1, 8)
se <- which(tag1==tag2)
data <- data[se,]
unique(substr(data$ID1, 1, 8))  %>% length() # 45 Patients
unique(c(data$ID1,data$ID2)) %>% length()  #145 regions
dim(data) # 194 4 
# c) mean of RNA-ITH per patient -----------------------------------------------

data$Patient <- substr(data$ID1, 1, 8)

myavg <- data %>%
  group_by(Patient) %>%
  summarise(avg = mean(myv)) %>%
  select(Patient, avg) %>%
  data.frame() %>% 
  arrange(-avg)

# d) Plot data -----------------------------------------------------------------

p.data <- data %>%
  merge(myavg, by= "Patient") %>%
  mutate(Patient = factor(Patient, levels = myavg$Patient)) %>%
  arrange(Patient)

# e) Dot plot ------------------------------------------------------------------

p <- ggplot() + 
  geom_point(data = p.data, aes(x=Patient, y=myv),
             size = 3, color = "#d73027", alpha = 0.6) + 
  geom_line(data = p.data, aes(x=Patient,  y=avg), size = 1, group = 1,
            linetype=1, col="#67000d") +     
  theme_bw(base_size = 12, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Patient ID",y="RNA-ITH",
       title="Patient-Specific RNA-ITH (TRACERx)") +
  theme(legend.position="none",
        aspect.ratio=1/2, 
        panel.border = element_rect(size = 1, fill = NA),
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "serif", size=12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",family = "serif", size=10,
                                   angle = 60,  vjust = 1, hjust=1),
        axis.text.y = element_text(face="bold", color="black", family = "serif",size=10),
        axis.title.x = element_text(face="bold", color="black",family = "serif", size=12),
        axis.title.y = element_text(face="bold",color="black",family = "serif", size=12)) +
  scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1, 0.25),expand = c(0,0)) 

pdf(myFig1, width= 10, height= 5)
print(p)
dev.off()





# [2] CPRIT-MIRA ##################################################################
# [2.1] Pairwise Similarity of Expression data =================================
rm(list=ls())
myinf1 <- "/rsrch3/scratch/genomic_med/cli15/JayZhang/geneset/CohortC_CPRIT-MIRA/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"

outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig1/"
myoutf1 <- paste0(outdir,"/Fig1_pairwise_distance_CPRIT-MIRA.rda")

# a) calculate similarity ------------------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)

xx <- apply(data>0, 1, sum)
se <- which(xx>0)
data <- data[se,]
dim(data)
range(data)

data <- log10(data+1)

mysam <- colnames(data)

nn <- length(mysam)*(length(mysam)-1)/2
res <- matrix(0, nn, 3)
colnames(res) <- c("ID1", "ID2","E.Dist")
res <- as.data.frame(res)

# Euclidean Distance of gene expression 
cnum <- ncol(data)
count <- 0
for(i in 1:(cnum-1)) {
  tmpx <- as.numeric(data[,i])
  for(j in (i+1):cnum)
  {
    cat("\r\r\r", i, "->", j)
    count <- count + 1
    tmpy <- as.numeric(data[,j])
    res[count,1] <- mysam[i]
    res[count,2] <- mysam[j]
    xx <- sqrt(sum((tmpx-tmpy)^2))
    res[count,3] <- xx
  }
}
# b)  remove normal ----------------------------------------------
data <- res

table(substr(data[,1], 11, 11) ) # N is normal
# N   T 
# 2918 3523 
se <- which(substr(data[,1], 11, 11) =="T"  &  substr(data[,2], 11, 11) =="T")
data <- data[se,]

exp.mat <- data

# c) Export Pairwise Similarity data for Fig.1 ---------------------------------
save(exp.mat, file=myoutf1)
# [2.2] Compare RNA inter & intra heterogeneity ================================
rm(list=ls())
library(ggplot2)
library(dplyr)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig1/"
myinf1 <- paste0(outdir,"/Fig1_pairwise_distance_CPRIT-MIRA.rda")
myFig1 <- paste0(outdir,"/Fig1_RNA-ITH_distribution_CPRIT-MIRA.pdf")


load(myinf1)
ls()
# a) Normalize -----------------------------------------------------------------
data <- exp.mat
data[1:3,]
xx <- data$E.Dist
tmp <- (xx-min(xx))/(max(xx)-min(xx))
range(tmp)

data$myv <- tmp
tag1 <- substr(data$ID1, 1, 8)
tag2 <- substr(data$ID2, 1, 8)
sum(tag1==tag2)		## 53
sum(tag1!=tag2)		##  1963

se <- which(tag1==tag2)
myx1 <- data$myv[se] 
se <- which(tag1!=tag2)
myx2 <- data$myv[se] 

wilcox.test(myx1, myx2, alternative="l")$p.value			##  1.632602e-24

# b) Plot data -----------------------------------------------------------------

p.data <- rbind(
  data.frame(Value = myx1, Feature = rep("Intratumor", length(myx1))),
  data.frame(Value = myx2, Feature = rep("Intertumor", length(myx2)))
) %>%
  mutate(Feature = factor(Feature, 
                          levels = c("Intratumor", "Intertumor")),
         Value = as.numeric(Value))
# c) Violin plot ---------------------------------------------------------------

p <- ggplot(p.data, aes(x = Feature, y = Value, fill = Feature)) + 
  geom_violin(alpha = 0.6) +
  geom_boxplot(alpha = 0.8, width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("#3288bd", "#5e4fa2"))  + 
  annotate("text", x = 1.0, y = 0.98, size = 7, 
           label = "p < 0.0001", family = "serif") +
  guides(fill = "none") +
  theme_bw(base_size = 20, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="",y="RNA Heterogeneity",
       title="CPRIT-MIRA",
       base_size = 20, base_family = "serif",face="bold") +
  theme(legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        legend.position="right",
        legend.box = "vertical",
        legend.direction= "vertical",
        panel.grid=element_blank(),
        legend.key.width = unit(0.5,"cm"),
        legend.title = element_text(face="bold", color="black",family = "serif", size=10),
        legend.text= element_text(face="bold", color="black",family = "serif", size=10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 4.5/3,
        axis.text.x = element_text(face="bold", color="black", size=18),
        axis.text.y = element_text(face="bold", color="black", size=15),
        axis.title.x = element_text(face="bold", color="black", size=18),
        axis.title.y = element_text(face="bold",color="black", size=17)) +
  scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1, 0.25),expand = c(0,0)) 

# d) Export --------------------------------------------------------------------

pdf(myFig1, width = 4, height = 6)
print(p)
dev.off()

# [2.3] Patient Specific RNA-ITH ===============================================
rm(list=ls())
library(ggplot2)
library(dplyr)
outdir <- "~/ChaoCheng/F_lung_heterogeneity/FB_region_compare/FB06_MS4/Fig1/"
myinf1 <- paste0(outdir,"/Fig1_pairwise_distance_CPRIT-MIRA.rda")
myFig1 <- paste0(outdir,"/Fig1_RNA-ITH_patient_CPRIT-MIRA.pdf")

# a) Normalize -----------------------------------------------------------------

load(myinf1)
data <- exp.mat
xx <- data$E.Dist
xx <- (xx-min(xx))/(max(xx)-min(xx))
data$myv <- xx
# b) Only Keep multi-region samples --------------------------------------------
xx <- unique(c(data$ID1, data$ID2))
length(xx) #64 regions
yy <- substr(xx, 1, 8)
length(unique(yy)) # 25 patients
sort(table(yy))
# yy
# 170547-T 191645-T 250813-T 317699-T 832461-T 834243-T 838180-T 840550-T 840674-T 
# 2        2        2        2        2        2        2        2        2 
# 852845-T 873708-T 213109-T 244136-T 250854-T 260856-T 270881-T 831868-T 833190-T 
# 2        2        3        3        3        3        3        3        3 
# 833625-T 841808-T 849616-T 858752-T 875099-T 886403-T 888712-T 
# 3        3        3        3        3        3        3 
# only keep the tumor has more than 1 region

tag1 <- substr(data$ID1, 1, 8)
tag2 <- substr(data$ID2, 1, 8)
se <- which(tag1==tag2)
data <- data[se,]
unique(substr(data$ID1, 1, 8))  %>% length() # 25 Patients
unique(c(data$ID1,data$ID2)) %>% length()  #64 regions

dim(data) # 53  4
# c) mean of RNA-ITH per patient -----------------------------------------------

data$Patient <- substr(data$ID1, 1, 8)

myavg <- data %>%
  group_by(Patient) %>%
  summarise(avg = mean(myv)) %>%
  select(Patient, avg) %>%
  data.frame() %>% 
  arrange(-avg)

# d) Plot data -----------------------------------------------------------------

p.data <- data %>%
  merge(myavg, by= "Patient") %>%
  mutate(Patient = factor(Patient, levels = myavg$Patient)) %>%
  arrange(Patient)

# e) Dot plot ------------------------------------------------------------------

p <- ggplot() + 
  geom_point(data = p.data, aes(x=Patient, y=myv),
             size = 3, color = "#3288bd", alpha = 0.6) + 
  geom_line(data = p.data, aes(x=Patient,  y=avg), size = 1, group = 1,
            linetype=1, col="#08306b") +     
  theme_bw(base_size = 12, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Patient ID",y="RNA-ITH",
       title="Patient-Specific RNA-ITH (CPRIT-MIRA)") +
  theme(legend.position="none",
        aspect.ratio=1/2, 
        panel.border = element_rect(size = 1, fill = NA),
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "serif", size=12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",family = "serif", size=10,
                                   angle = 60,  vjust = 1, hjust=1),
        axis.text.y = element_text(face="bold", color="black", family = "serif",size=10),
        axis.title.x = element_text(face="bold", color="black",family = "serif", size=12),
        axis.title.y = element_text(face="bold",color="black",family = "serif", size=12)) +
  scale_y_continuous(limits=c(-0.05,1.05), breaks=seq(0,1, 0.25),expand = c(0,0)) 

pdf(myFig1, width= 10, height= 5)
print(p)
dev.off()




