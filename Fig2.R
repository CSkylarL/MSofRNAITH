#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.2
# Author: Chao Cheng; Chenyang Li
# 12/13/2022
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [1] TRACERx ##################################################################

# [1.1] Correlation between intra and inter heterogeneity at the gene level ====
rm(list=ls())
myinf1 <- "~/Mydata/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/RNAseq_Sample_FullID.txt"

outdir <- "./Fig2/"
myoutf1 <- paste0(outdir,"/Fig2_Intra-inter_variation_genes_TRACERx.rda")

data <- read.table(myinf1, sep="\t", header=T, row.names=1,  check.names=F)

info <- read.table(myinf2, sep="\t", header=T, 
                   stringsAsFactors=F, quote="", comment.char="")
# a) map ID --------------------------------------------------------------------
map <- info[,1]
names(map) <- info[,2]
xx <- substr(colnames(data), 1, 15)
tmp <- map[xx]
colnames(data) <- tmp

# b) remove low expression -----------------------------------------------------
xx <- apply(data>0, 1, sum)
se <- which(xx>=20)
data <- data[se,]
dim(data) ## 25936   161

# c) remove same expression ----------------------------------------------------
data <- log10(data+1)
xx <- apply(data, 1, sd)
se <- which(xx>0)
data <- data[se,]
dim(data)				##  25936   161

# d) calculate the intra/inter tumor variation for each gene -------------------
mytag <- substr(colnames(data),1 , 8)
pid <- unique(mytag)

df <- ncol(data)

myavg <- apply(data,1, mean)
tot.var <- apply((data-myavg)^2, 1, sum)
in.var <- rep(0, nrow(data))

for(i in 1:length(pid)){
  se <- which(mytag==pid[i])
  xx <- data[, se, drop=F]
  tmp <- apply(xx,1, mean)
  in.var <- in.var + apply((xx-tmp)^2, 1, sum)
  for(k in 1:ncol(xx)){
    xx[,k] <- tmp
  }
  if(i==1){
    uni.tmp <- xx
  }else{
    uni.tmp <- cbind(uni.tmp, xx)
  }
}
tmp <- apply(uni.tmp, 1, mean)
out.var <- apply((uni.tmp-tmp)^2, 1, sum)

res <- data.frame(tot.var, in.var, out.var)


in.frac <- round(res[,2]/res[,1],3)
res <- cbind(res, in.frac)

norm.in= res[,2]/df
norm.out <- res[,3]/df
res <- cbind(res, norm.in, norm.out)

mydata <- res
save(mydata, file = myoutf1)


# e) Scatter pllot -------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggpubr)
outdir <- "./Fig2/"
myinf1 <- paste0(outdir,"/Fig2_Intra-inter_variation_genes_TRACERx.rda")
myFig1 = paste0(outdir,"/ScatterPlot_Intra-Inter_gene_TRACERx.pdf")

load(myinf1)
ls()

p.data <- mydata
lm(p.data$norm.out ~ p.data$norm.in)

p1 <- ggplot(p.data,mapping = aes(x = norm.in,y=norm.out)) +
  geom_point(color = "#d73027" , size=2,alpha=0.3) +
  geom_smooth( method = "lm", color = "#67000d",lty = 4,size = 1,se = T,
               formula = y ~ x, alpha = 0.2) +
  stat_cor(method = "pearson", color = "#67000d",
           label.x = 0.02, label.y = 1.2, family = "serif", size = 5,
           p.accuracy = 0.001, r.accuracy = 0.001) +
  theme_bw(base_size = 15, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="Intratumor Variations",y="Intertumor Variations",
       title="Gene Expression: Pearson Correlation \nBetween Intra and Inter tumor variations") +
  theme(legend.position="top",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18)) +
  scale_x_continuous(limits=c(-0.005,0.25), 
                     breaks=seq(0,0.25, 0.05),expand = c(0,0)) +
  scale_y_continuous(limits=c(-0.02,1.3), 
                     breaks=seq(0,1.5, 0.25),expand = c(0,0)) 
  


pdf(myFig1, width= 6, height= 5)
print(p1)
dev.off()


# [1.2] Differential gene expression with extreme value ========================

rm(list=ls())

library(survival)

myinf1 <- "~/Mydata/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

outdir <- "./Fig2/"
myoutf1 <- paste0(outdir,"/Fig2_gene_coxph_TRACERx.rda")

# a) Expression and survival info ----------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1)

mytag <- substr(colnames(data), 17, 25)

info <- read.csv(myinf2, header=T, row.names=1)

sum(mytag%in%row.names(info))
sum(!mytag%in%row.names(info))

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

# b) convert gene expression data into patient matrix --------------------------
p.tag <- mytag
pid <- unique(p.tag)

# remove low expression 
xx <- apply(data>0, 1, sum)
se <- which(xx>=20)
data <- data[se,]
dim(data)

# remove same expression
data <- log10(data+1)
xx <- apply(data, 1, sd)
se <- which(xx>0)
data <- data[se,]

# normalize
data <- (data-apply(data, 1, mean))/apply(data, 1, sd) 

# mean/max/min of gene expression
tmp <- matrix(0, nrow(data), length(pid))
row.names(tmp) <- row.names(data)
colnames(tmp) <- pid
max.v <- min.v <- avg.v  <- tmp

for(k in 1:length(pid))
{
  cat("\r", k)
  se <- which(p.tag==pid[k])
  tmp <- data[, se, drop=F]
  max.v[,k]  <- apply(tmp, 1, max)		
  min.v[,k]  <- apply(tmp, 1, min)			
  avg.v[,k]  <- apply(tmp, 1, mean)			
}

# c) filter info to be consistent with gene expression -------------------------
comxx <- intersect(colnames(max.v), row.names(info))
info <- info[comxx,]
max.v <- max.v[,comxx]
min.v <- min.v[,comxx]
avg.v <- avg.v[,comxx]

# d) survival analysis univariate (gene) ---------------------------------------

tmp <- matrix(0, nrow(data), 3)
row.names(tmp) <- row.names(data)
colnames(tmp) <- c("max", "min", "avg")
tmp <- as.data.frame(tmp)
pval <- hr <- ci <- tmp


for(k in 1:nrow(max.v)){
  cat("\r", k)
  mytf1 <- as.numeric(max.v[k,])
  mytf2 <- as.numeric(min.v[k,])
  mytf3 <- as.numeric(avg.v[k,])
  xx <- cbind(mytf1, mytf2, mytf3, info)
  xx <- xx[xx[, "t.surv"]>0,]
  # max
  mycox <- coxph(Surv(t.surv, e.surv)~mytf1, xx) 
  mycox <- summary(mycox)
  pval[k, 1] <- mycox$coefficients[5]
  tmp <- mycox$conf.int
  hr[k,1] <- tmp[1]
  ci[k,1] <- mycox$concordance[1]
  # min
  mycox <- coxph(Surv(t.surv, e.surv)~mytf2, xx) 
  mycox <- summary(mycox)
  pval[k, 2] <- mycox$coefficients[5]
  tmp <- mycox$conf.int
  hr[k,2] <- tmp[1]
  ci[k,2] <- mycox$concordance[1]
  # mean
  mycox <- coxph(Surv(t.surv, e.surv)~mytf3, xx) 
  mycox <- summary(mycox)
  pval[k, 3] <- mycox$coefficients[5]
  tmp <- mycox$conf.int
  hr[k,3] <- tmp[1]
  ci[k,3] <- mycox$concordance[1]	
}

res <- data.frame(hr,pval, ci)
colnames(res) <- c("HR.max", "HR.min", "HR.avg", "PV.max", "PV.min", 
                   "PV.avg", "CI.max", "CI.min", "CI.avg")
uni.cox.rs <- res

# e) survival analysis multivariate (gene+clincal factors) ---------------------
# stage
xx <- tolower(info$Stage)
se <- grep("iii", xx)
xx[se] <- "III"
se <- grep("ii", xx)
xx[se] <- "II"
se <- grep("i", xx)
xx[se] <- "I"
table(xx)
info$Stage <- xx

# subtype
xx <- info$Subtype
table(xx)
xx <- ifelse(xx=="adenocarcinoma", 1, 2)
info$Subtype <- xx

# result table
tmp <- matrix(0, nrow(data), 3)
row.names(tmp) <- row.names(data)
colnames(tmp) <- c("max", "min", "avg")
tmp <- as.data.frame(tmp)
pval <- hr <- ci <- tmp

# make factor
tmp <- info[c("t.surv", "e.surv","Age","Stage","Subtype")]
tmp$Gender <- as.factor(info$Gender..M.F. )
tmp$Stage <- as.factor(info$Stage)
tmp$Subtype <- as.factor(info$Subtype)
info <- tmp
str(info)

# cox regression
for(k in 1:nrow(max.v)){
  cat("\r", k)
  mytf1 <- as.numeric(max.v[k,])
  mytf2 <- as.numeric(min.v[k,])
  mytf3 <- as.numeric(avg.v[k,])
  xx <- cbind(mytf1, mytf2, mytf3, info)
  xx <- xx[xx[, "t.surv"]>0,]
  # max
  mycox <- coxph(Surv(t.surv, e.surv)~mytf1+Age+Gender+Stage+Subtype, xx) 
  mycox <- summary(mycox)
  pval[k, 1] <- mycox$coefficients[1, 5]
  tmp <- mycox$conf.int
  hr[k,1] <- tmp[1, 1]
  ci[k,1] <- mycox$concordance[1]
  #min
  mycox <- coxph(Surv(t.surv, e.surv)~mytf2+Age+Gender+Stage+Subtype, xx) 
  mycox <- summary(mycox)
  pval[k, 2] <- mycox$coefficients[1, 5]
  tmp <- mycox$conf.int
  hr[k,2] <- tmp[1, 1]
  ci[k,2] <- mycox$concordance[1]
  #mean
  mycox <- coxph(Surv(t.surv, e.surv)~mytf3+Age+Gender+Stage+Subtype, xx) 
  mycox <- summary(mycox)
  pval[k, 3] <- mycox$coefficients[1, 5]
  tmp <- mycox$conf.int
  hr[k,3] <- tmp[1,1]
  ci[k,3] <- mycox$concordance[1]	
}

res <- data.frame(hr,pval, ci)
colnames(res) <- c("HR.max", "HR.min", "HR.avg", "PV.max", "PV.min", 
                   "PV.avg", "CI.max", "CI.min", "CI.avg")
mul.cox.rs <- res

save(uni.cox.rs, mul.cox.rs, file= myoutf1)


# [1.3] Volcano plot ===========================================================
rm(list = ls())
library(ggplot2)
library(gridExtra)
outdir <- "./Fig2/"
myinf1 <- paste0(outdir,"/Fig2_gene_coxph_TRACERx.rda")
myFig1<- paste0(outdir,"/Volcano_unicox_TRACERx.pdf")

load(myinf1)
ls()

data <- uni.cox.rs


# a) label Good Poor NS group --------------------------------------------------
for (i in c("avg","max", "min")){
  nam1 <- paste0("Group.",i)
  nam2 <- paste0("PV.",i)
  nam3 <- paste0("HR.",i)
  data[nam1]  <- as.factor(
    ifelse(
      data[nam2] < 0.01,
      ifelse(data[nam3] < 1, "Good", "Poor"),
      "NS"
    )
  )
}
data$label <- rownames(data)

# set the threhold for volcano plot

data <- data %>%
  mutate(HR.max = replace(HR.max,log2(HR.max) > 5, 2^5)) %>%
  mutate(HR.max = replace(HR.max,log2(HR.max) < -5, 2^-5)) %>%
  mutate(PV.max = replace(PV.max,-log10(PV.max) >= 6, 10^-6))  %>%
  mutate(HR.avg = replace(HR.avg,log2(HR.avg) > 5, 2^5)) %>%
  mutate(HR.avg = replace(HR.avg,log2(HR.avg) < -5, 2^-5)) %>%
  mutate(PV.avg = replace(PV.avg,-log10(PV.avg) >= 6, 10^-6))  %>%
  mutate(HR.min = replace(HR.min,log2(HR.min) > 5, 2^5)) %>%
  mutate(HR.min = replace(HR.min,log2(HR.min) < -5, 2^-5)) %>%
  mutate(PV.min = replace(PV.min,-log10(PV.min) >= 6, 10^-6))

# b) mean ----------------------------------------------------------------------
t <- table(data$Group.avg)
t
p1 <- ggplot() + 
  geom_point(data = data,aes(x=log2(HR.avg), y=-log10(PV.avg), 
                                  color=Group.avg, alpha=Group.avg), size=2) +
  scale_color_manual( guide = "none",
                     limits=c( "Good","Poor", "NS"),
                     values=c("#cb181d","#08519c", "#bdbdbd")) +
  scale_alpha_manual(guide = "none",
                     limits=c( "Good","Poor", "NS"),
                     values=c(0.6,0.6,0.2))+
  geom_vline(xintercept=log2(1), lty=4, col="#bdbdbd", lwd=1) + 
  geom_hline(yintercept=-log10(0.01), lty=4, col="#bdbdbd", lwd=1)+
  annotate("text", x = -3, y = 5.5, size = 5, color = "#cb181d",
           label = paste0("Protective\n",t[1]), family = "serif") +
  annotate("text", x = 3, y = 5.5, size = 5, color = "#08519c",
           label = paste0("Hazardous\n",t[3]), family = "serif") +
  theme_bw(base_size = 15, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="log2 (Hazard Ratio)",y="-log10 (P-Value)",
       title="Average Expression") +
  theme(legend.position="top",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "horizontal",
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=9),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=9),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))  +
  scale_y_continuous(limits=c(0, 6), breaks=seq(0,6, 1),expand = c(0,0)) +
  scale_x_continuous(limits=c(-5, 5), breaks=seq(-4, 4,2),expand = c(0,0))
p1

# c) max -----------------------------------------------------------------------
t <- table(data$Group.max)
t
p2 <- ggplot() + 
  geom_point(data = data,aes(x=log2(HR.max), y=-log10(PV.max), 
                             color=Group.max, alpha=Group.max), size=2) +
  scale_color_manual( guide = "none",
                      limits=c( "Good","Poor", "NS"),
                      values=c("#cb181d","#08519c", "#bdbdbd")) +
  scale_alpha_manual(guide = "none",
                     limits=c( "Good","Poor", "NS"),
                     values=c(0.6,0.6,0.2))+
  geom_vline(xintercept=log2(1), lty=4, col="#bdbdbd", lwd=1) + 
  geom_hline(yintercept=-log10(0.01), lty=4, col="#bdbdbd", lwd=1)+
  annotate("text", x = -3, y = 5.5, size = 5, color = "#cb181d",
           label = paste0("Protective\n",t[1]), family = "serif") +
  annotate("text", x = 3, y = 5.5, size = 5, color = "#08519c",
           label = paste0("Hazardous\n",t[3]), family = "serif") +
  theme_bw(base_size = 15, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="log2 (Hazard Ratio)",y="-log10 (P-Value)",
       title="Maximal Expression") +
  theme(legend.position="top",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "horizontal",
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=9),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=9),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))  +
  scale_y_continuous(limits=c(0, 6), breaks=seq(0,6, 1),expand = c(0,0)) +
  scale_x_continuous(limits=c(-5, 5), breaks=seq(-4, 4,2),expand = c(0,0))
p2
# d) min ----------------------------------------------------------------------
t <- table(data$Group.min)
t
p3 <- ggplot() + 
  geom_point(data = data,aes(x=log2(HR.min), y=-log10(PV.min), 
                             color=Group.min, alpha=Group.min), size=2) +
  scale_color_manual( guide = "none",
                     limits=c( "Good","Poor", "NS"),
                     values=c("#cb181d","#08519c", "#bdbdbd")) +
  scale_alpha_manual(guide = "none",
                     limits=c( "Good","Poor", "NS"),
                     values=c(0.6,0.6,0.2))+
  geom_vline(xintercept=log2(1), lty=4, col="#bdbdbd", lwd=1) + 
  geom_hline(yintercept=-log10(0.01), lty=4, col="#bdbdbd", lwd=1)+
  annotate("text", x = -3, y = 5.5, size = 5, color = "#cb181d",
           label = paste0("Protective\n",t[1]), family = "serif") +
  annotate("text", x = 3, y = 5.5, size = 5, color = "#08519c",
           label = paste0("Hazardous\n",t[3]), family = "serif") +
  theme_bw(base_size = 15, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="log2 (Hazard Ratio)",y="-log10 (P-Value)",
       title="Minimal Expression") +
  theme(legend.position="top",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "horizontal",
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=9),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=9),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))  +
  scale_y_continuous(limits=c(0, 6), breaks=seq(0,6, 1),
                     expand = c(0,0)) +
  scale_x_continuous(limits=c(-5, 5), breaks=seq(-4, 4,2),
                     expand = c(0,0))
p3
# 7 NA
# Warning message:
#   Removed 7 rows containing missing values (geom_point). 
# export -----------------------------------------------------------------------
pdf(myFig1, width= 18, height= 5)
grid.arrange(p1, p2,p3, ncol=3)
dev.off()



# [1.4] Venn ===================================================================
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(ggvenn)
library(gridExtra)

outdir <- "./Fig2/"
myinf1 <- paste0(outdir,"/Fig2_gene_coxph_TRACERx.rda")
myFig1<- paste0(outdir,"/Venn_unicox_TRACERx.pdf")

load(myinf1)
ls()

data <- uni.cox.rs
# a) label Good Poor NS group --------------------------------------------------
listInput.good <- list()
listInput.poor<- list()

for (i in c("avg","max", "min")){
  nam1 <- paste0("Group.",i)
  nam2 <- paste0("PV.",i)
  nam3 <- paste0("HR.",i)
  data[nam1]  <- as.factor(
    ifelse(
      data[nam2] < 0.01,
      ifelse(data[nam3] < 1, "Good", "Poor"),
      "NS"
    )
  )
  se1 <- which(data[nam1] == "Good")
  listInput.good[[i]] <-  rownames(data)[se1]
  se2 <- which(data[nam1] == "Poor")
  listInput.poor[[i]] <- rownames(data)[se2]
}
data$label <- rownames(data)

# b) Venn ----------------------------------------------------------------------
names(listInput.good) <- c("Avg.", "Max.", "Min.")
names(listInput.poor) <- c("Avg.", "Max.", "Min.")
col <- brewer.pal(8,"Pastel2")[1:3]
p1 <- ggvenn(listInput.good, names(listInput.good), show_percentage = F,
             text_size = 5,
             fill_color =col)+
  labs(title = "Protective genes") +
  theme(plot.title = element_text(size = 18, hjust = 0.5))

p2 <- ggvenn(listInput.poor, names(listInput.good), show_percentage = F,
       text_size = 5,
       fill_color = col) +
  labs(title = "Hazardous genes") +
  theme(plot.title = element_text(size = 18, hjust = 0.5))


# b) export---------------------------------------------------------------------

pdf(myFig1, width = 8, height = 5)
grid.arrange(p1, p2, ncol=2)
dev.off()

# [1.5] Scatter plot ===================================================
rm(list = ls())
library(ggplot2)

outdir <- "./Fig2/"
myinf1 <- paste0(outdir,"/Fig2_gene_coxph_TRACERx.rda")
myFig1<- paste0(outdir,"/CI_Scatter_unicox_TRACERx.pdf")
myFig2 <- paste0(outdir,"/Pvalue_Scatter_unicox_TRACERx.pdf")

load(myinf1)
ls()

data <- uni.cox.rs


# a) label Good Poor NS group --------------------------------------------------

data$label <- rownames(data)
for (i in c("avg","max", "min")){
  nam1 <- paste0("Group.",i)
  nam2 <- paste0("PV.",i)
  nam3 <- paste0("HR.",i)
  data[nam1]  <- as.factor(
    ifelse(
      data[nam2] < 0.01,
      ifelse(data[nam3] < 1, "Good", "Poor"),
      "NS"
    )
  )
}
# b) C-index Hazardous Genes ---------------------------------------------------
nn <- data %>% filter(Group.avg =="Poor") %>% nrow()
x1 <- data %>% filter(Group.avg =="Poor") %>% 
  filter(CI.max > CI.avg ) %>% nrow()
x2 <- data %>% filter(Group.avg =="Poor") %>% 
  filter(CI.max < CI.avg ) %>% nrow()

p1 <- data %>% filter(Group.avg =="Poor") %>%
  ggplot() +
  geom_point(mapping = aes(x =  CI.avg ,y=CI.max),
             color = "#08519c",size=2,alpha=0.6) +
  geom_abline(slope = 1,intercept = 0,
              lty=4,lwd=1,color="grey")+
  annotate("text", x = 0.55, y = 0.75, size = 5, color = "black",
           label = paste0(round(x1/nn*100, 2), " %"), family = "serif") +
  annotate("text", x = 0.75, y = 0.55, size = 5, color = "black",
           label = paste0(round(x2/nn*100, 2), " %"), family = "serif") +
  theme_bw(base_size = 18, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="C-Index (Avg. Expr.)",y="C-Index (Max. Expr.)",
       title="Hazardous Genes") +
  theme(legend.position="none",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "vertical",
        legend.key = element_blank(),
        legend.background=element_blank(),
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))+
  scale_y_continuous(limits=c(0.48, 0.82), 
                     breaks=seq(0.5,1, 0.1),expand = c(0,0)) +
  scale_x_continuous(limits=c(0.48, 0.82), 
                     breaks=seq(0.5,1, 0.1),expand = c(0,0))


# c) C-index Protective Genes --------------------------------------------------
nn <- data %>% filter(Group.avg =="Good") %>% nrow()
x1 <- data %>% filter(Group.avg =="Good") %>% 
  filter(CI.min > CI.avg ) %>% nrow()
x2 <- data %>% filter(Group.avg =="Good") %>% 
  filter(CI.min < CI.avg ) %>% nrow()

p2 <- data %>% filter(Group.avg =="Good") %>%
  ggplot() +
  geom_point(mapping = aes(x =  CI.avg ,y=CI.min),
             color = "#cb181d",size=2,alpha=0.6) +
  geom_abline(slope = 1,intercept = 0,
              lty=4,lwd=1,color="grey")+
  annotate("text", x = 0.55, y = 0.75, size = 5, color = "black",
           label = paste0(round(x1/nn*100, 2), " %"), family = "serif") +
  annotate("text", x = 0.75, y = 0.55, size = 5, color = "black",
           label = paste0(round(x2/nn*100, 2), " %"), family = "serif") +
  theme_bw(base_size = 18, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="C-Index (Avg. Expr.)",y="C-Index (Min. Expr.)",
       title="Protective Genes") +
  theme(legend.position="none",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "vertical",
        legend.key = element_blank(),
        legend.background=element_blank(),
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))+
  scale_y_continuous(limits=c(0.48, 0.82), 
                     breaks=seq(0.5,1, 0.1),expand = c(0,0)) +
  scale_x_continuous(limits=c(0.48, 0.82), 
                     breaks=seq(0.5,1, 0.1),expand = c(0,0))
# d) P-Value Hazardous Genes ---------------------------------------------------
nn <- data %>% filter(Group.avg =="Poor") %>% nrow()
x1 <- data %>% filter(Group.avg =="Poor") %>% 
  filter(PV.max < PV.avg ) %>% nrow()
x2 <- data %>% filter(Group.avg =="Poor") %>% 
  filter(PV.max > PV.avg ) %>% nrow()

p3 <- data %>% filter(Group.avg =="Poor") %>%
  ggplot() +
  geom_point(mapping = aes(x =  -log10(PV.avg) ,y= -log10(PV.max)),
             color = "#08519c",size=2,alpha=0.6) +
  geom_abline(slope = 1,intercept = 0,
              lty=4,lwd=1,color="grey")+
  annotate("text", x = 1.5, y = 5, size = 5, color = "black",
           label = paste0(round(x1/nn*100, 2), " %"), family = "serif") +
  annotate("text", x = 5, y = 1, size = 5, color = "black",
           label = paste0(round(x2/nn*100, 2), " %"), family = "serif") +
  theme_bw(base_size = 18, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="-log10 (P-Value) (Avg. Expr.)",y="-log10 (P-Value) (Max. Expr.)",
       title="Hazardous Genes") +
  theme(legend.position="none",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "vertical",
        legend.key = element_blank(),
        legend.background=element_blank(),
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))+
  scale_y_continuous(limits=c(0, 6), 
                     breaks=seq(0,6, 1),expand = c(0,0)) +
  scale_x_continuous(limits=c(0, 6), 
                     breaks=seq(0,6, 1),expand = c(0,0))


# e) P-Value Protective Genes --------------------------------------------------
nn <- data %>% filter(Group.avg =="Good") %>% nrow()
x1 <- data %>% filter(Group.avg =="Good") %>% 
  filter(PV.min < PV.avg ) %>% nrow()
x2 <- data %>% filter(Group.avg =="Good") %>% 
  filter(PV.min > PV.avg ) %>% nrow()

p4 <- data %>% filter(Group.avg =="Good") %>%
  ggplot() +
  geom_point(mapping = aes(x =  -log10(PV.avg) ,y= -log10(PV.min)),
             color = "#cb181d",size=2,alpha=0.6) +
  geom_abline(slope = 1,intercept = 0,
              lty=4,lwd=1,color="grey")+
  annotate("text", x = 1.5, y = 5, size = 5, color = "black",
           label = paste0(round(x1/nn*100, 2), " %"), family = "serif") +
  annotate("text", x = 5, y = 1, size = 5, color = "black",
           label = paste0(round(x2/nn*100, 2), " %"), family = "serif") +
  theme_bw(base_size = 18, base_family = "serif") + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  labs(x="-log10 (P-Value) (Avg. Expr.)",y="-log10 (P-Value) (Min. Expr.)",
       title="Protective Genes") +
  theme(legend.position="none",
        aspect.ratio = 1/1,
        panel.border = element_rect(size = 1, fill = NA),
        legend.direction = "vertical",
        legend.key = element_blank(),
        legend.background=element_blank(),
        panel.grid=element_blank(),
        legend.title = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        legend.text= element_text(face="bold", color="black",
                                  family = "serif", size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black",
                                   family = "serif", size=15),
        axis.text.y = element_text(face="bold", color="black", 
                                   family = "serif",size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    family = "serif", size=18),
        axis.title.y = element_text(face="bold",color="black",
                                    family = "serif", size=18))+
  scale_y_continuous(limits=c(0, 6), 
                     breaks=seq(0,6, 1),expand = c(0,0)) +
  scale_x_continuous(limits=c(0, 6), 
                     breaks=seq(0,6, 1),expand = c(0,0))

# f) export---------------------------------------------------------------------

pdf(myFig1, width = 8, height = 5)
grid.arrange(p1, p2, ncol=2)
dev.off()

pdf(myFig2, width = 8, height = 5)
grid.arrange(p3, p4, ncol=2)
dev.off()

# g) example -------------------------------------------------------------------
data <- uni.cox.rs
mygen <- c("DDX5", "MDS2",# fusion/oncogene
           "FCAR", "TEAD4",# other hazardous
           "KAT6B","CEBPA") # Tumor suppressor

data[mygen,]
# HR.max    HR.min    HR.avg       PV.max      PV.min     PV.avg    CI.max    CI.min    CI.avg
# DDX5  2.1984090 1.2753137 1.7347619 0.0080173941 0.300514505 0.05842963 0.6754386 0.5263158 0.6085526
# MDS2  1.5202425 0.5593640 1.2692290 0.0075816195 0.218073941 0.32476070 0.6474781 0.5729167 0.5904605
# FCAR  2.4091648 1.0955060 1.8508726 0.0003006784 0.714043608 0.02061049 0.7280702 0.5290570 0.6567982
# TEAD4 1.5568360 1.3804516 1.5613625 0.0307685092 0.200304119 0.06786114 0.6650219 0.5986842 0.6282895
# KAT6B 0.7544117 0.4725118 0.4319028 0.2629006070 0.001133486 0.01189491 0.5893640 0.7445175 0.6776316
# CEBPA 0.6703537 0.2845609 0.4463638 0.0990378357 0.003712520 0.01246308 0.6025219 0.7209430 0.6743421

# [1.6] Example ================================================================
rm(list=ls())


library(survival)
library(survminer)

myinf1 <- "~/Mydata/CohortD_TRACERx/data/RNAseq_Genes_rsem_Symbol_FPKM.txt"
myinf2 <- "~/Mydata/CohortD_TRACERx/Clinical_Sample_info/CohortD_TRACERx_clincal_info.csv"

outdir <- "./Fig2/"
myFig1 <- paste0(outdir,"/KM_example_TRACERx.pdf")

# a) Expression and survival info ----------------------------------------------
data <- read.table(myinf1, sep="\t", header=T, row.names=1)


mytag <- substr(colnames(data), 17, 25)

info <- read.csv(myinf2, header=T, row.names=1)

sum(mytag%in%row.names(info))
sum(!mytag%in%row.names(info))

t.surv <- info$ORIGINAL.publsihed.Time_to_recurrence_or_death..months.
e.surv <- info$Recurrence.or.death
info <- cbind(t.surv, e.surv, info)

# b) convert gene expression data into patient matrix --------------------------
p.tag <- mytag
pid <- unique(p.tag)

# remove low expression 
xx <- apply(data>0, 1, sum)
se <- which(xx>=20)
data <- data[se,]
dim(data)

# remove same expression
data <- log10(data+1)
xx <- apply(data, 1, sd)
se <- which(xx>0)
data <- data[se,]

# normalize
data <- (data-apply(data, 1, mean))/apply(data, 1, sd) 

# example
mygen <- c( "MDS2",# fusion/oncogene
            "KAT6B") # Tumor suppressor
data <- data[mygen,]

# mean/max/min of gene expression
tmp <- matrix(0, nrow(data), length(pid))
row.names(tmp) <- row.names(data)
colnames(tmp) <- pid
max.v <- min.v <- avg.v  <- tmp

for(k in 1:length(pid)){
  cat("\r", k)
  se <- which(p.tag==pid[k])
  tmp <- data[, se, drop=F]
  max.v[,k]  <- apply(tmp, 1, max)		
  min.v[,k]  <- apply(tmp, 1, min)			
  avg.v[,k]  <- apply(tmp, 1, mean)			
}

# c) filter info to be consistent with gene expression -------------------------
comxx <- intersect(colnames(max.v), row.names(info))
info <- info[comxx,]
max.v <- max.v[,comxx]
min.v <- min.v[,comxx]
avg.v <- avg.v[,comxx]


# d) KM plot ------------------------------------------------------------------

all(rownames(info) == colnames(max.v))

p <- list()
count <- 0
for( k in mygen)  {
  kmData <- info[c("t.surv", "e.surv" )]
  data <- data.frame(max = max.v[k,],
                     min = min.v[k,] ,
                     avg = avg.v[k,])
  
  for(i in c("avg","max", "min")){
    kmData$mytag <- as.numeric(data[,i])
    cutoff <- median(kmData$mytag)
    # low (,median] high(median,)
    kmData$Exp <- as.factor(ifelse(test <- kmData$mytag > cutoff, 
                                   "high", "low"))
    
    
    
    ta <- table(kmData$Exp)
    km_fit <-  survfit(Surv(t.surv, e.surv) ~Exp, data = kmData)
    # log-rank test
    survdiff(Surv(t.surv, e.surv) ~Exp, data = kmData, rho = 0)
    count <- count+1
    p[[count]] <- ggsurvplot(
      km_fit,                     # survfit object with calculated statistics.
      pval = TRUE, 
      data = kmData,             # data used to fit survival curves.
      palette = c("red", "blue"),
      font.legend=18,
      legend = "top",
      legend.title="Exp",
      xlab = "Months",   # customize X axis label.
      break.time.by = 5,     # break X axis in time intervals by days.
      surv.median.line = "hv",  # add the median survival pointer.
      legend.labs = c(paste0("High (n=", ta[1], ")"),
                      paste0("Low (n=", ta[2], ")")),
      ggtheme = theme_classic2(base_size=25, base_family = "serif"),
      font_family= "serif",
      size=2,
      pval.method=T,
      pval.size=6,
      Exp.table = TRUE,
      title = paste0(k,"_",i),
      conf.int = T,
      conf.int.alpha = 0.2,
      conf.int.style="ribbon" #, "step"
    )
  }
  
  
}

p.list1 <- list()

for ( i in 1:length(p)){
  p.list1[[i]] <- p[[i]]$plot
}

# e) export --------------------------------------------------------------------

pdf(myFig1, width = 15,height = 10)
do.call("grid.arrange", c(plotlist = p.list1, ncol=3))
dev.off() 

