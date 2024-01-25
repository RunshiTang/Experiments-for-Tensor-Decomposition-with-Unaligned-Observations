{############# Real data analysis
rm(list=ls())
library(rTensor)
source('FTSVD.R')
library(MASS)
library(fdapace)
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(tidyverse)

setwd("C:/Users/97671/Box/HighD/RKHS/code/functional_tensor_svd-main")




##############
# 1. ECAM data
##############
counts0 <- read.csv("genus_count_cleaned.csv", row.names=1)
metadata <- read.csv("genus_metadata_cleaned.csv", row.names=1)
table(rownames(counts0)==rownames(metadata))
metauni <- unique(metadata[,c('studyid', 'delivery', 'diet', 
                              "mom_prenatal_abx",
                              "mom_ld_abx",
                              "sex")])
study.ID <- metauni$studyid
nsub <- nrow(metauni)

microbiome.data <- format_tfpca(counts0, metadata$day_of_life, metadata$studyid, threshold=0.9, pseudo_count=0.5)
# summarize the available month data...
count = rep(0,25)
p1 = nrow(metauni)
p2 = dim(microbiome.data[[1]])[1] - 1
for (i in 1:p1) {
  tn = round(microbiome.data[[i]][1,]/30)
  tn[tn==25] = 24
  tn = unique(tn)
  count[tn+1] = count[tn+1] + 1
}

# transform the data to a 42-50-25 tensor with missing values.

n = 19
microbiome.tensor = array(NA, dim= c(42,50,19))
for (i in 1:p1) {
  tn = round(microbiome.data[[i]][1,]/30)
  tn[tn==25] = 24
  for (k in 1:length(tn)) {
    if(tn[k]<=12)
      microbiome.tensor[i,,tn[k]+1] = microbiome.data[[i]][2:(p2+1),k]
    else
      microbiome.tensor[i,,ceiling(tn[k]/2)+7] = microbiome.data[[i]][2:(p2+1),k]
  }
}

data.interpolate = microbiome.tensor
for (i in 1:42) {
  for (j in 1:50) {
    if(all(is.na(microbiome.tensor[i,j,]))) data.interpolate[i,j,] = 0
    else data.interpolate[i,j,] = interpolate(microbiome.tensor[i,j,])
  }
}
data.interpolate = as.tensor(data.interpolate)

tn = c(0:12,14,16,18,20,22,24)

#res = FTSVD(data.interpolate, f.grid = c(list(NULL),list(NULL),list(tn/24)), 
#                  rank=3, penalty = 1e-5)
data.interpolate = data.interpolate - mean(data.interpolate@data)

lambda_max = 0
for (i in 1:3){
  lambda_max = max(lambda_max, svd(k_unfold(data.interpolate, i)@data)$d[1])
}
data.interpolate = data.interpolate / (lambda_max / sqrt(dim(data.interpolate)[3]))
}



silhouette_scores = data.frame(rank_R = 2:10, CP = NA, FTSVD = NA)

for (rank_R in 2:10) {
  {
    cpm = cp(data.interpolate, rank_R)
    
    all_PCs = paste("PC", 1:rank_R, sep = "")
    A.PC = cpm$U[[1]]
    colnames(A.PC) = all_PCs
    A.data <- metauni
    rownames(A.data) <- A.data$studyid
    A.data <- cbind(A.PC, A.data)
    npc <- ncol(A.PC)
    p_deliv_list <- vector("list", npc*(npc-1)/2)
    ij <- 1
    
    colnames(A.data)[colnames(A.data) == "mom_ld_abx"] = "Exposure to Antibiotic during Labor and Delivery of Mom"
    
    var = "diet"
    
    s = silhouette(x = as.numeric(as.factor(A.data[[var]])),
                   dist = dist(A.data[,all_PCs]))
    
    silhouette_scores[silhouette_scores$rank_R == rank_R, "CP"] = 
      mean(s[,3])
  }
  
  {
    # loading analysis.
    set.seed(111)
    res = FTSVD(data.interpolate, f.grid = c(list(NULL),list(NULL),list(tn/24)), 
                rank=rank_R, penalty = 0.00001)
    A.PC = res[[2]]
    
    all_PCs = paste("PC", 1:rank_R, sep = "")
    colnames(A.PC) = all_PCs
    
    A.data <- metauni
    rownames(A.data) <- A.data$studyid
    A.data <- cbind(A.PC, A.data)
    npc <- ncol(A.PC)
    p_deliv_list <- vector("list", npc*(npc-1)/2)
    ij <- 1
    
    colnames(A.data)[colnames(A.data) == "mom_ld_abx"] = "Exposure to Antibiotic during Labor and Delivery of Mom"
    colnames(A.data)[colnames(A.data) == "delivery"] = "Delivery"
    
    var = "diet"
    
    library(cluster)
    s = silhouette(x = as.numeric(as.factor(A.data[[var]])),
                   dist = dist(A.data[,all_PCs]))
    
    silhouette_scores[silhouette_scores$rank_R == rank_R, "FTSVD"] = mean(s[,3])
  }  
  
}


rm(list=setdiff(ls(), c("silhouette_scores")))
save.image("silhouette_scores.RData")





