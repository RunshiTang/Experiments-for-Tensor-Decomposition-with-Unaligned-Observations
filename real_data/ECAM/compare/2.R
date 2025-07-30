library(readr)
library(tidyverse)
set.seed(0)
setwd("real_data/ECAM/compare")

acc_frame = data.frame(R = 2:10, silhouette_score = NA, 
                       Mclust = NA, Kmeans = NA)

acc_GCP_frame = data.frame(R = 2:10, silhouette_score = NA, 
                       Mclust = NA, Kmeans = NA)

acc_GRKHS_frame = data.frame(R = 2:10, silhouette_score = NA, 
                           Mclust = NA, Kmeans = NA)

find_acc = function(R = 3){
  out_singular_values <- unlist(read_csv(paste("out_singular_values_R", R, ".csv", sep = ""), 
                                         col_names = FALSE))
  
  
  out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]
  
  
  out_A_all_ordered <- read_csv(paste("out_A_R", R, ".csv", sep = ""), 
                                col_names = FALSE)[,order(out_singular_values, decreasing = T)]
  
  
  tmp_name = paste("PC", 1:R, sep = "")
  colnames(out_A_all_ordered) = tmp_name
  
  demographic_data = read_csv("metauni.csv")
  
  demographic_data = cbind(demographic_data, out_A_all_ordered)
  
  colnames(demographic_data)[2] = "Delivery"
  
  colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"
  
  colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"
  
  library(ggpubr)
  
  var = "Diet"
  
  library(cluster)
  s = silhouette(x = as.numeric(as.factor(demographic_data[[var]])),
                 dist = dist(demographic_data[,tmp_name]))
  
  library(caret)
  k = kmeans(demographic_data[,tmp_name], 2, iter.max = 100)
  a_k = confusionMatrix(as.factor(s[,1]), as.factor(k$cluster))
  
  library(mclust)
  m = Mclust(demographic_data[,tmp_name], G = 2, modelNames = "VII")
  a_m = confusionMatrix(as.factor(s[,1]), as.factor(m$classification))
  
  return(c(mean(s[,3]), a_m$overall[1], a_k$overall[1]))
}

find_acc_GRKHS = function(R = 4){
  out_singular_values <- unlist(read_csv(paste("GRKHS/out_singular_values_R", R, ".csv", sep = ""), 
                                         col_names = FALSE))
  
  
  out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]
  
  
  out_A_all_ordered <- read_csv(paste("GRKHS/out_A_R", R, ".csv", sep = ""), 
                                col_names = FALSE)[,order(out_singular_values, decreasing = T)]
  
  
  tmp_name = paste("PC", 1:R, sep = "")
  colnames(out_A_all_ordered) = tmp_name
  
  demographic_data = read_csv("metauni.csv")
  
  demographic_data = cbind(demographic_data, out_A_all_ordered)
  
  colnames(demographic_data)[2] = "Delivery"
  
  colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"
  
  colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"
  
  library(ggpubr)

  
  var = "Exposure to Antibiotic during Labor and Delivery of Mom"
  var = "Delivery"
  var = "Diet"
  
  library(cluster)
  s = silhouette(x = as.numeric(as.factor(demographic_data[[var]])),
                 dist = dist(demographic_data[,tmp_name]))
  
  library(caret)
  k = kmeans(demographic_data[,tmp_name], 2, iter.max = 100)
  a_k = confusionMatrix(as.factor(s[,1]), as.factor(k$cluster))
  
  library(mclust)
  m = Mclust(demographic_data[,tmp_name], G = 2, modelNames = "VII")
  a_m = confusionMatrix(as.factor(s[,1]), as.factor(m$classification))
  
  return(c(mean(s[,3]), a_m$overall[1], a_k$overall[1]))
}


find_acc_GCP = function(R = 4){
  out_A_all_ordered <- read_csv(paste("GCP/GCP_R", R, ".csv", sep = ""), 
                                col_names = FALSE)
  tmp_name = paste("PC", 1:R, sep = "")
  colnames(out_A_all_ordered) = tmp_name
  
  demographic_data = read_csv("metauni.csv")
  
  demographic_data = cbind(demographic_data, out_A_all_ordered)
  
  colnames(demographic_data)[2] = "Delivery"
  
  colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"
  
  colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"
  
  library(ggpubr)
  
  var = "Diet"
  
  library(cluster)
  s = silhouette(x = as.numeric(as.factor(demographic_data[[var]])),
                 dist = dist(demographic_data[,tmp_name]))
  
  library(caret)
  k = kmeans(demographic_data[,tmp_name], 2, iter.max = 100)
  a_k = confusionMatrix(as.factor(s[,1]), as.factor(k$cluster))
  
  library(mclust)
  m = Mclust(demographic_data[,tmp_name], G = 2, modelNames = "VII")
  a_m = confusionMatrix(as.factor(s[,1]), as.factor(m$classification))
  
  return(c(mean(s[,3]), a_m$overall[1], a_k$overall[1]))
}



set.seed(111)
for (R in 2:10) {
  acc_frame[acc_frame$R == R, -1] = find_acc(R)
  acc_GCP_frame[acc_GCP_frame$R == R, -1] = find_acc_GCP(R)
  acc_GRKHS_frame[acc_GRKHS_frame$R == R, -1] = find_acc_GRKHS(R)
}

load("silhouette_scores.RData")
 
silhouette_scores[["S-RKHS-TD"]] = acc_frame[[2]]
silhouette_scores[["GCP"]] = acc_GCP_frame[[2]]
silhouette_scores[["S-GRKHS-TD"]] = acc_GRKHS_frame[[2]]

silhouette_scores = silhouette_scores %>% pivot_longer(!rank_R, names_to = "Algorithm", values_to = "Silhouette Score")

p = ggplot(silhouette_scores, aes(x = rank_R, y = `Silhouette Score`, color = Algorithm, shape = Algorithm)) + 
  geom_point() + geom_line(alpha = 0.3) + xlab("Rank")

ggsave("ECAM_comparison_s_score_all_R.pdf", p, width = 6, height = 3)



R = 3
out_singular_values <- unlist(read_csv(paste("GRKHS/out_singular_values_R", R, ".csv", sep = ""), 
                                       col_names = FALSE))


out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]


out_A_all_ordered <- read_csv(paste("GRKHS/out_A_R", R, ".csv", sep = ""), 
                              col_names = FALSE)[,order(out_singular_values, decreasing = T)]



tmp_name = paste("PC", 1:R, sep = "")
colnames(out_A_all_ordered) = tmp_name

demographic_data = read_csv("metauni.csv")

demographic_data = cbind(demographic_data, out_A_all_ordered)

colnames(demographic_data)[2] = "Delivery"

colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"

colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"

library(ggpubr)

var = "Exposure to Antibiotic during Labor and Delivery of Mom"
var = "Delivery"
var = "Diet"

p1 = ggplot(demographic_data)+geom_point(aes(x = PC1, y = PC2, color = Diet)) + labs(color="Diet") + ggtitle("S-GRKHS-TD") + xlim(c(0.025,0.063))
p2 = ggplot(demographic_data)+geom_point(aes(x = PC1, y = PC3, color = Diet)) + labs(color="Diet") + ggtitle("S-GRKHS-TD") + xlim(c(0.025,0.063))
p3 = ggplot(demographic_data)+geom_point(aes(x = PC2, y = PC3, color = Diet)) + labs(color="Diet") + ggtitle("S-GRKHS-TD")

p11 = ggarrange(p1,p2,p3,
              ncol = 3, nrow = 1,
              widths = c(1,1,1),
              common.legend = TRUE, legend = "bottom") 

#ggsave("ECAM_loading_SGRKHSTD_r3.pdf", p, width = 30/4, height = 3)


R = 3
out_singular_values <- unlist(read_csv(paste("out_singular_values_R", R, ".csv", sep = ""), 
                                       col_names = FALSE))


out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]


out_A_all_ordered <- read_csv(paste("out_A_R", R, ".csv", sep = ""), 
                              col_names = FALSE)[,order(out_singular_values, decreasing = T)]



tmp_name = paste("PC", 1:R, sep = "")
colnames(out_A_all_ordered) = tmp_name

demographic_data = read_csv("metauni.csv")

demographic_data = cbind(demographic_data, out_A_all_ordered)

colnames(demographic_data)[2] = "Delivery"

colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"

colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"

library(ggpubr)

var = "Exposure to Antibiotic during Labor and Delivery of Mom"
var = "Delivery"
var = "Diet"

p4 = ggplot(demographic_data)+geom_point(aes(x = PC1, y = PC2, color = Diet)) + labs(color="Diet")+ ggtitle("S-RKHS-TD") + xlim(c(0.085,0.215))
p5 = ggplot(demographic_data)+geom_point(aes(x = PC1, y = PC3, color = Diet)) + labs(color="Diet")+ ggtitle("S-RKHS-TD") + xlim(c(0.085,0.215))
p6 = ggplot(demographic_data)+geom_point(aes(x = PC2, y = PC3, color = Diet)) + labs(color="Diet")+ ggtitle("S-RKHS-TD")

p22 = ggarrange(p1,p2,p3, p4,p5,p6,
              ncol = 3, nrow = 2,
              widths = c(1,1,1),
              common.legend = TRUE, legend = "right")



ggsave("ECAM_loading_r3.pdf", p22, width = 31/4, height = 5.5)
