library(readr)
library(tidyverse)
set.seed(0)
setwd("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare")

acc_frame = data.frame(R = 2:10, silhouette_score = NA, 
                       Mclust = NA, Kmeans = NA)


find_acc = function(R = 3){
  out_singular_values <- unlist(read_csv(paste("out_singular_values_R", R, ".csv", sep = ""), 
                                         col_names = FALSE))
  
  
  out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]
  
  
  out_A_all_ordered <- read_csv(paste("out_A_R", R, ".csv", sep = ""), 
                                col_names = FALSE)[,order(out_singular_values, decreasing = T)]
  
  
  
  tmp_name = paste("PC", 1:R, sep = "")
  colnames(out_A_all_ordered) = tmp_name
  
  demographic_data = read_csv("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare/metauni.csv")
  
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
set.seed(111)
for (R in 2:10) {
  acc_frame[acc_frame$R == R, -1] = find_acc(R)
}

load("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare/silhouette_scores.RData")
 
silhouette_scores[["SRKHS-TD"]] = acc_frame[[2]]

silhouette_scores = silhouette_scores %>% pivot_longer(!rank_R, names_to = "Algorithm", values_to = "Silhouette Score")

p = ggplot(silhouette_scores, aes(x = rank_R, y = `Silhouette Score`, color = Algorithm)) + 
  geom_point() + geom_line(alpha = 0.5) + xlab("Rank")

ggsave("ECAM_comparison_s_score_all_R.pdf", p, width = 6, height = 3)



R = 3
out_singular_values <- unlist(read_csv(paste("out_singular_values_R", R, ".csv", sep = ""), 
                                       col_names = FALSE))


out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]


out_A_all_ordered <- read_csv(paste("out_A_R", R, ".csv", sep = ""), 
                              col_names = FALSE)[,order(out_singular_values, decreasing = T)]



tmp_name = paste("PC", 1:R, sep = "")
colnames(out_A_all_ordered) = tmp_name

demographic_data = read_csv("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare/metauni.csv")

demographic_data = cbind(demographic_data, out_A_all_ordered)

colnames(demographic_data)[2] = "Delivery"

colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"

colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"

library(ggpubr)

var = "Exposure to Antibiotic during Labor and Delivery of Mom"
var = "Delivery"
var = "Diet"

p1 = ggplot(demographic_data)+geom_point(aes(x = PC1, y = PC2, color = Diet)) + labs(color="Diet")
p2 = ggplot(demographic_data)+geom_point(aes(x = PC1, y = PC3, color = Diet)) + labs(color="Diet")
p3 = ggplot(demographic_data)+geom_point(aes(x = PC2, y = PC3, color = Diet)) + labs(color="Diet")

p = ggarrange(p1,p2,p3,
              ncol = 3, nrow = 1,
              widths = c(1,1,1),
              common.legend = TRUE, legend = "bottom")

ggsave("ECAM_loading_SRKHSTD_r3.pdf", p, width = 30/4, height = 3)

library(factoextra)
ps00 = fviz_silhouette(s) + ggtitle(paste("SRKHS-TD\n Silhouette score =", round(mean(s[,3]), digits = 3)))+
  scale_color_discrete(labels = c("bd", "fd")) + scale_fill_discrete(labels = c("bd", "fd"))

demographic_data0 = demographic_data[,c("Diet", tmp_name)] %>% pivot_longer(!Diet, names_to = "PC", values_to = "Loading")

p00 = ggplot(demographic_data0)+geom_boxplot(aes(x = PC, y = Loading, fill = Diet)) + labs(fill="Diet")+
  ggtitle("SRKHS-TD")

p00 = ggplot(demographic_data0)+
  geom_point(aes(x = 1:length(Loading), y = Loading, color = Diet))+
  facet_wrap( ~ PC, nrow = 4)+
  labs(fill="Diet")+
  xlab("Subject")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+
  ylab(" ")+
  ggtitle("SRKHS-TD")



p01 = ggplot(demographic_data)+geom_point(aes(x = PC3, y = PC4, color = Diet)) + labs(color="Diet")+
  ggtitle("SRKHS-TD")
load("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare/1.RData")
p = ggarrange(p22, p11 + ylab(" "), p00, p01, 
              ncol = 4, nrow = 1,
              widths = c(1,1,1,1.5),
              common.legend = TRUE, legend = "bottom")

ggsave("ECAM_comparison.pdf", p, width = 10, height = 6)


p = ggarrange(ps22 + ylim(-0.4,0.5) +
                scale_color_discrete(labels = c("bd", "fd")) + 
                scale_fill_discrete(labels = c("bd", "fd")), 
              ps11 + ylab(" ") + ylim(-0.4,0.5) +
                scale_color_discrete(labels = c("bd", "fd")) + 
                scale_fill_discrete(labels = c("bd", "fd")) , ps00+ ylab(" ") + ylim(-0.4,0.5),
              ncol = 3, nrow = 1,
              widths = c(1,1,1),
              common.legend = TRUE, legend = "bottom") + 
  scale_fill_discrete(labels = c("bd", "fd"))

ggsave("ECAM_comparison_s_score.pdf", p, width = 10, height = 4)

