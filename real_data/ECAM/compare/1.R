library(readr)
library(tidyverse)

setwd("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare")

load("1.RData")

out_singular_values <- unlist(read_csv("out_singular_values_R4.csv", 
                                       col_names = FALSE))


out_singular_values_ordered = out_singular_values[order(out_singular_values, decreasing = T)]


out_A_all_ordered <- read_csv("out_A_R4.csv", 
                              col_names = FALSE)[,order(out_singular_values, decreasing = T)]

colnames(out_A_all_ordered)[1:2] = c("PC1", "PC2")
colnames(out_A_all_ordered)[3:4] = c("PC3", "PC4")


demographic_data = read_csv("//sscwin/dfsroot/users/rtang56/Desktop/ECAM_compare/metauni.csv")

demographic_data = cbind(demographic_data, out_A_all_ordered)

colnames(demographic_data)[2] = "Delivery"

colnames(demographic_data)[5] = "Exposure to Antibiotic during Labor and Delivery of Mom"

colnames(demographic_data)[colnames(demographic_data) == "diet"] = "Diet"

library(ggpubr)

var = "Exposure to Antibiotic during Labor and Delivery of Mom"
var = "Delivery"
var = "Diet"

p = ggscatter(demographic_data, x = "PC3", y = "PC4",
              color = var, palette = "jco",
              size = 3, alpha = 0.6)+
  border()
  

xp = ggboxplot(demographic_data, x = var, y = "PC3", fill = var,
               palette = "jco", 
               orientation = "horizontal")+ clean_theme() 

yp = ggboxplot(demographic_data, x = var, y = "PC4", fill = var,
               palette = "jco") + clean_theme() 

p1 = ggarrange(yp, p, NULL, xp, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(0.7, 2), heights = c(2, 1),
          common.legend = TRUE)

p1 = annotate_figure(p1,
                     top = text_grob("RKHS-TD", face = "bold"),
)

ggsave("ECAM_cluster.pdf", p1, width = 6,height = 4.5)

library(cluster)
s = silhouette(x = as.numeric(as.factor(demographic_data[[var]])),
               dist = dist(demographic_data[,c("PC1", "PC2", "PC3", "PC4")]))

mean(s[,3])


library(GGally)
ggpairs(demographic_data, columns = 7:10, ggplot2::aes(colour = diet))














demographic_data0 = demographic_data[,c(3,7:10)] %>% pivot_longer(!Diet, names_to = "PC", values_to = "Loading")

p00 = ggplot(demographic_data0)+geom_boxplot(aes(x = PC, y = Loading, fill = Diet)) + labs(fill="Diet")+
  ggtitle("SRKHS-TD")
p01 = ggplot(demographic_data)+geom_point(aes(x = PC3, y = PC4, color = Diet)) + labs(color="Diet")+
  ggtitle("SRKHS-TD")

p = ggarrange(p22, p11, p00, p01, 
               ncol = 4, nrow = 1,
               widths = c(1,1,1,1),
               common.legend = TRUE, legend = "bottom")

ggsave("ECAM_comparison.pdf", p, width = 10, height = 3)




