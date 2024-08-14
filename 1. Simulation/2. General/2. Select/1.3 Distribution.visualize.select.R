## qqplot
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

options(stringsAsFactors=F)
options(scipen = 200)

tableDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Table/Distribution/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Figure/Distribution/'

norm_var_dist_1_11 = read.table(paste0(tableDir,"var_dist_norm_0.0909090909090909_h.txt"), head =T)
#norm_var_dist_1_11
norm_var_dist_1_11$sigma02_z = sqrt(norm_var_dist_1_11$n) * (norm_var_dist_1_11$sigma02-0.6)/sqrt(norm_var_dist_1_11$v02)
norm_var_dist_1_11$sigma02_z <- norm_var_dist_1_11$sigma02_z[order(norm_var_dist_1_11$n, norm_var_dist_1_11$p, norm_var_dist_1_11$sigma02_z)]
norm_var_dist_1_11$sigmaa2_z = sqrt(norm_var_dist_1_11$n) * (norm_var_dist_1_11$sigmaa2-0.06)/sqrt(norm_var_dist_1_11$v12)
norm_var_dist_1_11$sigmaa2_z <- norm_var_dist_1_11$sigmaa2_z[order(norm_var_dist_1_11$n, norm_var_dist_1_11$p, norm_var_dist_1_11$sigmaa2_z)]

norm_var_dist_1_11$n <- paste0("n = ",norm_var_dist_1_11$n)
norm_var_dist_1_11$n <- ordered(norm_var_dist_1_11$n, levels = c("n = 2000","n = 5000","n = 10000","n = 20000"))
norm_var_dist_1_11$p <- paste0("p = ",norm_var_dist_1_11$p)
norm_var_dist_1_11$p <- ordered(norm_var_dist_1_11$p, levels = c("p = 1000"))
norm_var_dist_1_11$z <- rep(qnorm(ppoints(50)),4)

gg_norm_var_dist_1_11_sigma02 <- ggplot(data=norm_var_dist_1_11, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under normal assumption with Heritability = 1/11")
png(paste0(figureDir, '/norm_var_dist_1_11_sigma02.png'), width = 600, height = 300)
gg_norm_var_dist_1_11_sigma02
dev.off()

gg_norm_var_dist_1_11_sigmaa2 <- ggplot(data=norm_var_dist_1_11, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under normal assumption with Heritability = 1/11")
png(paste0(figureDir, '/norm_var_dist_1_11_sigmaa2.png'), width = 600, height = 300)
gg_norm_var_dist_1_11_sigmaa2
dev.off()