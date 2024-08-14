#\sigma_{0}^{2}: n=p=500 and n=p=1000;
#\sigma_{\alpha}^{2}: n=p=2000 and n=p=5000.

## qqplot
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

options(stringsAsFactors=F)
options(scipen = 200)

tableDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Table/Distribution/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Figure/Distribution/'

#norm_var_dist_0.25 = read.table(paste0(tableDir,"fix_var_dist_norm_0.25_h.txt"), head =T)
norm_var_dist_0.25 = read.table(paste0(tableDir,"var_dist_norm_0.25_h.txt"), head =T)
n_sample = 1
norm_var_dist_0.25 = norm_var_dist_0.25[order(norm_var_dist_0.25$n,norm_var_dist_0.25$p),]

real_sigma02 = 0.6
real_sigmaa2 = 0.2

#norm_var_dist_0.25
norm_var_dist_0.25_sigma02 <- norm_var_dist_0.25[which((norm_var_dist_0.25$n == 1000 & norm_var_dist_0.25$p == 1000)),]
norm_var_dist_0.25_sigma02$sigma02_z <- (norm_var_dist_0.25_sigma02$sigma02-real_sigma02)/sqrt(norm_var_dist_0.25_sigma02$v02)
norm_var_dist_0.25_sigma02$sigma02_z <- norm_var_dist_0.25_sigma02$sigma02_z[order(norm_var_dist_0.25_sigma02$n, norm_var_dist_0.25_sigma02$p, norm_var_dist_0.25_sigma02$sigma02_z)]
norm_var_dist_0.25_sigma02$sigmaa2_z = (norm_var_dist_0.25_sigma02$sigmaa2-real_sigmaa2)/sqrt(norm_var_dist_0.25_sigma02$v12)
norm_var_dist_0.25_sigma02$sigmaa2_z <- norm_var_dist_0.25_sigma02$sigmaa2_z[order(norm_var_dist_0.25_sigma02$n, norm_var_dist_0.25_sigma02$p, norm_var_dist_0.25_sigma02$sigmaa2_z)]

norm_var_dist_0.25_sigma02$n <- formatC(norm_var_dist_0.25_sigma02$n,format = "e", digits = 0)
norm_var_dist_0.25_sigma02$n <- paste0("n = ",norm_var_dist_0.25_sigma02$n)
norm_var_dist_0.25_sigma02$n <- ordered(norm_var_dist_0.25_sigma02$n, levels = c("n = 1e+03"))

norm_var_dist_0.25_sigma02$p <- formatC(norm_var_dist_0.25_sigma02$p,format = "e", digits = 0)
norm_var_dist_0.25_sigma02$p <- paste0("p = ",norm_var_dist_0.25_sigma02$p)
norm_var_dist_0.25_sigma02$p <- ordered(norm_var_dist_0.25_sigma02$p, levels = c("p = 1e+03"))
norm_var_dist_0.25_sigma02$z <- rep(qnorm(ppoints(100)),n_sample)

gg_norm_var_dist_0.25_sigma02 <- ggplot(data=norm_var_dist_0.25_sigma02, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_wrap(p~n, ncol = 2, scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 QQ-plot")

png(paste0(figureDir, '/sigma02_2.png'), width = 400, height = 400)
gg_norm_var_dist_0.25_sigma02
dev.off()

#norm_var_dist_0.25
norm_var_dist_0.25_sigmaa2 <- norm_var_dist_0.25[which((norm_var_dist_0.25$n == 20000)),]
norm_var_dist_0.25_sigmaa2$sigma02_z <- (norm_var_dist_0.25_sigmaa2$sigma02-real_sigmaa2)/sqrt(norm_var_dist_0.25_sigmaa2$v02)
norm_var_dist_0.25_sigmaa2$sigma02_z <- norm_var_dist_0.25_sigmaa2$sigma02_z[order(norm_var_dist_0.25_sigmaa2$n, norm_var_dist_0.25_sigmaa2$p, norm_var_dist_0.25_sigmaa2$sigma02_z)]
norm_var_dist_0.25_sigmaa2$sigmaa2_z = (norm_var_dist_0.25_sigmaa2$sigmaa2-real_sigmaa2)/sqrt(norm_var_dist_0.25_sigmaa2$v12)
norm_var_dist_0.25_sigmaa2$sigmaa2_z <- norm_var_dist_0.25_sigmaa2$sigmaa2_z[order(norm_var_dist_0.25_sigmaa2$n, norm_var_dist_0.25_sigmaa2$p, norm_var_dist_0.25_sigmaa2$sigmaa2_z)]

#norm_var_dist_0.25
norm_var_dist_0.25$sigma02_z = (norm_var_dist_0.25$sigma02-real_sigma02)/sqrt(norm_var_dist_0.25$v02)
norm_var_dist_0.25$sigma02_z <- norm_var_dist_0.25$sigma02_z[order(norm_var_dist_0.25$n, norm_var_dist_0.25$p, norm_var_dist_0.25$sigma02_z)]
norm_var_dist_0.25$sigmaa2_z = (norm_var_dist_0.25$sigmaa2-real_sigmaa2)/sqrt(norm_var_dist_0.25$v12)
norm_var_dist_0.25$sigmaa2_z <- norm_var_dist_0.25$sigmaa2_z[order(norm_var_dist_0.25$n, norm_var_dist_0.25$p, norm_var_dist_0.25$sigmaa2_z)]

norm_var_dist_0.25_sigmaa2$n <- formatC(norm_var_dist_0.25_sigmaa2$n,format = "e", digits = 0)
norm_var_dist_0.25_sigmaa2$n <- paste0("n = ",norm_var_dist_0.25_sigmaa2$n)
norm_var_dist_0.25_sigmaa2$n <- ordered(norm_var_dist_0.25_sigmaa2$n, levels = c("n = 2e+04"))

norm_var_dist_0.25_sigmaa2$p <- formatC(norm_var_dist_0.25_sigmaa2$p,format = "e", digits = 0)
norm_var_dist_0.25_sigmaa2$p <- paste0("p = ",norm_var_dist_0.25_sigmaa2$p)
norm_var_dist_0.25_sigmaa2$p <- ordered(norm_var_dist_0.25_sigmaa2$p, levels = c("p = 5e+03"))
norm_var_dist_0.25_sigmaa2$z <- rep(qnorm(ppoints(100)),n_sample)

gg_norm_var_dist_0.25_sigmaa2 <- ggplot(data=norm_var_dist_0.25_sigmaa2, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_wrap(p~n, ncol = 2, scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 QQ-plot")

png(paste0(figureDir, 'sigmaa2_3.png'), width = 400, height = 400)
gg_norm_var_dist_0.25_sigmaa2
dev.off()

gg_norm_var_dist_0.25 <- ggarrange(gg_norm_var_dist_0.25_sigma02, gg_norm_var_dist_0.25_sigmaa2)
png(paste0(figureDir, '/sigma02_sigmaa2.png'), width = 1600, height = 400)
gg_norm_var_dist_0.25
dev.off()