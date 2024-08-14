## qqplot
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

options(stringsAsFactors=F)
options(scipen = 200)

tableDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Table/Distribution/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Figure/Distribution/'

real_sigma02 = 0.6
real_sigmaa2 = 0.2

mix_var_dist_0.25 = read.table(paste0(tableDir,"var_dist_mix_0.25_h.txt"), head =T)
n_sample = 10

#mix_var_dist_0.25
mix_var_dist_0.25$sigma02_z = (mix_var_dist_0.25$sigma02-real_sigma02)/sqrt(mix_var_dist_0.25$v02)
mix_var_dist_0.25$sigma02_z <- mix_var_dist_0.25$sigma02_z[order(mix_var_dist_0.25$n, mix_var_dist_0.25$p, mix_var_dist_0.25$sigma02_z)]
mix_var_dist_0.25$sigmaa2_z = (mix_var_dist_0.25$sigmaa2-real_sigmaa2)/sqrt(mix_var_dist_0.25$v12)
mix_var_dist_0.25$sigmaa2_z <- mix_var_dist_0.25$sigmaa2_z[order(mix_var_dist_0.25$n, mix_var_dist_0.25$p, mix_var_dist_0.25$sigmaa2_z)]

mix_var_dist_0.25$n <- formatC(mix_var_dist_0.25$n,format = "e", digits = 0)
mix_var_dist_0.25$n <- paste0("n = ",mix_var_dist_0.25$n)
mix_var_dist_0.25$n <- ordered(mix_var_dist_0.25$n, levels = c("n = 5e+02","n = 1e+03","n = 2e+03","n = 5e+03",
                                                                 "n = 1e+04","n = 2e+04","n = 5e+04","n = 1e+05",
                                                                 "n = 2e+05", "n = 5e+05"))

mix_var_dist_0.25$p <- formatC(mix_var_dist_0.25$p,format = "e", digits = 0)
mix_var_dist_0.25$p <- paste0("p = ",mix_var_dist_0.25$p)
mix_var_dist_0.25$p <- ordered(mix_var_dist_0.25$p, levels = c("p = 5e+03"))
mix_var_dist_0.25$z <- rep(qnorm(ppoints(100)),n_sample)

gg_mix_var_dist_0.25_sigma02 <- ggplot(data=mix_var_dist_0.25, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_wrap(p~n, ncol = 2, scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 QQ-plot")

png(paste0(figureDir, '/mix_var_dist_0.25_sigma02.png'), width = 400, height = 600)
gg_mix_var_dist_0.25_sigma02
dev.off()

gg_mix_var_dist_0.25_sigmaa2 <- ggplot(data=mix_var_dist_0.25, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_wrap(p~n,ncol = 2, scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 QQ-plot")

png(paste0(figureDir, '/mix_var_dist_0.25_sigmaa2.png'), width = 400, height = 600)
gg_mix_var_dist_0.25_sigmaa2
dev.off()

gg_mix_var_dist_0.25 <- ggarrange(gg_mix_var_dist_0.25_sigma02, gg_mix_var_dist_0.25_sigmaa2)
png(paste0(figureDir, '/mix_var_dist_0.25.png'), width = 800, height = 1000)
gg_mix_var_dist_0.25
dev.off()