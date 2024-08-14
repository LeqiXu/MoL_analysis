## qqplot
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)

options(stringsAsFactors=F)
options(scipen = 200)

tableDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Table/Distribution/more/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Figure/Distribution/more/'

norm_var_dist_1_11 = read.table(paste0(tableDir,"var_dist_norm_0.0909090909090909_h.txt"), head =T)
norm_var_dist_0.25 = read.table(paste0(tableDir,"var_dist_norm_0.25_h.txt"), head =T)
norm_var_dist_0.2 = read.table(paste0(tableDir,"var_dist_norm_0.2_h.txt"), head =T)
norm_var_dist_0.6 = read.table(paste0(tableDir,"var_dist_norm_0.6_h.txt"), head =T)

mix_var_dist_1_11 = read.table(paste0(tableDir,"var_dist_mix_0.0909090909090909_h.txt"), head =T)
mix_var_dist_0.25 = read.table(paste0(tableDir,"var_dist_mix_0.25_h.txt"), head =T)
mix_var_dist_0.2 = read.table(paste0(tableDir,"var_dist_mix_0.2_h.txt"), head =T)
mix_var_dist_0.6 = read.table(paste0(tableDir,"var_dist_mix_0.6_h.txt"), head =T)

##norm
#norm_var_dist_0.25
norm_var_dist_0.25$sigma02_z = sqrt(norm_var_dist_0.25$n) * (norm_var_dist_0.25$sigma02-0.6)/sqrt(norm_var_dist_0.25$v02)
norm_var_dist_0.25$sigma02_z <- norm_var_dist_0.25$sigma02_z[order(norm_var_dist_0.25$n, norm_var_dist_0.25$p, norm_var_dist_0.25$sigma02_z)]
norm_var_dist_0.25$sigmaa2_z = sqrt(norm_var_dist_0.25$n) * (norm_var_dist_0.25$sigmaa2-0.2)/sqrt(norm_var_dist_0.25$v12)
norm_var_dist_0.25$sigmaa2_z <- norm_var_dist_0.25$sigmaa2_z[order(norm_var_dist_0.25$n, norm_var_dist_0.25$p, norm_var_dist_0.25$sigmaa2_z)]

norm_var_dist_0.25$n <- paste0("n = ",norm_var_dist_0.25$n)
norm_var_dist_0.25$n <- ordered(norm_var_dist_0.25$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
norm_var_dist_0.25$p <- paste0("p = ",norm_var_dist_0.25$p)
norm_var_dist_0.25$p <- ordered(norm_var_dist_0.25$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
norm_var_dist_0.25$z <- rep(qnorm(ppoints(100)),25)

gg_norm_var_dist_0.25_sigma02 <- ggplot(data=norm_var_dist_0.25, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under normal assumption with Heritability = 0.25")
png(paste0(figureDir, '/norm_var_dist_0.25_sigma02.png'), width = 800, height = 1000)
gg_norm_var_dist_0.25_sigma02
dev.off()

gg_norm_var_dist_0.25_sigmaa2 <- ggplot(data=norm_var_dist_0.25, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under normal assumption with Heritability = 0.25")
png(paste0(figureDir, '/norm_var_dist_0.25_sigmaa2.png'), width = 800, height = 1000)
gg_norm_var_dist_0.25_sigmaa2
dev.off()

#norm_var_dist_1_11
norm_var_dist_1_11$sigma02_z = sqrt(norm_var_dist_1_11$n) * (norm_var_dist_1_11$sigma02-0.6)/sqrt(norm_var_dist_1_11$v02)
norm_var_dist_1_11$sigma02_z <- norm_var_dist_1_11$sigma02_z[order(norm_var_dist_1_11$n, norm_var_dist_1_11$p, norm_var_dist_1_11$sigma02_z)]
norm_var_dist_1_11$sigmaa2_z = sqrt(norm_var_dist_1_11$n) * (norm_var_dist_1_11$sigmaa2-0.06)/sqrt(norm_var_dist_1_11$v12)
norm_var_dist_1_11$sigmaa2_z <- norm_var_dist_1_11$sigmaa2_z[order(norm_var_dist_1_11$n, norm_var_dist_1_11$p, norm_var_dist_1_11$sigmaa2_z)]

norm_var_dist_1_11$n <- paste0("n = ",norm_var_dist_1_11$n)
norm_var_dist_1_11$n <- ordered(norm_var_dist_1_11$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000", "n = 30000"))
norm_var_dist_1_11$p <- paste0("p = ",norm_var_dist_1_11$p)
norm_var_dist_1_11$p <- ordered(norm_var_dist_1_11$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
norm_var_dist_1_11$z <- rep(qnorm(ppoints(100)),25)

gg_norm_var_dist_1_11_sigma02 <- ggplot(data=norm_var_dist_1_11, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under normal assumption with Heritability = 1/11")
png(paste0(figureDir, '/norm_var_dist_1_11_sigma02.png'), width = 800, height = 1000)
gg_norm_var_dist_1_11_sigma02
dev.off()

gg_norm_var_dist_1_11_sigmaa2 <- ggplot(data=norm_var_dist_1_11, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under normal assumption with Heritability = 1/11")
png(paste0(figureDir, '/norm_var_dist_1_11_sigmaa2.png'), width = 800, height = 1000)
gg_norm_var_dist_1_11_sigmaa2
dev.off()

#norm_var_dist_0.6
norm_var_dist_0.6$sigma02_z = sqrt(norm_var_dist_0.6$n) * (norm_var_dist_0.6$sigma02-0.2)/sqrt(norm_var_dist_0.6$v02)
norm_var_dist_0.6$sigma02_z <- norm_var_dist_0.6$sigma02_z[order(norm_var_dist_0.6$n, norm_var_dist_0.6$p, norm_var_dist_0.6$sigma02_z)]
norm_var_dist_0.6$sigmaa2_z = sqrt(norm_var_dist_0.6$n) * (norm_var_dist_0.6$sigmaa2-0.3)/sqrt(norm_var_dist_0.6$v12)
norm_var_dist_0.6$sigmaa2_z <- norm_var_dist_0.6$sigmaa2_z[order(norm_var_dist_0.6$n, norm_var_dist_0.6$p, norm_var_dist_0.6$sigmaa2_z)]

norm_var_dist_0.6$n <- paste0("n = ",norm_var_dist_0.6$n)
norm_var_dist_0.6$n <- ordered(norm_var_dist_0.6$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
norm_var_dist_0.6$p <- paste0("p = ",norm_var_dist_0.6$p)
norm_var_dist_0.6$p <- ordered(norm_var_dist_0.6$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
norm_var_dist_0.6$z <- rep(qnorm(ppoints(100)),25)

gg_norm_var_dist_0.6_sigma02 <- ggplot(data=norm_var_dist_0.6, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under normal assumption with Heritability = 0.6")
png(paste0(figureDir, '/norm_var_dist_0.6_sigma02.png'), width = 800, height = 1000)
gg_norm_var_dist_0.6_sigma02
dev.off()

gg_norm_var_dist_0.6_sigmaa2 <- ggplot(data=norm_var_dist_0.6, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under normal assumption with Heritability = 0.6")
png(paste0(figureDir, '/norm_var_dist_0.6_sigmaa2.png'), width = 800, height = 1000)
gg_norm_var_dist_0.6_sigmaa2
dev.off()

#norm_var_dist_0.2
norm_var_dist_0.2$sigma02_z = sqrt(norm_var_dist_0.2$n) * (norm_var_dist_0.2$sigma02-0.2)/sqrt(norm_var_dist_0.2$v02)
norm_var_dist_0.2$sigma02_z <- norm_var_dist_0.2$sigma02_z[order(norm_var_dist_0.2$n, norm_var_dist_0.2$p, norm_var_dist_0.2$sigma02_z)]
norm_var_dist_0.2$sigmaa2_z = sqrt(norm_var_dist_0.2$n) * (norm_var_dist_0.2$sigmaa2-0.05)/sqrt(norm_var_dist_0.2$v12)
norm_var_dist_0.2$sigmaa2_z <- norm_var_dist_0.2$sigmaa2_z[order(norm_var_dist_0.2$n, norm_var_dist_0.2$p, norm_var_dist_0.2$sigmaa2_z)]

norm_var_dist_0.2$n <- paste0("n = ",norm_var_dist_0.2$n)
norm_var_dist_0.2$n <- ordered(norm_var_dist_0.2$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
norm_var_dist_0.2$p <- paste0("p = ",norm_var_dist_0.2$p)
norm_var_dist_0.2$p <- ordered(norm_var_dist_0.2$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
norm_var_dist_0.2$z <- rep(qnorm(ppoints(100)),25)

gg_norm_var_dist_0.2_sigma02 <- ggplot(data=norm_var_dist_0.2, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under normal assumption with Heritability = 0.2")
png(paste0(figureDir, '/norm_var_dist_0.2_sigma02.png'), width = 800, height = 1000)
gg_norm_var_dist_0.2_sigma02
dev.off()

gg_norm_var_dist_0.2_sigmaa2 <- ggplot(data=norm_var_dist_0.2, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under normal assumption with Heritability = 0.2")
png(paste0(figureDir, '/norm_var_dist_0.2_sigmaa2.png'), width = 800, height = 1000)
gg_norm_var_dist_0.2_sigmaa2
dev.off()

## mix
#mix_var_dist_0.25
mix_var_dist_0.25$sigma02_z = sqrt(mix_var_dist_0.25$n) * (mix_var_dist_0.25$sigma02-0.6)/sqrt(mix_var_dist_0.25$v02)
mix_var_dist_0.25$sigma02_z <- mix_var_dist_0.25$sigma02_z[order(mix_var_dist_0.25$n, mix_var_dist_0.25$p, mix_var_dist_0.25$sigma02_z)]
mix_var_dist_0.25$sigmaa2_z = sqrt(mix_var_dist_0.25$n) * (mix_var_dist_0.25$sigmaa2-0.2)/sqrt(mix_var_dist_0.25$v12)
mix_var_dist_0.25$sigmaa2_z <- mix_var_dist_0.25$sigmaa2_z[order(mix_var_dist_0.25$n, mix_var_dist_0.25$p, mix_var_dist_0.25$sigmaa2_z)]

mix_var_dist_0.25$n <- paste0("n = ",mix_var_dist_0.25$n)
mix_var_dist_0.25$n <- ordered(mix_var_dist_0.25$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
mix_var_dist_0.25$p <- paste0("p = ",mix_var_dist_0.25$p)
mix_var_dist_0.25$p <- ordered(mix_var_dist_0.25$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
mix_var_dist_0.25$z <- rep(qnorm(ppoints(100)),25)

gg_mix_var_dist_0.25_sigma02 <- ggplot(data=mix_var_dist_0.25, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under beta-binomial assumption with Heritability = 0.25")
png(paste0(figureDir, '/mix_var_dist_0.25_sigma02.png'), width = 800, height = 1000)
gg_mix_var_dist_0.25_sigma02
dev.off()

gg_mix_var_dist_0.25_sigmaa2 <- ggplot(data=mix_var_dist_0.25, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under beta-binomial assumption with Heritability = 0.25")
png(paste0(figureDir, '/mix_var_dist_0.25_sigmaa2.png'), width = 800, height = 1000)
gg_mix_var_dist_0.25_sigmaa2
dev.off()

#mix_var_dist_1_11
mix_var_dist_1_11$sigma02_z = sqrt(mix_var_dist_1_11$n) * (mix_var_dist_1_11$sigma02-0.6)/sqrt(mix_var_dist_1_11$v02)
mix_var_dist_1_11$sigma02_z <- mix_var_dist_1_11$sigma02_z[order(mix_var_dist_1_11$n, mix_var_dist_1_11$p, mix_var_dist_1_11$sigma02_z)]
mix_var_dist_1_11$sigmaa2_z = sqrt(mix_var_dist_1_11$n) * (mix_var_dist_1_11$sigmaa2-0.06)/sqrt(mix_var_dist_1_11$v12)
mix_var_dist_1_11$sigmaa2_z <- mix_var_dist_1_11$sigmaa2_z[order(mix_var_dist_1_11$n, mix_var_dist_1_11$p, mix_var_dist_1_11$sigmaa2_z)]

mix_var_dist_1_11$n <- paste0("n = ",mix_var_dist_1_11$n)
mix_var_dist_1_11$n <- ordered(mix_var_dist_1_11$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
mix_var_dist_1_11$p <- paste0("p = ",mix_var_dist_1_11$p)
mix_var_dist_1_11$p <- ordered(mix_var_dist_1_11$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
mix_var_dist_1_11$z <- rep(qnorm(ppoints(100)),25)

gg_mix_var_dist_1_11_sigma02 <- ggplot(data=mix_var_dist_1_11, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under beta-binomial assumption with Heritability = 1/11")
png(paste0(figureDir, '/mix_var_dist_1_11_sigma02.png'), width = 800, height = 1000)
gg_mix_var_dist_1_11_sigma02
dev.off()

gg_mix_var_dist_1_11_sigmaa2 <- ggplot(data=mix_var_dist_1_11, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under beta-binomial assumption with Heritability = 1/11")
png(paste0(figureDir, '/mix_var_dist_1_11_sigmaa2.png'), width = 800, height = 1000)
gg_mix_var_dist_1_11_sigmaa2
dev.off()

#mix_var_dist_0.6
mix_var_dist_0.6$sigma02_z = sqrt(mix_var_dist_0.6$n) * (mix_var_dist_0.6$sigma02-0.2)/sqrt(mix_var_dist_0.6$v02)
mix_var_dist_0.6$sigma02_z <- mix_var_dist_0.6$sigma02_z[order(mix_var_dist_0.6$n, mix_var_dist_0.6$p, mix_var_dist_0.6$sigma02_z)]
mix_var_dist_0.6$sigmaa2_z = sqrt(mix_var_dist_0.6$n) * (mix_var_dist_0.6$sigmaa2-0.3)/sqrt(mix_var_dist_0.6$v12)
mix_var_dist_0.6$sigmaa2_z <- mix_var_dist_0.6$sigmaa2_z[order(mix_var_dist_0.6$n, mix_var_dist_0.6$p, mix_var_dist_0.6$sigmaa2_z)]

mix_var_dist_0.6$n <- paste0("n = ",mix_var_dist_0.6$n)
mix_var_dist_0.6$n <- ordered(mix_var_dist_0.6$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
mix_var_dist_0.6$p <- paste0("p = ",mix_var_dist_0.6$p)
mix_var_dist_0.6$p <- ordered(mix_var_dist_0.6$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
mix_var_dist_0.6$z <- rep(qnorm(ppoints(100)),25)

gg_mix_var_dist_0.6_sigma02 <- ggplot(data=mix_var_dist_0.6, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under beta-binomial assumption with Heritability = 0.6")
png(paste0(figureDir, '/mix_var_dist_0.6_sigma02.png'), width = 800, height = 1000)
gg_mix_var_dist_0.6_sigma02
dev.off()

gg_mix_var_dist_0.6_sigmaa2 <- ggplot(data=mix_var_dist_0.6, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under beta-binomial assumption with Heritability = 0.6")
png(paste0(figureDir, '/mix_var_dist_0.6_sigmaa2.png'), width = 800, height = 1000)
gg_mix_var_dist_0.6_sigmaa2
dev.off()

#mix_var_dist_0.2
mix_var_dist_0.2$sigma02_z = sqrt(mix_var_dist_0.2$n) * (mix_var_dist_0.2$sigma02-0.2)/sqrt(mix_var_dist_0.2$v02)
mix_var_dist_0.2$sigma02_z <- mix_var_dist_0.2$sigma02_z[order(mix_var_dist_0.2$n, mix_var_dist_0.2$p, mix_var_dist_0.2$sigma02_z)]
mix_var_dist_0.2$sigmaa2_z = sqrt(mix_var_dist_0.2$n) * (mix_var_dist_0.2$sigmaa2-0.05)/sqrt(mix_var_dist_0.2$v12)
mix_var_dist_0.2$sigmaa2_z <- mix_var_dist_0.2$sigmaa2_z[order(mix_var_dist_0.2$n, mix_var_dist_0.2$p, mix_var_dist_0.2$sigmaa2_z)]

mix_var_dist_0.2$n <- paste0("n = ",mix_var_dist_0.2$n)
mix_var_dist_0.2$n <- ordered(mix_var_dist_0.2$n, levels = c("n = 10000", "n = 15000", "n = 20000", "n = 25000","n = 30000"))
mix_var_dist_0.2$p <- paste0("p = ",mix_var_dist_0.2$p)
mix_var_dist_0.2$p <- ordered(mix_var_dist_0.2$p, levels = c("p = 1000", "p = 2000", "p = 5000", "p = 10000", "p = 20000"))
mix_var_dist_0.2$z <- rep(qnorm(ppoints(100)),25)

gg_mix_var_dist_0.2_sigma02 <- ggplot(data=mix_var_dist_0.2, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under beta-binomial assumption with Heritability = 0.2")
png(paste0(figureDir, '/mix_var_dist_0.2_sigma02.png'), width = 800, height = 1000)
gg_mix_var_dist_0.2_sigma02
dev.off()

gg_mix_var_dist_0.2_sigmaa2 <- ggplot(data=mix_var_dist_0.2, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("sigmaa2 distribution under beta-binomial assumption with Heritability = 0.2")
png(paste0(figureDir, '/mix_var_dist_0.2_sigmaa2.png'), width = 800, height = 1000)
gg_mix_var_dist_0.2_sigmaa2
dev.off()