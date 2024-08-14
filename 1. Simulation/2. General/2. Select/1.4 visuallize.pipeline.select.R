library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

options(stringsAsFactors=F)
options(scipen = 0)


tableDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Table/Heritability/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Figure/Heritability/'

mix_hreal_1_11 = read.table(paste0(tableDir,"mix_0.0909090909090909_heritability.txt"), head =T)
mix_hreal_0.25 = read.table(paste0(tableDir,"mix_0.25_heritability.txt"), head =T)
mix_hreal_0.2 = read.table(paste0(tableDir,"mix_0.2_heritability.txt"), head =T)
mix_hreal_0.6 = read.table(paste0(tableDir,"mix_0.6_heritability.txt"), head =T)

mix_hreal_1_11 <- na.omit(mix_hreal_1_11)
mix_hreal_0.25 <- na.omit(mix_hreal_0.25)
mix_hreal_0.2 <- na.omit(mix_hreal_0.2)
mix_hreal_0.6 <- na.omit(mix_hreal_0.6)

mix_hreal_1_11$real_h2 = 1/11
mix_hreal_0.25$real_h2 = 0.25
mix_hreal_0.2$real_h2 = 0.2
mix_hreal_0.6$real_h2 = 0.6

mix_hreal_1_11_mse = mix_hreal_1_11 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  select(real_h2, method, n, p, mse, bias2, variance) %>% unique()
mix_hreal_0.25_mse = mix_hreal_0.25 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  select(real_h2, method, n, p, mse, bias2, variance) %>% unique()
mix_hreal_0.2_mse = mix_hreal_0.2 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  select(real_h2, method, n, p, mse, bias2, variance) %>% unique()
mix_hreal_0.6_mse = mix_hreal_0.6 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  select(real_h2, method, n, p, mse, bias2, variance) %>% unique()

mix_hreal_1_11_mse = mix_hreal_1_11_mse[which(mix_hreal_1_11_mse$method == "MoL"),]
mix_hreal_1_11_mse$p <- ordered(mix_hreal_1_11_mse$p, levels=c("500", "1000", "2000","5000","10000"))

mix_hreal_0.25_mse = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "MoL"),]
mix_hreal_0.25_mse$p <- ordered(mix_hreal_0.25_mse$p, levels=c("500", "1000", "2000","5000","10000"))

mix_hreal_0.2_mse = mix_hreal_0.2_mse[which(mix_hreal_0.2_mse$method == "MoL"),]
mix_hreal_0.2_mse$p <- ordered(mix_hreal_0.2_mse$p, levels=c("500", "1000", "2000","5000","10000"))

mix_hreal_0.6_mse = mix_hreal_0.6_mse[which(mix_hreal_0.6_mse$method == "MoL"),]
mix_hreal_0.6_mse$p <- ordered(mix_hreal_0.6_mse$p, levels=c("500", "1000", "2000","5000","10000"))

plabel = unique(paste0("p = ", mix_hreal_1_11_mse$p))
names(plabel) = unique(mix_hreal_1_11_mse$p)

gg_mix_1_11_mse = ggplot(data = mix_hreal_1_11_mse, aes(x = n,y = mse)) + 
  geom_line(alpha = 0.7, colour = "#2e6cd1") + 
  ylab("MoL mse / PQL mse") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") + labs(title = "Heritability = 1/11")

gg_mix_0.25_mse = ggplot(data = mix_hreal_0.25_mse, aes(x = n,y = mse)) + 
  geom_line(alpha = 0.7, colour = "#2e6cd1") + 
  ylab("MoL mse / PQL mse") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") + labs(title = "Heritability = 0.25")

gg_mix_0.2_mse = ggplot(data = mix_hreal_0.2_mse, aes(x = n,y = mse)) + 
  geom_line(alpha = 0.7, colour = "#2e6cd1") + 
  ylab("MoL mse / PQL mse") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") + labs(title = "Heritability = 0.2")

gg_mix_0.6_mse = ggplot(data = mix_hreal_0.6_mse, aes(x = n,y = mse)) + 
  geom_line(alpha = 0.7, colour = "#2e6cd1") + 
  ylab("MoL mse / PQL mse") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") + labs(title = "Heritability = 0.6")

png(paste0(figureDir, '/Line_mix_mse.png'), width = 1600, height = 1200)
gg_mix_mse = ggarrange(gg_mix_0.25_mse, gg_mix_1_11_mse, gg_mix_0.6_mse, gg_mix_0.2_mse, ncol = 1, nrow = 4)
gg_mix_mse
dev.off()

mix_hreal_1_11$method <- ordered(mix_hreal_1_11$method, levels=c("MoL", "PQLseq"))
mix_hreal_1_11$n <- ordered(mix_hreal_1_11$n, levels=c("200", "500", "1000", "2000","5000"))
mix_hreal_1_11$p <- ordered(mix_hreal_1_11$p, levels=c("500", "1000", "2000","5000","10000"))

mix_hreal_0.25$method <- ordered(mix_hreal_0.25$method, levels=c("MoL", "PQLseq"))
mix_hreal_0.25$n <- ordered(mix_hreal_0.25$n, levels=c("200", "500", "1000", "2000","5000"))
mix_hreal_0.25$p <- ordered(mix_hreal_0.25$p, levels=c("500", "1000", "2000","5000","10000"))

mix_hreal_0.2$method <- ordered(mix_hreal_0.2$method, levels=c("MoL", "PQLseq"))
mix_hreal_0.2$n <- ordered(mix_hreal_0.2$n, levels=c("200", "500", "1000", "2000","5000"))
mix_hreal_0.2$p <- ordered(mix_hreal_0.2$p, levels=c("500", "1000", "2000","5000","10000"))

mix_hreal_0.6$method <- ordered(mix_hreal_0.6$method, levels=c("MoL", "PQLseq"))
mix_hreal_0.6$n <- ordered(mix_hreal_0.6$n, levels=c("200", "500", "1000", "2000","5000"))
mix_hreal_0.6$p <- ordered(mix_hreal_0.6$p, levels=c("500", "1000", "2000","5000","10000"))

gg_mix_1_11 = ggplot(data = mix_hreal_1_11, aes(x = n,y = heritability)) + 
  geom_boxplot(aes(fill = method), outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 1/11, linetype = "dashed", color = "red") +
  scale_fill_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e")) + 
  ylab("estimated heritability") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") +  labs(title = "Heritability = 1/11")

gg_mix_0.25 = ggplot(data = mix_hreal_0.25, aes(x = n,y = heritability)) + 
  geom_boxplot(aes(fill = method), outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "red") +
  scale_fill_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e")) + 
  ylab("estimated heritability") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") + labs(title = "Heritability = 0.25")

gg_mix_0.2 = ggplot(data = mix_hreal_0.2, aes(x = n,y = heritability)) + 
  geom_boxplot(aes(fill = method), outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
  scale_fill_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e")) + 
  ylab("estimated heritability") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") + labs(title = "Heritability = 0.2")

gg_mix_0.6 = ggplot(data = mix_hreal_0.6, aes(x = n,y = heritability)) + 
  geom_boxplot(aes(fill = method), outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "red") +
  scale_fill_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e"))+ 
  ylab("estimated heritability") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel), scales = "free") +  labs(title = "Heritability = 0.6")

png(paste0(figureDir, '/Boxplot_mix.png'), width = 1600, height = 1200)
gg_mix = ggarrange(gg_mix_0.25, gg_mix_1_11, gg_mix_0.6, gg_mix_0.2, ncol = 1, nrow = 4)
gg_mix
dev.off()

mix = rbind(mix_hreal_1_11, mix_hreal_0.25, mix_hreal_0.2, mix_hreal_0.6)
mix = na.omit(mix)
mix_time = aggregate(mix$comp_time, list(mix$real_h2, mix$method, mix$n,mix$p), FUN=mean)
colnames(mix_time) = c("real_h2", "method", "n", "p", "mean_comp_time")
mix_time$n = as.integer(mix_time$n)
mix_time$n = ifelse(mix_time$n == 1, 200, 
                    ifelse(mix_time$n == 2, 500, 
                           ifelse(mix_time$n == 3, 1000, ifelse(mix_time$n == 4, 2000, 5000))))
mix_time$n = as.factor(mix_time$n)
mix_time_1_11 = mix_time[which(mix_time$real_h2 == 1/11), ]
mix_time_0.25 = mix_time[which(mix_time$real_h2 == 0.25), ]
mix_time_0.2 = mix_time[which(mix_time$real_h2 == 0.2), ]
mix_time_0.6 = mix_time[which(mix_time$real_h2 == 0.6), ]


gg_mix_time_1_11 = ggplot(data = mix_time_1_11, mapping = aes(x=n, y=mean_comp_time, colour=method, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_colour_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e"))  + 
  ylab("mean time") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel)) +labs(title = "Heritability = 1/11")

gg_mix_time_0.25 = ggplot(data = mix_time_0.25, mapping = aes(x=n, y=mean_comp_time, colour=method, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_colour_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e")) + 
  ylab("mean time") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel)) +labs(title = "Heritability = 0.25")

gg_mix_time_0.2 = ggplot(data = mix_time_0.2, mapping = aes(x=n, y=mean_comp_time, colour=method, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_colour_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e")) + 
  ylab("mean time") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel)) +labs(title = "Heritability = 0.2")

gg_mix_time_0.6 =ggplot(data = mix_time_0.6, mapping = aes(x=n, y=mean_comp_time, colour=method, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_colour_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e"))  + 
  ylab("mean time") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel)) +labs(title = "Heritability = 0.6")

png(paste0(figureDir, '/Bar_mix_time.png'), width = 1500, height = 1200)
gg_mix_time = ggarrange(gg_mix_time_0.25, gg_mix_time_1_11, gg_mix_time_0.6, gg_mix_time_0.2, ncol = 1, nrow = 4)
gg_mix_time
dev.off()

all = rbind(norm,mix)
aggregate(all$comp_time,by = list(all$method,all$n), mean)
aggregate(all$comp_time,by = list(all$method,all$p), mean)

all_time = aggregate(all$comp_time, list(all$method, all$n,all$p), FUN=mean)
colnames(all_time) = c("method", "n", "p", "mean_comp_time")
all_time$n = as.integer(all_time$n)
all_time$n = ifelse(all_time$n == 1, 200, 
                    ifelse(all_time$n == 2, 500, 
                           ifelse(all_time$n == 3, 1000, ifelse(all_time$n == 4, 2000, 5000))))
all_time$n = as.factor(all_time$n)

plabel = unique(paste0("p = ", all_time$p))
names(plabel) = unique(all_time$p)
gg_all_time = ggplot(data = all_time, mapping = aes(x=n, y=log(mean_comp_time), colour=method, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_colour_manual(labels = c("MoL", "PQLseq"),values=c("#2e6cd1", "#d1972e")) + 
  ylab("mean time") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) +
  facet_wrap(~ p, ncol = 5, labeller = labeller(p = plabel))

png(paste0(figureDir, '/Bar_log_all_time.png'), width = 1500, height = 800)
gg_all_time
dev.off()