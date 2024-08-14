library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

options(stringsAsFactors=F)
options(scipen = 0)

tableDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Table/Heritability/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Figure/Heritability/'

norm_hreal_1_11 = read.table(paste0(tableDir,"norm_0.0909090909090909_heritability.txt"), head =T)

norm_hreal_1_11 <- na.omit(norm_hreal_1_11)
norm_hreal_1_11$real_h2 = 1/11
norm_hreal_1_11_mse = norm_hreal_1_11 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  select(real_h2, method, n, p, mse, bias2, variance) %>% unique()
norm_hreal_1_11_mse$mse_ratio = rep(norm_hreal_1_11_mse$mse[which(norm_hreal_1_11_mse$method == "MoL")] / 
                                      norm_hreal_1_11_mse$mse[which(norm_hreal_1_11_mse$method == "PQLseq")], 2)

gg_norm_1_11_mse = ggplot(data = norm_hreal_1_11_mse, aes(x = n,y = mse_ratio)) + 
  geom_line(alpha = 0.7, colour = "#2e6cd1") + 
  ylab("MoL mse / PQL mse") + xlab("sample size (n)") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20),
        plot.title = element_text(size = 24), 
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20)) + 
  geom_hline(aes(yintercept = 1), colour = "#d1972e") + labs(title = "Heritability = 1/11") +
  scale_x_continuous(breaks = c(2000,5000,10000,20000))

png(paste0(figureDir, '/Line_norm_mse.png'), width = 600, height = 300)
gg_norm_1_11_mse
dev.off()