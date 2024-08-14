library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

options(stringsAsFactors=F)
options(scipen = 0)

tableDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Table/Heritability/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Figure/Heritability/'

mix_hreal_0.25 = read.table(paste0(tableDir,"mix_0.25_heritability.txt"), head =T)

mix_hreal_0.25 <- na.omit(mix_hreal_0.25)
mix_hreal_0.25$real_h2 = 0.25
mix_hreal_0.25_mse = mix_hreal_0.25 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  mutate(time = mean(comp_time)) %>%
  select(real_h2, method, n, p, mse, bias2, variance,time) %>% unique()

# mse ratio
#mix_hreal_0.25_mse_mol = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "MoL" & ! mix_hreal_0.25_mse$n %in% c(2000, 100000)),]
mix_hreal_0.25_mse_mol = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "MoL"),]
mix_hreal_0.25_mse_pql = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "PQLseq"),]

mix_hreal_0.25_mse_ratio = rbind(mix_hreal_0.25_mse_mol,mix_hreal_0.25_mse_mol,
                                  mix_hreal_0.25_mse_mol,mix_hreal_0.25_mse_mol)[,c("n","p")]
colnames(mix_hreal_0.25_mse_ratio) = c("MoL_n","p")
mix_hreal_0.25_mse_ratio$PQL_n = rep(mix_hreal_0.25_mse_pql$n, each = length(mix_hreal_0.25_mse_mol$n))
mix_hreal_0.25_mse_ratio$var_ratio = mix_hreal_0.25_mse_mol$variance / rep(mix_hreal_0.25_mse_pql$variance, each = length(mix_hreal_0.25_mse_mol$variance))
mix_hreal_0.25_mse_ratio$mse_ratio = mix_hreal_0.25_mse_mol$mse / rep(mix_hreal_0.25_mse_pql$mse, each = length(mix_hreal_0.25_mse_mol$mse))
mix_hreal_0.25_mse_ratio$time_ratio = mix_hreal_0.25_mse_mol$time / rep(mix_hreal_0.25_mse_pql$time, each = length(mix_hreal_0.25_mse_mol$time))
mix_hreal_0.25_mse_ratio = merge(mix_hreal_0.25_mse_ratio, mix_hreal_0.25_mse_pql[,c("n", "variance")], 
                                  by.x = c("PQL_n"), by.y = c("n"))
mix_hreal_0.25_mse_ratio = merge(mix_hreal_0.25_mse_ratio, mix_hreal_0.25_mse_mol[,c("n", "variance")], 
                                  by.x = c("MoL_n"), by.y = c("n"))
names(mix_hreal_0.25_mse_ratio)[7:8] = c("PQL_variance", "MoL_variance")

mix_hreal_0.25_mse_ratio_p = mix_hreal_0.25_mse_ratio[which(mix_hreal_0.25_mse_ratio$MoL_n
                                                              >= mix_hreal_0.25_mse_ratio$PQL_n &
                                                                mix_hreal_0.25_mse_ratio$MoL_n
                                                              <= mix_hreal_0.25_mse_ratio$PQL_n * 10),]


mix_hreal_0.25_mse_ratio_p$PQL_n <- formatC(mix_hreal_0.25_mse_ratio_p$PQL_n,format = "e", digits = 0)
mix_hreal_0.25_mse_ratio_p$PQL_n <- paste0("n = ",mix_hreal_0.25_mse_ratio_p$PQL_n)
mix_hreal_0.25_mse_ratio_p$PQL_n <- ordered(mix_hreal_0.25_mse_ratio_p$PQL_n, 
                                             levels = c("n = 5e+02","n = 1e+03","n = 2e+03","n = 5e+03"))

mix_hreal_0.25_mse_ratio_p$p <- formatC(mix_hreal_0.25_mse_ratio_p$p,format = "e", digits = 0)
mix_hreal_0.25_mse_ratio_p$p <- paste0("p = ",mix_hreal_0.25_mse_ratio_p$p)
mix_hreal_0.25_mse_ratio_p$p <- ordered(mix_hreal_0.25_mse_ratio_p$p, levels = c("p = 5e+03"))

gg_mix_0.25_mse_n = ggplot(data = mix_hreal_0.25_mse_ratio_p, aes(x = MoL_n,y = mse_ratio)) + 
  geom_line(alpha = 0.7) + facet_wrap(p~PQL_n, scales = "free") + 
  ylab("MoL mse / PQL mse") + xlab("MoL sample size (n)") + labs(color = "PQL sample size (n)") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        plot.title = element_text(size = 16), 
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        strip.text = element_text(size=14)) + 
  geom_hline(aes(yintercept = 1), colour = "#DF536B") + labs(title = "Heritability = 0.25") +
  scale_x_continuous(breaks = unique(mix_hreal_0.25_mse_ratio$MoL_n), labels = scientific)

png(paste0(figureDir, '/gg_mix_0.25_mse_n.png'), width = 1000, height = 400)
gg_mix_0.25_mse_n
dev.off()

gg_mix_0.25_mse_time = ggplot(data = mix_hreal_0.25_mse_ratio_p, aes(x = MoL_n,y = time_ratio)) + 
  geom_line(alpha = 0.7) + facet_wrap(p~PQL_n, scales = "free") + 
  ylab("MoL mse / PQL mse") + xlab("MoL sample size (n)") + labs(color = "PQL sample size (n)") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
        plot.title = element_text(size = 16), 
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        strip.text = element_text(size=14)) + 
  geom_hline(aes(yintercept = 1), colour = "#DF536B") + labs(title = "Heritability = 0.25") +
  scale_x_continuous(breaks = unique(mix_hreal_0.25_mse_ratio$MoL_n), labels = scientific)

png(paste0(figureDir, '/gg_mix_0.25_mse_time.png'), width = 1000, height = 400)
gg_mix_0.25_mse_time
dev.off()

mix_hreal_0.25_mse_ratio = mix_hreal_0.25_mse_ratio[,c(3,2,7,1,8,4,5,6)]
mix_hreal_0.25_mse_ratio = mix_hreal_0.25_mse_ratio[order(mix_hreal_0.25_mse_ratio$PQL_n),]
mix_hreal_0.25_mse_ratio$p <- formatC(mix_hreal_0.25_mse_ratio$p,format = "e", digits = 0)
mix_hreal_0.25_mse_ratio$PQL_n <- formatC(mix_hreal_0.25_mse_ratio$PQL_n,format = "e", digits = 0)
mix_hreal_0.25_mse_ratio$MoL_n <- formatC(mix_hreal_0.25_mse_ratio$MoL_n,format = "e", digits = 0)

write.csv(mix_hreal_0.25_mse_ratio, paste0(tableDir,"mix_hreal_0.25_mse_ratio.csv"), quote=F, row.names = F)

## ase ratio
mix_hreal_0.25_ase_pql = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "PQLseq"),]
mix_hreal_0.25_ase_mol = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "MoL" & 
                                              mix_hreal_0.25_mse$n %in% mix_hreal_0.25_ase_pql$n),]
mix_hreal_0.25_ase_ratio = data.frame(p = mix_hreal_0.25_ase_mol$p, n = mix_hreal_0.25_ase_mol$n)
mix_hreal_0.25_ase_ratio$ase_ratio = (mix_hreal_0.25_ase_pql$mse * mix_hreal_0.25_ase_pql$time) /
  (mix_hreal_0.25_ase_mol$mse * mix_hreal_0.25_ase_mol$time)
mix_hreal_0.25_ase_ratio$n = formatC(mix_hreal_0.25_ase_ratio$n,format = "e", digits = 0)
mix_hreal_0.25_ase_ratio$p = formatC(mix_hreal_0.25_ase_ratio$p,format = "e", digits = 0)

write.csv(mix_hreal_0.25_ase_ratio, paste0(tableDir,"mix_hreal_0.25_ase_ratio.csv"), quote=F, row.names = F)
