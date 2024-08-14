library(ggpubr)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

options(stringsAsFactors=F)
options(scipen = 0)

tableDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Table/Heritability/'
figureDir= '/Users/xuleqi/Desktop/Research/Project/2.1 Count_Data_GWAS_Project：MoL/Result/Simulation/General/Select/Figure/Heritability/'

mix_hreal_0.25 = read.table(paste0(tableDir,"fix_mix_0.25_heritability.txt"), head =T)

mix_hreal_0.25 <- na.omit(mix_hreal_0.25)
mix_hreal_0.25$real_h2 = 0.25
mix_hreal_0.25_mse = mix_hreal_0.25 %>% group_by(method,n,p) %>% mutate(mse = 1 / 100 * sum((heritability - real_h2)^2)) %>% 
  mutate(bias2 = 1 / 100 *sum((mean(heritability) - real_h2)^2)) %>% 
  mutate(variance = 1 / 100 * sum((heritability - mean(heritability))^2)) %>%
  mutate(time = mean(comp_time)) %>%
  select(real_h2, method, n, p, mse, bias2, variance,time) %>% unique()

# mse ratio
mix_hreal_0.25_mse_mol = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "MoL"),]
mix_hreal_0.25_mse_pql = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "PQLseq"),]

mix_hreal_0.25_mse_ratio = mix_hreal_0.25_mse_mol[,c("n","p")]
colnames(mix_hreal_0.25_mse_ratio) = c("MoL_n","p")
mix_hreal_0.25_mse_ratio$PQL_n = rep(mix_hreal_0.25_mse_pql$n, 3)
mix_hreal_0.25_mse_ratio$var_ratio = mix_hreal_0.25_mse_mol$variance / rep(mix_hreal_0.25_mse_pql$variance, 3)
mix_hreal_0.25_mse_ratio$mse_ratio = mix_hreal_0.25_mse_mol$mse / rep(mix_hreal_0.25_mse_pql$mse, 3)
mix_hreal_0.25_mse_ratio$time_ratio = mix_hreal_0.25_mse_mol$time / rep(mix_hreal_0.25_mse_pql$time, 3)
mix_hreal_0.25_mse_ratio = merge(mix_hreal_0.25_mse_ratio, mix_hreal_0.25_mse_pql[,c("n", "variance")], 
                                  by.x = c("PQL_n"), by.y = c("n"))
mix_hreal_0.25_mse_ratio = merge(mix_hreal_0.25_mse_ratio, mix_hreal_0.25_mse_mol[,c("n","p","variance")], 
                                  by.x = c("MoL_n","p"), by.y = c("n","p"))
names(mix_hreal_0.25_mse_ratio)[7:8] = c("PQL_variance", "MoL_variance")

mix_hreal_0.25_mse_ratio = mix_hreal_0.25_mse_ratio[,c(2,3,7,1,8,4,5,6)]
mix_hreal_0.25_mse_ratio = mix_hreal_0.25_mse_ratio[order(mix_hreal_0.25_mse_ratio$MoL_n, 
                                                            -1 * (mix_hreal_0.25_mse_ratio$p)),]
mix_hreal_0.25_mse_ratio$p <- formatC(mix_hreal_0.25_mse_ratio$p,format = "e", digits = 0)
mix_hreal_0.25_mse_ratio$PQL_n <- formatC(mix_hreal_0.25_mse_ratio$PQL_n,format = "e", digits = 0)
mix_hreal_0.25_mse_ratio$MoL_n <- formatC(mix_hreal_0.25_mse_ratio$MoL_n,format = "e", digits = 0)

write.csv(mix_hreal_0.25_mse_ratio, paste0(tableDir,"fix_mix_hreal_0.25_mse_ratio.csv"), quote=F, row.names = F)

## ase ratio
mix_hreal_0.25_ase_pql = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "PQLseq"),]
mix_hreal_0.25_ase_mol = mix_hreal_0.25_mse[which(mix_hreal_0.25_mse$method == "MoL" & 
                                              paste0(mix_hreal_0.25_mse$n,mix_hreal_0.25_mse$p) %in% 
                                              paste0(mix_hreal_0.25_ase_pql$n,mix_hreal_0.25_ase_pql$p)),]
mix_hreal_0.25_ase_ratio = data.frame(p = mix_hreal_0.25_ase_mol$p, n = mix_hreal_0.25_ase_mol$n)
mix_hreal_0.25_ase_ratio$ase_ratio = (mix_hreal_0.25_ase_pql$mse * mix_hreal_0.25_ase_pql$time) /
  (mix_hreal_0.25_ase_mol$mse * mix_hreal_0.25_ase_mol$time)
mix_hreal_0.25_ase_ratio$n = formatC(mix_hreal_0.25_ase_ratio$n,format = "e", digits = 0)
mix_hreal_0.25_ase_ratio$p = formatC(mix_hreal_0.25_ase_ratio$p,format = "e", digits = 0)

write.csv(mix_hreal_0.25_ase_ratio, paste0(tableDir,"fix_mix_hreal_0.25_ase_ratio.csv"), quote=F, row.names = F)
