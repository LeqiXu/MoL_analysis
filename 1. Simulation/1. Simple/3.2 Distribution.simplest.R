library(data.table)
library(ggplot2)

file_path = "/Users/xuleqi/Desktop/Research/Project/Count_Data_GWAS_Project：MoL/Result/Simulation/Simple/Table"
setwd(file_path)
file_name = dir(pattern = "mol")
sum_stat = lapply(file_name,fread)

n = c(rep(1000,10),rep(200,2),rep(2000,2),rep(500,2),rep(5000,2))
p = c(rep(c(1000,10000,2000,500,5000), each = 2), rep(2000,8))
h = c(rep(c(0.3,0.6),9))
mu = 0
N = 10
w = 1
sigma02 = 1 - h
sigmaa2 = h
sigma2 = 1

v02 = 4 * exp(sigma2) * (sigma2^2 + 3 * sigmaa2 + 1 - (2 * sigmaa2 + 1) * exp(sigma2)) + exp(4 * sigma2) - 1 + 4 / (exp(0.5 * sigma2) * N) * (exp(sigma2) - sigmaa2 - 1) + 2 / (exp(2 * sigma2) * N^2)
va2 = 2 * sigmaa2^2 * (n/p + 2 * (exp(sigma2) + (1 + exp(1.5 * sigma2))/(sigmaa2 * exp(0.5 * sigma2) * N^3)))

total_sum = data.table(true_heritability = rep(h,each = 100),
                       n = rep(n,each = 100), p = rep(p,each = 100),
                       sigma02 = rep(NA,1800), v02 = rep(NA,1800),
                       sigmaa2 = rep(NA,1800), va2 = rep(NA,1800),
                       z = rep(qnorm(ppoints(100)),18))

for(k in 1:length(sum_stat)){
  total_sum$sigma02[(100 * (k-1) + 1): (100 * k)] = sum_stat[[k]]$sigma02[order(sum_stat[[k]]$sigma02)]
  total_sum$sigmaa2[(100 * (k-1) + 1): (100 * k)] = sum_stat[[k]]$sigma12[order(sum_stat[[k]]$sigma12)]
  total_sum$v02[(100 * (k-1) + 1): (100 * k)] = v02[k]
  total_sum$va2[(100 * (k-1) + 1): (100 * k)] = va2[k]
}

total_sum$sigma02_z = sqrt(total_sum$n) * (total_sum$sigma02 - (1 - total_sum$true_heritability)) / sqrt(total_sum$v02)
total_sum$sigmaa2_z = sqrt(total_sum$n) * (total_sum$sigmaa2 - total_sum$true_heritability) / sqrt(total_sum$va2)

total_sum$n = as.factor(total_sum$n)
total_sum$p = as.factor(total_sum$p)

p1 <- ggplot(total_sum, aes(x=z, y=sigma02_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(true_heritability~p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigma02 distribution under normal assumption")
print(p1)

p2 <- ggplot(total_sum, aes(x=z, y=sigmaa2_z)) + 
  geom_abline(slope = 1, intercept = c(0,0), size=0.8) +
  geom_point(size=3, alpha=0.7, col="#ed7864") + 
  theme(plot.title = element_text(size=18, hjust=0.5, vjust=0.5)) + 
  xlab("expected") + ylab("observed") + 
  facet_grid(true_heritability~p~n,scales = "free") +
  theme(axis.title=element_text(size=16), axis.text = element_text(size=12),
        strip.text = element_text(size=14)) +
  ggtitle("Sigmaa2 distribution under normal assumption")
print(p2)