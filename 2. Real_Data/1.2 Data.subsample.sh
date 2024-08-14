## A small data set sampling for PQL calculation:
setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data")
library(data.table)

N = c(5000,10000,20000)
size = c("s","l","ll")
for (trait in c("week_champagne","week_red_wine")){
for(i in 1:length(N)){
for (j in 1:100){
keep_trait = fread(paste0("keep_",trait,".txt"))
keep_id_sample = sample(1:nrow(keep_trait), N[i])
keep_trait_sample = keep_trait[keep_id_sample,]
fwrite(keep_trait_sample,paste0("keep_",trait,"_",size[i],"_",j,".txt"), sep = "\t")
}
}
}
