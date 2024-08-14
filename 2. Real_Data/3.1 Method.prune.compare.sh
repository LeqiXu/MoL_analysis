setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/result/prune")
library(data.table)

## organize the table
dist=100
thres=0.001

total_df_all = c()

for(trait in c("week_champagne","week_red_wine")){
    PQLseq_prune_s = lapply(paste0("PQLseq_prune_",dist,"_",thres,"_",trait,"_s_",1:100,"_result.txt"), fread)
    PQLseq_prune_l = lapply(paste0("PQLseq_prune_",dist,"_",thres,"_",trait,"_l_",1:100,"_result.txt"), fread)

    h2_s = rep(0,100)
    h2_l = rep(0,100)

    for(i in 1:100){
    h2_s[i] = PQLseq_prune_s[[i]]$h2
    h2_l[i] = PQLseq_prune_l[[i]]$h2
    }

    total_df_s = data.table(trait=trait,distance=dist,threshold=thres,sample_size=5000,heritability=h2_s)
    total_df_l = data.table(trait=trait,distance=dist,threshold=thres,sample_size=10000,heritability=h2_l)
    total_df = rbind(total_df_s,total_df_l)
    total_df_all = rbind(total_df_all,total_df)
}

write.table(total_df_all, file="/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/result/PQLseq_prune_table.tsv", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# plot1
setwd("/Users/xuleqi/Desktop/Research/Project/3.1 Count_Data_GWAS_Method：MoL/Result/Real_data/")
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

trait = c("week_champagne","week_red_wine")
PQLseq_prune_total = fread("PQLseq_prune_table.tsv")
MoL_prune_total = fread("MoL_prune_table.tsv")
MoL_prune_total = MoL_prune_total[which(MoL_prune_total$threshold == 0.001),]

PQLseq_prune_total$type = paste0(PQLseq_prune_total$trait,"_",PQLseq_prune_total$distance,"_",PQLseq_prune_total$threshold)
MoL_prune_total$type = paste0(MoL_prune_total$trait,"_",MoL_prune_total$distance,"_",MoL_prune_total$threshold)

ggplot(PQLseq_prune_total, aes(x = sample_size, y = heritability, group = sample_size)) + theme(legend.position = "bottom") +
                    geom_boxplot() + scale_fill_brewer(palette = "Set2") + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
                    geom_hline(data = MoL_prune_total, aes(yintercept = heritability), color = "red", linetype = "dashed") + facet_wrap(~type)

write.table(MoL_prune_total, file="/Users/xuleqi/Desktop/Research/Project/3.1 Count_Data_GWAS_Method：MoL/Result/Real_data/MoL_prune_table_0.001.tsv", append=F, sep="\t", quote=F, row.names=F, col.names=T)

# plot2
setwd("/Users/xuleqi/Desktop/Research/Project/3.1 Count_Data_GWAS_Method：MoL/Result/Real_data/")

library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

trait = c("week champagne","week red wine")
PQLseq_prune_total = fread("PQLseq_prune_table.tsv")
MoL_prune_total = fread("MoL_prune_table.tsv")
MoL_prune_total = MoL_prune_total[which(MoL_prune_total$threshold == 0.001),]

PQLseq_prune_total = PQLseq_prune_total[which(PQLseq_prune_total$distance==100),]
MoL_prune_total = MoL_prune_total[which(MoL_prune_total$distance == 100),]

## MoL table
MoL_table = MoL_prune_total[,c("trait","sample_size","heritability")]
MoL_table$trait = str_replace_all(MoL_table$trait,"_"," ")
colnames(MoL_table)[3] = "h2"
MoL_table$method = "MoL"
MoL_table$std_h2 = sqrt(MoL_prune_total$v_h2)
MoL_table = MoL_table[,c("method","trait","sample_size","h2","std_h2")]

## PQL_table
PQL_table = PQLseq_prune_total[,.("h2" = mean(heritability), "std_h2" = sd(heritability)), by = c("trait","sample_size")]
PQL_table$trait = str_replace_all(PQL_table$trait,"_"," ")
PQL_table$method = "PQLseq"
PQL_table = PQL_table[,c("method","trait","sample_size","h2","std_h2")]

final_table = rbind(MoL_table,PQL_table)
final_table$sample_size = factor(final_table$sample_size,levels=c(5000,10000,111350,133609))

ggplot(final_table, aes(x=sample_size, y=h2, group=method, color=method)) + 
  geom_point()+
  geom_errorbar(aes(ymin=h2-1.96*std_h2, ymax=h2+1.96*std_h2), width=.2) + 
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~trait,scales = "free_x") + theme_bw() + theme(strip.background = element_rect(color = "white", fill = "white"))