library(ggplot2)

options(stringsAsFactors = F)

nlist = rep(c(5000, 2000, rep(1000, 5), 500, 200), 2)
plist = rep(c(rep(2000, 3), c(10000, 5000, 1000, 500), rep(2000, 2)), 2)
hlist = c(rep(0.3, 9), rep(0.6, 9))
N = 10

path = "/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/"

for(i in 1:length(nlist)){
  n = nlist[i]
  p = plist[i]
  h = hlist[i]
  mol_df = data.frame()
  pql_df = data.frame()
  for(j in 1:100){
    curr_filename = paste0(path, "MOL/n_",  n, "_p_", p, "_h_", h, "_N_", N, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    mol_df = rbind(mol_df, df)
    curr_filename = paste0(path, "PQL/n_",  n, "_p_", p, "_h_", h, "_N_", N, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    pql_df = rbind(pql_df, df)
  }
  
  write.table(mol_df, paste0(path, "MOL/n_",  n, "_p_", p, "_h_", h, "_N_", N, 
                             "/n_", n, "_p_", p, "_h_", h, "_N_", N, "_summary.txt"),
              quote=F, row.names = F, col.names=T)
  
  write.table(pql_df, paste0(path, "PQL/n_",  n, "_p_", p, "_h_", h, "_N_", N, 
                             "/n_", n, "_p_", p, "_h_", h, "_N_", N, "_summary.txt"),
              quote=F, row.names = F, col.names=T)
}

##########################################################################

job = c()
for(i in 1:length(nlist)){
  n = nlist[i]
  p = plist[i]
  h = hlist[i]
  job = c(job, paste0("scp yz738@farnam.hpc.yale.edu:",
                      path, "MOL/n_",  n, "_p_", p, "_h_", h, "_N_", N, "/",
                      "n_", n, "_p_", p, "_h_", h, "_N_", N, "_summary.txt",
                      " ./mol_",
                      "n_", n, "_p_", p, "_h_", h, "_N_", N, "_summary.txt"))
  job = c(job, paste0("scp yz738@farnam.hpc.yale.edu:",
                      path, "PQL/n_",  n, "_p_", p, "_h_", h, "_N_", N, "/",
                      "n_", n, "_p_", p, "_h_", h, "_N_", N, "_summary.txt",
                      " ./pql_",
                      "n_", n, "_p_", p, "_h_", h, "_N_", N, "_summary.txt"))
}

write.table(job[20:36], "/Users/yilinagzhang/Documents/MolGWAS/Results/simpleSimulation/copy.txt",
            quote=F, col.names=F, row.names = F)

##########################################################################

library(ggplot2)
library(ggthemes)
library(scales)
options(stringsAsFactors=F)
options(scipen = 200)

mainDir= '/Users/yilinagzhang/Documents/MolGWAS/Figure/Simple'
resDir = '/Users/yilinagzhang/Documents/MolGWAS/Results/simpleSimulation'

# boxplot h 0.3 p 2000

mol_df_n_200_p_2000_h_0.3 = read.table(paste0(resDir, "/mol_n_200_p_2000_h_0.3_N_10_summary.txt"),
                                       head =T)
pql_df_n_200_p_2000_h_0.3 = read.table(paste0(resDir, "/pql_n_200_p_2000_h_0.3_N_10_summary.txt"),
                                       head =T)
mol_df_n_500_p_2000_h_0.3 = read.table(paste0(resDir, "/mol_n_500_p_2000_h_0.3_N_10_summary.txt"),
                                       head =T)
pql_df_n_500_p_2000_h_0.3 = read.table(paste0(resDir, "/pql_n_500_p_2000_h_0.3_N_10_summary.txt"),
                                       head =T)
mol_df_n_1000_p_2000_h_0.3 = read.table(paste0(resDir, "/mol_n_1000_p_2000_h_0.3_N_10_summary.txt"),
                                       head =T)
pql_df_n_1000_p_2000_h_0.3 = read.table(paste0(resDir, "/pql_n_1000_p_2000_h_0.3_N_10_summary.txt"),
                                       head =T)
mol_df_n_2000_p_2000_h_0.3 = read.table(paste0(resDir, "/mol_n_2000_p_2000_h_0.3_N_10_summary.txt"),
                                        head =T)
pql_df_n_2000_p_2000_h_0.3 = read.table(paste0(resDir, "/pql_n_2000_p_2000_h_0.3_N_10_summary.txt"),
                                        head =T)
mol_df_n_5000_p_2000_h_0.3 = read.table(paste0(resDir, "/mol_n_5000_p_2000_h_0.3_N_10_summary.txt"),
                                        head =T)
pql_df_n_5000_p_2000_h_0.3 = read.table(paste0(resDir, "/pql_n_5000_p_2000_h_0.3_N_10_summary.txt"),
                                        head =T)

data_to_plot = data.frame(n=c(rep("200", 100), rep("500", 100), rep("1000", 100),
                              rep("2000", 100), rep("5000", 100),
                              rep("200", 100), rep("500", 100), rep("1000", 100),
                              rep("2000", 100), rep("5000", 100)),
                          method=c(rep("MOL", 500), rep("PQLseq", 500)),
                          estimate=c(mol_df_n_200_p_2000_h_0.3$heritability,
                                     mol_df_n_500_p_2000_h_0.3$heritability,
                                     mol_df_n_1000_p_2000_h_0.3$heritability,
                                     mol_df_n_2000_p_2000_h_0.3$heritability,
                                     mol_df_n_5000_p_2000_h_0.3$heritability,
                                     pql_df_n_200_p_2000_h_0.3$h2,
                                     pql_df_n_500_p_2000_h_0.3$h2,
                                     pql_df_n_1000_p_2000_h_0.3$h2,
                                     pql_df_n_2000_p_2000_h_0.3$h2,
                                     pql_df_n_5000_p_2000_h_0.3$h2))

data_to_plot$method <- ordered(data_to_plot$method, levels=c("MOL", "PQLseq"))
data_to_plot$n <- ordered(data_to_plot$n, levels=c("200", "500", "1000", "2000",
                                                   "5000"))
pdf(paste0(mainDir, '/Boxplot_h_0.3_N_10_p_2000.pdf'), width = 15, height = 18)
gg = ggplot(data_to_plot, aes(n, estimate))
gg = gg + geom_boxplot(aes(fill=method), outlier.shape = NA, alpha=0.7)
gg = gg + theme_minimal(base_family="Helvetica")
gg = gg + scale_fill_manual(labels=c("MOL", "PQLseq"),values=c("#2e6cd1", "#d1972e"))
gg = gg + ylab("estimated heritability") + xlab("sample size")
gg = gg + theme(axis.text=element_text(size=42), axis.title=element_text(size=48))
gg = gg + geom_hline(yintercept = 0.3, linetype="dashed", color="red")
# gg = gg + theme(legend.title=element_text(size=42), legend.text=element_text(size=42), legend.position=c(1,0), legend.justification=c(1,0), legend.text.align = 0)
# gg = gg + scale_y_continuous(breaks=seq(-1, 1, 0.002))
# gg = gg + coord_cartesian(ylim = c(-0.008,0.014))
gg
dev.off()

##########################################################################

library(ggplot2)
library(ggthemes)
library(scales)
options(stringsAsFactors=F)
options(scipen = 200)

mainDir= '/Users/yilinagzhang/Documents/MolGWAS/Figure/Simple'
resDir = '/Users/yilinagzhang/Documents/MolGWAS/Results/simpleSimulation'

# boxplot h 0.3 n 1000

mol_df_n_1000_p_500_h_0.3 = read.table(paste0(resDir, "/mol_n_1000_p_500_h_0.3_N_10_summary.txt"),
                                       head =T)
pql_df_n_1000_p_500_h_0.3 = read.table(paste0(resDir, "/pql_n_1000_p_500_h_0.3_N_10_summary.txt"),
                                       head =T)
mol_df_n_1000_p_1000_h_0.3 = read.table(paste0(resDir, "/mol_n_1000_p_1000_h_0.3_N_10_summary.txt"),
                                       head =T)
pql_df_n_1000_p_1000_h_0.3 = read.table(paste0(resDir, "/pql_n_1000_p_1000_h_0.3_N_10_summary.txt"),
                                       head =T)
mol_df_n_1000_p_2000_h_0.3 = read.table(paste0(resDir, "/mol_n_1000_p_2000_h_0.3_N_10_summary.txt"),
                                        head =T)
pql_df_n_1000_p_2000_h_0.3 = read.table(paste0(resDir, "/pql_n_1000_p_2000_h_0.3_N_10_summary.txt"),
                                        head =T)
mol_df_n_1000_p_5000_h_0.3 = read.table(paste0(resDir, "/mol_n_1000_p_5000_h_0.3_N_10_summary.txt"),
                                        head =T)
pql_df_n_1000_p_5000_h_0.3 = read.table(paste0(resDir, "/pql_n_1000_p_5000_h_0.3_N_10_summary.txt"),
                                        head =T)
mol_df_n_1000_p_10000_h_0.3 = read.table(paste0(resDir, "/mol_n_1000_p_10000_h_0.3_N_10_summary.txt"),
                                        head =T)
pql_df_n_1000_p_10000_h_0.3 = read.table(paste0(resDir, "/pql_n_1000_p_10000_h_0.3_N_10_summary.txt"),
                                        head =T)

data_to_plot = data.frame(p=c(rep("500", 100), rep("1000", 100), rep("2000", 100),
                              rep("5000", 100), rep("10000", 100),
                              rep("500", 100), rep("1000", 100), rep("2000", 100),
                              rep("5000", 100), rep("10000", 100)),
                          method=c(rep("MOL", 500), rep("PQLseq", 500)),
                          estimate=c(mol_df_n_1000_p_500_h_0.3$heritability,
                                     mol_df_n_1000_p_1000_h_0.3$heritability,
                                     mol_df_n_1000_p_2000_h_0.3$heritability,
                                     mol_df_n_1000_p_5000_h_0.3$heritability,
                                     mol_df_n_1000_p_10000_h_0.3$heritability,
                                     pql_df_n_1000_p_500_h_0.3$h2,
                                     pql_df_n_1000_p_1000_h_0.3$h2,
                                     pql_df_n_1000_p_2000_h_0.3$h2,
                                     pql_df_n_1000_p_5000_h_0.3$h2,
                                     pql_df_n_1000_p_10000_h_0.3$h2))

data_to_plot$method <- ordered(data_to_plot$method, levels=c("MOL", "PQLseq"))
data_to_plot$p <- ordered(data_to_plot$p, levels=c("500", "1000", "2000", "5000",
                                                   "10000"))
pdf(paste0(mainDir, '/Boxplot_h_0.3_N_10_n_1000.pdf'), width = 15, height = 18)
gg = ggplot(data_to_plot, aes(p, estimate))
gg = gg + geom_boxplot(aes(fill=method), outlier.shape = NA, alpha=0.7)
gg = gg + theme_minimal(base_family="Helvetica")
gg = gg + scale_fill_manual(labels=c("MOL", "PQLseq"),values=c("#2e6cd1", "#d1972e"))
gg = gg + ylab("estimated heritability") + xlab("number of markers")
gg = gg + theme(axis.text=element_text(size=42), axis.title=element_text(size=48))
gg = gg + geom_hline(yintercept = 0.3, linetype="dashed", color="red")
# gg = gg + theme(legend.title=element_text(size=42), legend.text=element_text(size=42), legend.position=c(1,0), legend.justification=c(1,0), legend.text.align = 0)
# gg = gg + scale_y_continuous(breaks=seq(-1, 1, 0.002))
# gg = gg + coord_cartesian(ylim = c(-0.008,0.014))
gg
dev.off()
