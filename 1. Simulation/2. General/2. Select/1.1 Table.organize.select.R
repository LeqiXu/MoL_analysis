## d = mix, n = c(2000, 5000, 10000), p = 1000
# Find out first 100 not NA lines
options(stringsAsFactors = F)

dlist = rep("mix",4)
nlist = c(2000, 5000, 10000, 20000)
plist = rep(1000,4)
wlist = rep(0.15,4)
hlist = rep(0.4,4)
N = 10
mu = 0.2

miss = c(3,5,17,90,109,136,174,184,200)
all = c(1:200)
rest = setdiff(all,miss)

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"

simulation_index = list()

for(i in 1:length(dlist)){
  d = dlist[i]
  n = nlist[i]
  p = plist[i]
  w = wlist[i]
  h = hlist[i]
  pql_df = data.frame()
  for(j in 1:200){
    if(!(j %in% miss)){
      curr_filename = paste0(path, "PQL_result/",  d, "_n_", n, "_p_", p, "_w_", w,
                             "_h_", h, "_N_", N, "_mu_", mu, "/")
      df = read.table(paste0(curr_filename, j, ".txt"), head=T)
      pql_df = rbind(pql_df, df) 
    } else{
      pql_df = rbind(pql_df,NA)
    }
      
  }
  simulation_index[[i]] = which(pql_df$converged == TRUE)
  simulation_index[[i]] = simulation_index[[i]][1:50]
}


# MoL and PQL results summary
# p = 500, n = 50000

options(stringsAsFactors = F)

dlist = rep("mix",4)
nlist = c(50000,100000,200000,500000)
plist = rep(5000,4)
wlist = rep(0.5,4)
hlist = rep(0.4,4)
N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"

## check
for(i in 1:length(dlist)){
  d = dlist[i]
  n = nlist[i]
  p = plist[i]
  w = wlist[i]
  h = hlist[i]
  print(i)
  mol_df = data.frame()
  #for(j in simulation_index[[i]]){
  for(j in 1:100){
    curr_filename = paste0(path, "MoL_result/general/",  d, "_n_", n, "_p_", p, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    mol_df = rbind(mol_df, df)
  }
}

## code
for(i in 1:length(dlist)){
  d = dlist[i]
  n = nlist[i]
  p = plist[i]
  w = wlist[i]
  h = hlist[i]
  mol_df = data.frame()
  #for(j in simulation_index[[i]]){
  for(j in 1:100){
    curr_filename = paste0(path, "MoL_result/general/",  d, "_n_", n, "_p_", p, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    mol_df = rbind(mol_df, df)
  }
  
  write.table(mol_df, paste0(path, "MoL_result/general/", d, "_n_", n, "_p_", p, 
                             "_w_", w, "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt"),
              quote=F, row.names = F, col.names=T)
}


for(i in 1:length(dlist)){
  d = dlist[i]
  n = nlist[i]
  p = plist[i]
  w = wlist[i]
  h = hlist[i]
  pql_df = data.frame()
  for(j in simulation_index[[i]]){
    curr_filename = paste0(path, "PQL_result/",  d, "_n_", n, "_p_", p, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    pql_df = rbind(pql_df, df)
  }
  write.table(pql_df, paste0(path, "PQL_result/select/", d, "_n_", n, "_p_", p, 
                             "_w_", w, "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt"),
              quote=F, row.names = F, col.names=T)
}

# MoL and PQL Heritability and computational time
library(data.table)

options(stringsAsFactors = F)
d_l = c("mix")

w_l = c(0.15)
h_l = c(0.4)
hreal_l = c(1/11)
sigma02_l = c(0.6)
sigmaa2_l = c(0.06)

n_l = c(2000, 5000, 10000,20000)
p_l = rep(1000,4)

N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"
MoL_result_list = list()
PQL_result_list = list()

for (i in 1:length(d_l)){
  d = d_l[i]
  MoL_result_list[[i]] = list()
  PQL_result_list[[i]] = list()
  for (j in 1:length(hreal_l)){
    w = w_l[j]
    h = h_l[j]
    hreal = hreal_l[j]
    MoL_file_list = paste0(path, "MoL_result/general/select/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    MoL_summary_list = lapply(MoL_file_list, fread)
    MoL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,200), method = rep("MoL",200),
                                           n = rep(n_l,each = 50), p = rep(p_l,each = 50),
                                           heritability = rep(NA,200), comp_time = rep(NA,200))
    for (k in 1:length(MoL_summary_list)){
      MoL_current_file = MoL_summary_list[[k]]
      MoL_result_list[[i]][[j]]$heritability[(50 * (k-1) + 1): (50 * k)] = MoL_current_file$heritability
      MoL_result_list[[i]][[j]]$comp_time[(50 * (k-1) + 1): (50 * k)] = MoL_current_file$time
    }
    
    PQL_file_list = paste0(path, "PQL_result/select/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    PQL_summary_list = lapply(PQL_file_list, fread)
    PQL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,200), method = rep("PQLseq",200),
                                           n = rep(n_l,each = 50), p = rep(p_l,each = 50),
                                           heritability = rep(NA,200), comp_time = rep(NA,200))
    for (k in 1:length(PQL_summary_list)){
      PQL_current_file = PQL_summary_list[[k]]
      PQL_result_list[[i]][[j]]$heritability[(50 * (k-1) + 1): (50 * k)] = PQL_current_file$h2
      PQL_result_list[[i]][[j]]$comp_time[(50 * (k-1) + 1): (50 * k)] = PQL_current_file$time
    }
    
    overall_result = rbind(MoL_result_list[[i]][[j]], PQL_result_list[[i]][[j]])
    
    write.table(overall_result, paste0(path, "overall_result/select/", d, "_", hreal, "_heritability.txt"),
                quote=F, row.names = F, col.names=T)
  }
}

# MoL sigma and its variance summary
library(data.table)

options(stringsAsFactors = F)
d_l = c("mix")

w_l = c(0.15)
h_l = c(0.4)
hreal_l = c(1/11)
sigma02_l = c(0.6)
sigmaa2_l = c(0.06)
n_l = c(2000, 5000, 10000,20000)
p_l = rep(1000,4)

N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"
MoL_result_list = list()

## code
for (i in 1:length(d_l)){
  d = d_l[i]
  MoL_result_list[[i]] = list()
  
  for (j in 1:length(hreal_l)){
    w = w_l[j]
    h = h_l[j]
    sigma02 = sigma02_l[j]
    sigmaa2 = sigmaa2_l[j]
    sigma2 = sigma02 + sigmaa2
    hreal = hreal_l[j]
    MoL_file_list = paste0(path, "MoL_result/general/select/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    #MoL_file_list = paste0(path, "MoL_result/general/original/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
    #"_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    MoL_summary_list = lapply(MoL_file_list, fread)
    MoL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,200),
                                           n = rep(n_l,each = 50), p = rep(p_l,each = 50),
                                           sigma02 = rep(NA,200), v02 = rep(NA,200),
                                           sigmaa2 = rep(NA,200), v12 = rep(NA,200))
    for (k in 1:length(MoL_summary_list)){
      MoL_current_file = MoL_summary_list[[k]]
      MoL_result_list[[i]][[j]]$sigma02[(50 * (k-1) + 1): (50 * k)] = MoL_current_file$sigma02
      MoL_result_list[[i]][[j]]$sigmaa2[(50 * (k-1) + 1): (50 * k)] = MoL_current_file$sigmaa2
    }
    b1 = exp(mu + 0.5 * sigma2) * N
    b2 = sigmaa2 * exp(2 * mu + sigma2) * N^2
    b3 = exp(2 * (mu + sigma2)) * (N^2 + N)
    real_n = MoL_result_list[[i]][[j]]$n
    real_p = MoL_result_list[[i]][[j]]$p
    MoL_result_list[[i]][[j]]$v02 = rep(4 * (((sigmaa2 + 1)^2 + sigmaa2)* exp(sigma2) - 1) * (1 + 1/N) - 
                                          4 * ((2 * sigmaa2 + 1) * exp(2 * sigma2) - 1) * ((N^2 + 3*N + 1)/(N^2 + N)) + 
                                          (exp(4 * sigma2) - 1) * ((N^4 + 6*N^3 + 7*N^2 + N)/(N^2 + N)^2) + 
                                          (4/b1) * (exp(sigma2) * (N^2 + 3*N + 1) / (N+1)^2 - (sigmaa2 + 1)) + 2/b3,200)
    MoL_result_list[[i]][[j]]$v12 = sigmaa2^2 * (real_n/real_p * (3/w - 1) + 4 * (exp(sigma2) * (1 + 1/N) + (b1 + b3)/b2))
    
    write.table(MoL_result_list[[i]][[j]], paste0(path, "MoL_result/general/select/var_dist_", d, "_", hreal, "_h.txt"),
                quote=F, row.names = F, col.names=T)
  }
}