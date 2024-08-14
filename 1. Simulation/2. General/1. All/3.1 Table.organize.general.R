# Find out first 100 not NA lines
options(stringsAsFactors = F)

dlist = c(rep("norm",100),rep("mix",100))
nlist = rep(c(rep(200,20), rep(500,20), rep(1000,20), rep(2000,20),rep(5000,20)),2)
plist = rep(c(rep(500, 4), rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

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
    curr_filename = paste0(path, "PQL_result/",  d, "_n_", n, "_p_", p, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    pql_df = rbind(pql_df, df)
  }
  simulation_index[[i]] = which(pql_df$converged == TRUE)
  simulation_index[[i]] = simulation_index[[i]][1:100]
}


# MoL and PQL results summary
## original results
options(stringsAsFactors = F)

dlist = c(rep("norm",100),rep("mix",100))
nlist = rep(c(rep(200,20), rep(500,20), rep(1000,20), rep(2000,20),rep(5000,20)), 2)
plist = rep(c(rep(500, 4), rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"

## generate more for MoL
options(stringsAsFactors = F)
dlist = c(rep("norm",100), rep("mix",100))
nlist = rep(c(rep(10000,20), rep(15000,20), rep(20000,20), rep(25000,20),rep(30000,20)), 2)
plist = rep(c(rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4), rep(20000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
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
  for(j in simulation_index[[i]]){
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
  for(j in simulation_index[[i]]){
    curr_filename = paste0(path, "MoL_result/general/",  d, "_n_", n, "_p_", p, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "/")
    df = read.table(paste0(curr_filename, j, ".txt"), head=T)
    mol_df = rbind(mol_df, df)
  }
  
  write.table(mol_df, paste0(path, "MoL_result/general/original/", d, "_n_", n, "_p_", p, 
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
  write.table(pql_df, paste0(path, "PQL_result/original/", d, "_n_", n, "_p_", p, 
                           "_w_", w, "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt"),
              quote=F, row.names = F, col.names=T)
}

# MoL and PQL Heritability and computational time
library(data.table)

options(stringsAsFactors = F)
d_l = c("norm", "mix")

w_l = c(0.15, 0.5, 0.0625, 0.375)
h_l = c(0.4, 0.4, 0.8, 0.8)
hreal_l = c(1/11, 0.25, 0.2, 0.6)
sigma02_l = c(0.6,0.6,0.2,0.2)
sigmaa2_l = c(0.06,0.2,0.05,0.3)

n_l = c(rep(200,5), rep(500,5), rep(1000,5), rep(2000,5),rep(5000,5))
p_l = rep(c(500,1000,2000,5000,10000),5)

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
    MoL_file_list = paste0(path, "MoL_result/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    MoL_summary_list = lapply(MoL_file_list, fread)
    MoL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,2500), method = rep("MoL",2500),
                                           n = rep(n_l,each = 100), p = rep(p_l,each = 100),
                                           heritability = rep(NA,2500), comp_time = rep(NA,2500))
    for (k in 1:length(MoL_summary_list)){
      MoL_current_file = MoL_summary_list[[k]]
      MoL_result_list[[i]][[j]]$heritability[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$heritability
      MoL_result_list[[i]][[j]]$comp_time[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$time
    }
    
    PQL_file_list = paste0(path, "PQL_result/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    PQL_summary_list = lapply(PQL_file_list, fread)
    PQL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,2500), method = rep("PQLseq",2500),
                                           n = rep(n_l,each = 100), p = rep(p_l,each = 100),
                                           heritability = rep(NA,2500), comp_time = rep(NA,2500))
    for (k in 1:length(PQL_summary_list)){
      PQL_current_file = PQL_summary_list[[k]]
      PQL_result_list[[i]][[j]]$heritability[(100 * (k-1) + 1): (100 * k)] = PQL_current_file$h2
      PQL_result_list[[i]][[j]]$comp_time[(100 * (k-1) + 1): (100 * k)] = PQL_current_file$time
    }
    
    overall_result = rbind(MoL_result_list[[i]][[j]], PQL_result_list[[i]][[j]])
    
    write.table(overall_result, paste0(path, "overall_result/", d, "_", hreal, "_heritability.txt"),
                quote=F, row.names = F, col.names=T)
  }
}

# MoL sigma12, sigmaa2 and w
## original results
library(data.table)

options(stringsAsFactors = F)
d_l = c("norm","mix")

w_l = c(0.15, 0.5, 0.0625, 0.375)
h_l = c(0.4, 0.4, 0.8, 0.8)
sigma02_l = c(0.6,0.6,0.2,0.2)
sigmaa2_l = c(0.06,0.2,0.05,0.3)

n_l = c(rep(200,5), rep(500,5), rep(1000,5), rep(2000,5),rep(5000,5))
p_l = c(500,1000,2000,5000,10000)

N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"
result_list = list()

## generate more for MoL
library(data.table)

options(stringsAsFactors = F)
d_l = c("norm","mix")

w_l = c(0.15, 0.5, 0.0625, 0.375)
h_l = c(0.4, 0.4, 0.8, 0.8)
sigma02_l = c(0.6,0.6,0.2,0.2)
sigmaa2_l = c(0.06,0.2,0.05,0.3)

n_l = c(rep(10000,5), rep(15000,5), rep(20000,5), rep(25000,5),rep(30000,5))
p_l = c(1000,2000,5000,10000,20000)

N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"
result_list = list()

## code
for (i in 1:length(d_l)){
  d = d_l[i]
  result_list[[i]] = list()
  for (j in 1:length(w_l)){
    w = w_l[j]
    h = h_l[j]
    file_list = paste0(path, "MoL_result/general/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    summary_list = lapply(file_list, fread)
    result_list[[i]][[j]] = data.table(sigma02 = rep(NA,25), sigma02_sd = rep(NA,25), sigma02_coverage = rep(NA,25),
                                       sigmaa2 = rep(NA,25), sigmaa2_sd = rep(NA,25), sigmaa2_coverage = rep(NA,25),
                                       omega = rep(NA,25), omega_sd = rep(NA,25), 
                                       omega_trunc = rep(NA,25), omega_trunc_sd = rep(NA,25))
    for (k in 1:length(summary_list)){
      current_file = summary_list[[k]]
      result_list[[i]][[j]]$sigma02[k] = mean(current_file$sigma02)
      result_list[[i]][[j]]$sigma02_sd[k] = sqrt(var(current_file$sigma02))
      result_list[[i]][[j]]$sigma02_coverage[k] = mean(current_file$coverage02,na.rm=TRUE)
      result_list[[i]][[j]]$sigmaa2[k] = mean(current_file$sigmaa2)
      result_list[[i]][[j]]$sigmaa2_sd[k] = sqrt(var(current_file$sigmaa2))
      result_list[[i]][[j]]$sigmaa2_coverage[k] = mean(current_file$coveragea2,na.rm=TRUE)
      result_list[[i]][[j]]$omega[k] = mean(current_file$omega)
      result_list[[i]][[j]]$omega_sd[k] = sqrt(var(current_file$omega))
      for(o in 1:length(current_file$omega)){
        if(current_file$omega[o] > 1){
          current_file$omega[o] = 1
        } else if(current_file$omega[o] < 0){
          current_file$omega[o] = 0
        }
      }
      result_list[[i]][[j]]$omega_trunc[k] = mean(current_file$omega)
      result_list[[i]][[j]]$omega_trunc_sd[k] = sqrt(var(current_file$omega))
    }
    write.csv(result_list[[i]][[j]], paste0(path, "MoL_result/general/more/more_", d, "_w_", w, "_variance.csv"),
                quote=F, row.names = F, col.names=T)
    # write.csv(result_list[[i]][[j]], paste0(path, "MoL_result/general/orginal/", d, "_w_", w, "_variance.csv"),
    # quote=F, row.names = F, col.names=T)
  }
}

# MoL sigma and its variance summary
library(data.table)

options(stringsAsFactors = F)
d_l = c("norm", "mix")

w_l = c(0.15, 0.5, 0.0625, 0.375)
h_l = c(0.4, 0.4, 0.8, 0.8)
hreal_l = c(1/11, 0.25, 0.2, 0.6)
sigma02_l = c(0.6,0.6,0.2,0.2)
sigmaa2_l = c(0.06,0.2,0.05,0.3)

n_l = c(rep(200,5), rep(500,5), rep(1000,5), rep(2000,5),rep(5000,5))
p_l = rep(c(500,1000,2000,5000,10000),5)

N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"
MoL_result_list = list()

## generate more for MoL
library(data.table)

options(stringsAsFactors = F)
d_l = c("norm", "mix")

w_l = c(0.15, 0.5, 0.0625, 0.375)
h_l = c(0.4, 0.4, 0.8, 0.8)
hreal_l = c(1/11, 0.25, 0.2, 0.6)
sigma02_l = c(0.6,0.6,0.2,0.2)
sigmaa2_l = c(0.06,0.2,0.05,0.3)

n_l = c(rep(10000,5), rep(15000,5), rep(20000,5), rep(25000,5),rep(30000,5))
p_l = rep(c(1000,2000,5000,10000,20000),5)

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
    MoL_file_list = paste0(path, "MoL_result/general/more/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    #MoL_file_list = paste0(path, "MoL_result/general/original/",  d, "_n_", n_l, "_p_", p_l, "_w_", w,
    #"_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    MoL_summary_list = lapply(MoL_file_list, fread)
    MoL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,2500),
                                           n = rep(n_l,each = 100), p = rep(p_l,each = 100),
                                           sigma02 = rep(NA,2500), v02 = rep(NA,2500),
                                           sigmaa2 = rep(NA,2500), v12 = rep(NA,2500))
    for (k in 1:length(MoL_summary_list)){
      MoL_current_file = MoL_summary_list[[k]]
      MoL_result_list[[i]][[j]]$sigma02[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$sigma02
      MoL_result_list[[i]][[j]]$sigmaa2[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$sigmaa2
    }
    b1 = exp(mu + 0.5 * sigma2) * N
    b2 = sigmaa2 * exp(2 * mu + sigma2) * N^2
    b3 = exp(2 * (mu + sigma2)) * (N^2 + N)
    real_n = MoL_result_list[[i]][[j]]$n
    real_p = MoL_result_list[[i]][[j]]$p
    MoL_result_list[[i]][[j]]$v02 = rep(4 * (((sigmaa2 + 1)^2 + sigmaa2)* exp(sigma2) - 1) * (1 + 1/N) - 
      4 * ((2 * sigmaa2 + 1) * exp(2 * sigma2) - 1) * ((N^2 + 3*N + 1)/(N^2 + N)) + 
      (exp(4 * sigma2) - 1) * ((N^4 + 6*N^3 + 7*N^2 + N)/(N^2 + N)^2) + 
      (4/b1) * (exp(sigma2) * (N^2 + 3*N + 1) / (N+1)^2 - (sigmaa2 + 1)) + 2/b3,2500)
    MoL_result_list[[i]][[j]]$v12 = sigmaa2^2 * (real_n/real_p * (3/w - 1) + 4 * (exp(sigma2) * (1 + 1/N) + (b1 + b3)/b2))
 
    write.table(MoL_result_list[[i]][[j]], paste0(path, "MoL_result/general/more/var_dist_", d, "_", hreal, "_h.txt"),
                quote=F, row.names = F, col.names=T)
  }
}