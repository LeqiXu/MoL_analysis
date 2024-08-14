options(stringsAsFactors = F)
nlist = c(500)
dlist = rep("mix",length(nlist))
plist = rep(5000,length(nlist))
wlist = rep(0.5,length(nlist))
hlist = rep(0.4,length(nlist))
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

## code
for(i in 1:length(dlist)){
  d = dlist[i]
  n = nlist[i]
  p = plist[i]
  w = wlist[i]
  h = hlist[i]
  mol_df = data.frame()
  #for(j in 1:100){
  for(j in simulation_index[[i]]){
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
  write.table(pql_df, paste0(path, "PQL_result/", d, "_n_", n, "_p_", p, 
                             "_w_", w, "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt"),
              quote=F, row.names = F, col.names=T)
}

# MoL sigma table
# MoL sigma and its variance summary
library(data.table)

options(stringsAsFactors = F)
d_l = c("mix")

w_l = c(0.5)
h_l = c(0.4)
hreal_l = c(0.25)
sigma02_l = c(0.6)
sigmaa2_l = c(0.2)

#n_l_mol = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)
#p_l_mol = rep(5000, length(n_l_mol))

n_l_mol = c(500, 1000, 2000, 5000, 5000, 10000, 20000, 50000, 50000, 100000, 200000, 500000)
p_l_mol = rep(c(500,1000,2000,5000),3)

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
    MoL_file_list = paste0(path, "MoL_result/general/",  d, "_n_", n_l_mol, "_p_", p_l_mol, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    MoL_summary_list = lapply(MoL_file_list, fread)
    MoL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,100 * length(n_l_mol)),
                                           n = rep(n_l_mol,each = 100), p = rep(p_l_mol,each = 100),
                                           sigma02 = rep(NA,100 * length(n_l_mol)), v02 = rep(NA,100 * length(n_l_mol)),
                                           sigmaa2 = rep(NA,100 * length(n_l_mol)), v12 = rep(NA,100 * length(n_l_mol)))
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
                                          (4/b1) * (exp(sigma2) * (N^2 + 3*N + 1) / (N+1)^2 - (sigmaa2 + 1)) + 2/b3, 100 * length(n_l_mol))
    MoL_result_list[[i]][[j]]$v02 = 1 / real_n * MoL_result_list[[i]][[j]]$v02
    MoL_result_list[[i]][[j]]$v12 = sigmaa2^2 * (real_n/real_p * (3/w - 1) + 4 * (exp(sigma2) * (1 + 1/N) + (b1 + b3)/b2))
    MoL_result_list[[i]][[j]]$v12 = 1 / real_n * MoL_result_list[[i]][[j]]$v12
    #write.table(MoL_result_list[[i]][[j]], paste0(path, "MoL_result/general/select/var_dist_", d, "_", hreal, "_h.txt"),
    #            quote=F, row.names = F, col.names=T)
    write.table(MoL_result_list[[i]][[j]], paste0(path, "MoL_result/general/select/fix_var_dist_", d, "_", hreal, "_h.txt"),
                quote=F, row.names = F, col.names=T)
  }
}

## generate summary table for MoL
library(data.table)

options(stringsAsFactors = F)
d_l = c("mix")

w_l = c(0.5)
h_l = c(0.4)
hreal_l = c(0.25)
sigma02_l = c(0.6)
sigmaa2_l = c(0.2)

#n_l_mol = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)
#p_l_mol = rep(5000, length(n_l_mol))

n_l_mol = c(500, 1000, 2000, 5000, 5000, 10000, 20000, 50000, 50000, 100000, 200000, 500000)
p_l_mol = rep(c(500,1000,2000,5000),3)

N = 10
mu = 0.2

path = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/"
result_list = list()

## code
for (i in 1:length(d_l)){
  d = d_l[i]
  result_list[[i]] = list()
  
  #var_file = fread(paste0(path, "MoL_result/general/select/var_dist_", d, "_", hreal_l[i], "_h.txt"))
  var_file = fread(paste0(path, "MoL_result/general/select/fix_var_dist_", d, "_", hreal_l[i], "_h.txt"))
  var_file = unique(var_file[,c("n","p","v02","v12")])
  for (j in 1:length(w_l)){
    w = w_l[j]
    h = h_l[j]
    file_list = paste0(path, "MoL_result/general/",  d, "_n_", n_l_mol, "_p_", p_l_mol, "_w_", w,
                       "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    summary_list = lapply(file_list, fread)
    result_list[[i]][[j]] = data.table(sigma02 = rep(NA,length(n_l_mol)), sigma02_est.sd = rep(NA,length(n_l_mol)),
                                       sigma02_emp.sd = rep(NA,length(n_l_mol)), sigma02_coverage = rep(NA,length(n_l_mol)), 
                                       sigmaa2 = rep(NA,length(n_l_mol)), sigmaa2_est.sd = rep(NA,length(n_l_mol)),
                                       sigmaa2_emp.sd = rep(NA,length(n_l_mol)), sigmaa2_coverage = rep(NA,length(n_l_mol)),
                                       omega = rep(NA,length(n_l_mol)), omega_emp.sd = rep(NA,length(n_l_mol)), 
                                       omega_trunc = rep(NA,length(n_l_mol)), omega_trunc_emp.sd = rep(NA,length(n_l_mol)))
    for (k in 1:length(summary_list)){
      current_file = summary_list[[k]]
      current_file$coverage02[which(is.na(current_file$coverage02))] = 0
      current_file$coveragea2[which(is.na(current_file$coveragea2))] = 0
      result_list[[i]][[j]]$sigma02[k] = mean(current_file$sigma02)
      result_list[[i]][[j]]$sigma02_est.sd[k] = sqrt(var_file$v02[k])
      result_list[[i]][[j]]$sigma02_emp.sd[k] = sqrt(var(current_file$sigma02))
      result_list[[i]][[j]]$sigma02_coverage[k] = mean(current_file$coverage02)
      result_list[[i]][[j]]$sigmaa2[k] = mean(current_file$sigmaa2)
      result_list[[i]][[j]]$sigmaa2_est.sd[k] = sqrt(var_file$v12[k])
      result_list[[i]][[j]]$sigmaa2_emp.sd[k] = sqrt(var(current_file$sigmaa2))
      result_list[[i]][[j]]$sigmaa2_coverage[k] = mean(current_file$coveragea2)
      result_list[[i]][[j]]$omega[k] = mean(current_file$omega)
      result_list[[i]][[j]]$omega_emp.sd[k] = sqrt(var(current_file$omega))
      for(o in 1:length(current_file$omega)){
        if(current_file$omega[o] > 1){
          current_file$omega[o] = 1
        } else if(current_file$omega[o] < 0){
          current_file$omega[o] = 0
        }
      }
      result_list[[i]][[j]]$omega_trunc[k] = mean(current_file$omega)
      result_list[[i]][[j]]$omega_trunc_emp.sd[k] = sqrt(var(current_file$omega))
    }
    
    #write.csv(result_list[[i]][[j]], paste0(path, "MoL_result/general/select/", d, "_w_", w, "_variance.csv"),
    #          quote=F, row.names = F, col.names=T)
    write.csv(result_list[[i]][[j]], paste0(path, "MoL_result/general/select/fix_", d, "_w_", w, "_variance.csv"),
             quote=F, row.names = F, col.names=T)
  }
}

# MoL and PQL Heritability and computational time
library(data.table)

options(stringsAsFactors = F)
d_l = c("mix")

w_l = c(0.5)
h_l = c(0.4)
hreal_l = c(0.25)
sigma02_l = c(0.6)
sigmaa2_l = c(0.2)

#n_l_pql = c(500, 1000, 2000, 5000)
#n_l_mol = c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)
#p_l_pql = rep(5000, length(n_l_pql))
#p_l_mol = rep(5000, length(n_l_mol))

n_l_pql = c(500, 1000, 2000, 5000)
n_l_mol = c(500, 1000, 2000, 5000, 5000, 10000, 20000, 50000, 50000, 100000, 200000, 500000)
p_l_pql = rep(c(500,1000,2000,5000),1)
p_l_mol = rep(c(500,1000,2000,5000),3)
  
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
    MoL_file_list = paste0(path, "MoL_result/general/",  d, "_n_", n_l_mol, "_p_", p_l_mol, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    MoL_summary_list = lapply(MoL_file_list, fread)
    MoL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,100 * length(n_l_mol)), method = rep("MoL",100 * length(n_l_mol)),
                                           n = rep(n_l_mol,each = 100), p = rep(p_l_mol,each = 100),
                                           heritability = rep(NA,100 * length(n_l_mol)), comp_time = rep(NA,100 * length(n_l_mol)))
    for (k in 1:length(MoL_summary_list)){
      MoL_current_file = MoL_summary_list[[k]]
      MoL_result_list[[i]][[j]]$heritability[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$heritability
      if(c("time") %in% names(MoL_current_file)){
        MoL_result_list[[i]][[j]]$comp_time[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$time 
      } else{
        MoL_result_list[[i]][[j]]$comp_time[(100 * (k-1) + 1): (100 * k)] = MoL_current_file$time1
      }
    }
    
    PQL_file_list = paste0(path, "PQL_result/",  d, "_n_", n_l_pql, "_p_", p_l_pql, "_w_", w,
                           "_h_", h, "_N_", N, "_mu_", mu, "_summary.txt")
    PQL_summary_list = lapply(PQL_file_list, fread)
    PQL_result_list[[i]][[j]] = data.table(true_heritability = rep(hreal,100 * length(n_l_pql)), method = rep("PQLseq",100 * length(n_l_pql)),
                                           n = rep(n_l_pql,each = 100), p = rep(p_l_pql,each = 100),
                                           heritability = rep(NA,100 * length(n_l_pql)), comp_time = rep(NA,100 * length(n_l_pql)))
    for (k in 1:length(PQL_summary_list)){
      PQL_current_file = PQL_summary_list[[k]]
      PQL_result_list[[i]][[j]]$heritability[(100 * (k-1) + 1): (100 * k)] = PQL_current_file$h2
      PQL_result_list[[i]][[j]]$comp_time[(100 * (k-1) + 1): (100 * k)] = PQL_current_file$time
    }
    
    overall_result = rbind(MoL_result_list[[i]][[j]], PQL_result_list[[i]][[j]])
    
    #write.table(overall_result, paste0(path, "overall_result/select/", d, "_", hreal, "_heritability.txt"),
    #            quote=F, row.names = F, col.names=T)
    write.table(overall_result, paste0(path, "overall_result/select/fix_", d, "_", hreal, "_heritability.txt"),
                quote=F, row.names = F, col.names=T)
  }
}