######################################## MoL ########################################
#record: 5e+05, 2000: 50GB
#record: 3e+05, 10000: 70GB
#record: 3e+05, 2000: 20GB
#record: 2e+05, 10000: 40GB
#record: 1e+05, 10000: 40GB
#record: 1e+05, 5000: 20GB

# generate original
dlist = c(rep("norm",100), rep("mix",100))
nlist = rep(c(rep(200,20), rep(500,20), rep(1000,20), rep(2000,20),rep(5000,20)), 2)
plist = rep(c(rep(500, 4), rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

# generate more
dlist = c(rep("norm",100), rep("mix",100))
nlist = rep(c(rep(10000,20), rep(15000,20), rep(20000,20), rep(25000,20),rep(30000,20)), 2)
plist = rep(c(rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4), rep(20000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

# generate most
dlist = c(rep("norm",100), rep("mix",100))
nlist = rep(c(rep(50000,20), rep(100000,20), rep(200000,20), rep(300000,20),rep(500000,20)), 2)
plist = rep(c(rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4), rep(20000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

# generate one
# p = 500, n = 50000
dlist = c("mix")
nlist = c(50000)
plist = c(500)
wlist = c(0.5)
hlist = c(0.4)
N = 10
mu = 0.2

for(j in 1:length(dlist)){
  joblist = c()
  joblist = c(joblist, paste0("Rscript --vanilla General_MoL.R ",
                              dlist[j], " ", nlist[j], " ", 
                              plist[j], " ", wlist[j], " ",
                              hlist[j], " ", 1:200))
  write.table(joblist, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/MoL_", 
                              dlist[j], "_n_", nlist[j], "_p_", plist[j], "_w_", wlist[j],
                              "_h_", hlist[j], "_N_", N, "_mu_", mu, ".sh"),
              quote=F, col.names = F, row.names = F)
}

job = c()
for(j in 1:length(dlist)){
  file_name = paste0("MoL_", dlist[j], "_n_", nlist[j], "_p_", plist[j], "_w_", wlist[j],
                     "_h_", hlist[j], "_N_", N, "_mu_", mu)
  joblist = c("module load GCC;", "partition=scavenge;", "module load dSQ",
              paste0("dSQ --jobfile ", file_name, ".sh -p ${partition} -n 1 ",
                     "--constraint avx2 ",
                     "--mem-per-cpu=180g -t 3-00:00:00 --mail-type=ALL --batch-file ",
                     file_name, ".pbs"), paste0("sbatch ", file_name, ".pbs"))
  write.table(joblist, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/",
                              file_name, ".submit"),
              quote=F, col.names = F, row.names = F)
  job = c(job, paste0("sh ", file_name, ".submit"))
}

for(j in 1:length(dlist)){
  file_name = paste0("MoL_", dlist[j], "_n_", nlist[j], "_p_", plist[j], "_w_", wlist[j],
                     "_h_", hlist[j], "_N_", N, "_mu_", mu)
  joblist = c("module load GCC;", "partition=bigmem;", "module load dSQ",
              paste0("dSQ --jobfile ", file_name, ".sh -p ${partition} -n 1 ",
                     "--constraint avx2 ",
                     "--mem-per-cpu=500g -t 1-00:00:00 --mail-type=ALL --batch-file ",
                     file_name, ".pbs"), paste0("sbatch ", file_name, ".pbs"))
  write.table(joblist, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/",
                              file_name, ".submit"),
              quote=F, col.names = F, row.names = F)
  job = c(job, paste0("sh ", file_name, ".submit"))
}

write.table(job, "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/submit_all_MoL_most.sh",
            quote=F, row.names = F, col.names = F)


######################################## PQL seq ########################################
#record: 10000, 1000: 15GB 10h
#record: 10000, 2000: 15GB 5h
#record: 20000, 1000: 80GB 5 CPU 21h
# We are going to run the results with 10 CPUs for all the jobs
# now only 20000
# generate original
dlist = c(rep("norm",100), rep("mix",100))
nlist = rep(c(rep(200,20), rep(500,20), rep(1000,20), rep(2000,20),rep(5000,20)), 2)
plist = rep(c(rep(500, 4), rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

# generate more
dlist = c(rep("norm",100), rep("mix",100))
nlist = rep(c(rep(10000,20), rep(15000,20), rep(20000,20), rep(25000,20),rep(30000,20)), 2)
plist = rep(c(rep(1000, 4), rep(2000, 4), rep(5000, 4), rep(10000, 4), rep(20000, 4)), 10)
wlist = rep(c(0.15, 0.5, 0.0625, 0.375),50)
hlist = rep(c(0.4, 0.4, 0.8, 0.8), 50)
N = 10
mu = 0.2

# generate most
dlist = rep("norm",2)
nlist = c(20000)
#nlist = c(50000,100000)
plist = rep(1000,2)
wlist = c(0.15,2)
hlist = c(0.4,2)
N = 10
mu = 0.2

for(j in 1:length(nlist)){
  joblist = c()
  joblist = c(joblist, paste0("Rscript --vanilla General_PQLseq.R ",
                              dlist[j], " ", nlist[j], " ", 
                              plist[j], " ", wlist[j], " ",
                              hlist[j], " ", 1:200))
  write.table(joblist, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/PQL_",
                              dlist[j], "_n_", nlist[j], "_p_", plist[j], "_w_", wlist[j],
                              "_h_", hlist[j], "_N_", N, "_mu_", mu, ".sh"), 
              quote=F, col.names = F, row.names = F)
}

job = c()
for(j in 1:length(dlist)){
  file_name = paste0("PQL_", dlist[j], "_n_", nlist[j], "_p_", plist[j], "_w_", wlist[j],
                     "_h_", hlist[j], "_N_", N, "_mu_", mu)
  joblist = c("module load GCC;", "partition=scavenge,pi_zhao,bigmem;", "module load dSQ",
              paste0("dSQ --jobfile ", file_name, ".sh -p ${partition} -n 1 ",
                     "--constraint avx2 ",
                     "--mem=150g --cpus-per-task=12 -t 3-00:00:00 --mail-type=ALL --batch-file ",
                     file_name, ".pbs"), paste0("sbatch ", file_name, ".pbs"))
  write.table(joblist, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/",
                              file_name, ".submit"),
              quote=F, col.names = F, row.names = F)
  job = c(job, paste0("sh ", file_name, ".submit"))
}
write.table(job, "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/joblist/submit_all_PQL_most.sh",
            quote=F, row.names = F, col.names = F)
