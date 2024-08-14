nlist = rep(c(10000, 7500, rep(5000, 5), 2500, 2000), 2)
plist = rep(c(rep(2000, 3), c(5000, 3500, 1500, 1000), rep(2000, 2)), 2)
hlist = c(rep(0.3, 9), rep(0.6, 9))

for(j in 1:length(nlist)){
  joblist = c()
  joblist = c(joblist, paste0("Rscript --vanilla Simplest_Mol.R ",
                              nlist[j], " ", plist[j], " ",
                              hlist[j], " ", 1:100))
  write.table(joblist, paste0("/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/joblist/n_",
                              nlist[j], "_p_", plist[j], "_h_", hlist[j], ".sh"), 
              quote=F, col.names = F, row.names = F)
}

job = c()
for(j in 1:length(nlist)){
  file_name = paste0("n_", nlist[j], "_p_", plist[j], "_h_", hlist[j])
  joblist = c("module load GCC;", "partition=scavenge;", "module load dSQ",
              paste0("dSQ --jobfile ", file_name, ".sh -p ${partition} -n 1 ",
                     "--constraint avx2 ",
              "--mem-per-cpu=32g -t 24:00:00 --mail-type=ALL --batch-file ",
              file_name, ".pbs"), paste0("sbatch ", file_name, ".pbs"))
  write.table(joblist, paste0("/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/joblist/",
                              file_name, ".submit"),
              quote=F, col.names = F, row.names = F)
  job = c(job, paste0("sh ", file_name, ".submit"))
}
write.table(job, "/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/joblist/submit_all.sh",
            quote=F, row.names = F, col.names = F)


################################################################################

nlist = rep(c(10000, 7500, rep(5000, 5), 2500, 2000), 2)
plist = rep(c(rep(2000, 3), c(5000, 3500, 1500, 1000), rep(2000, 2)), 2)
hlist = c(rep(0.3, 9), rep(0.6, 9))
N = 10
for(j in 1:length(nlist)){
  joblist = c()
  joblist = c(joblist, paste0("if [ ! -f \"/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/PQL/",
                              "n_", nlist[j], "_p_", plist[j], "_h_", hlist[j], 
                              "_N_", N, "/", 1:100, ".txt\" ];then ",
                              "Rscript --vanilla Simplest_PQLseq.R ",
                              nlist[j], " ", plist[j], " ",
                              hlist[j], " ", 1:100, ";fi;"))
  write.table(joblist, paste0("/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/joblist/pql_n_",
                              nlist[j], "_p_", plist[j], "_h_", hlist[j], ".sh"), 
              quote=F, col.names = F, row.names = F)
}

job = c()
for(j in 1:length(nlist)){
  file_name = paste0("pql_n_", nlist[j], "_p_", plist[j], "_h_", hlist[j])
  joblist = c("module load GCC;", "partition=scavenge;", "module load dSQ",
              paste0("dSQ --jobfile ", file_name, ".sh -p ${partition} -n 1 ",
                     "--constraint avx2 ",
                     "--mem-per-cpu=32g -t 24:00:00 --mail-type=ALL --batch-file ",
                     file_name, ".pbs"), paste0("sbatch ", file_name, ".pbs"))
  write.table(joblist, paste0("/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/joblist/",
                              file_name, ".submit"),
              quote=F, col.names = F, row.names = F)
  job = c(job, paste0("sh ", file_name, ".submit"))
}
write.table(job, "/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/joblist/submit_all_pql.sh",
            quote=F, row.names = F, col.names = F)
