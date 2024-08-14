args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(PQLseq)
options(stringsAsFactors=F)
output_folder = "/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/PQL"
n = as.numeric(args[1])
p = as.numeric(args[2])
h = as.numeric(args[3])
i = as.numeric(args[4])
N = 10

setwd(output_folder)
newfolder = paste0("./n_", n, "_p_", p, "_h_", h, "_N_", N)
dir.create(newfolder)
setwd(newfolder)

set.seed(i)

z = matrix(rnorm(n * p), nrow=n)
alpha = rnorm(p, mean = 0, sd = sqrt(h / p))
epsilon = rnorm(n, mean = 0, sd = sqrt(1 - h))
eta = z%*%alpha + epsilon
y = rpois(n, exp(eta)*N)
dfy = data.frame(t(y))
predictor = rnorm(n)
tildez = scale(z)
relatednessmatrix = tildez %*% t(tildez) / sum(diag(tildez %*% t(tildez))) * n
totalcount = data.frame(t(rep(N, n)))
  
start_time = Sys.time()
model_DNA = pqlseq(RawCountDataSet = dfy, Phenotypes = predictor,
                   RelatednessMatrix = relatednessmatrix,
                   LibSize = totalcount)
end_time = Sys.time()
total_time = as.numeric(end_time - start_time, unit="secs")
model_DNA$time = total_time

write.table(model_DNA, paste0(i, ".txt"), quote=F, col.names = T, row.names = F)
