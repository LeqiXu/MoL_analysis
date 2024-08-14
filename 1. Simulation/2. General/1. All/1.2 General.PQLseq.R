args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(PQLseq)
options(stringsAsFactors=F)
output_folder = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/PQL_result"
d = as.character(args[1])
n = as.numeric(args[2])
p = as.numeric(args[3])
w = as.numeric(args[4])
h = as.numeric(args[5])
i = as.numeric(args[6])
N = 10
mu = 0.2

setwd(output_folder)
newfolder = paste0("./", d, "_n_", n, "_p_", p, "_w_", w, "_h_", h, "_N_", N, "_mu_", mu)
dir.create(newfolder)
setwd(newfolder)

set.seed(i)

if(d == "norm"){
  z = matrix(rnorm(n * p, mean = 0), nrow=n)
} else{
  prob = rbeta(n*p, 0.5, 0.5)
  prob = 0.5 * prob
  z = rep(0, n*p)
  for (entry in 1:length(z)){
    z[entry] = rbinom(1,size = 2, prob = prob[entry])
  }
  z = (z - 0.5) / sqrt(0.4375)
  z = matrix(z, nrow=n)
}

b = rbinom(p, size = 1, prob = w)
xi = rnorm(p, mean = 0, sd = sqrt(h / p))
alpha = b * xi
epsilon = rnorm(n, mean = 0, sd = sqrt(1 - h))
eta = z%*%alpha + epsilon
N_i = rpois(n,N)
y = rep(N,n)
for (j in 1:n){
  y[j] = rpois(1, exp(eta[j] + mu)*N_i[j])
}

dfy = data.frame(t(y))
covariate = data.frame(intercept = rep(1,n))
predictor = rnorm(n)
tildez = scale(z)
relatednessmatrix = tildez %*% t(tildez) / sum(diag(tildez %*% t(tildez))) * n
totalcount = data.frame(t(N_i))

start_time = Sys.time()
model_DNA = pqlseq(RawCountDataSet = dfy, Covariates = covariate,
                   Phenotypes = predictor,
                   RelatednessMatrix = relatednessmatrix,
                   LibSize = totalcount)
end_time = Sys.time()
total_time = as.numeric(end_time - start_time, unit="secs")
model_DNA$time = total_time

write.table(model_DNA, paste0(i, ".txt"), quote=F, col.names = T, row.names = F)
