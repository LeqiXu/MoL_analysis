args = commandArgs(trailingOnly=TRUE)
library(data.table)
options(stringsAsFactors=F)
output_folder = "/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Simulation/General/MoL_result/general/"
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

start_time = Sys.time()
zy = t(z) %*% y
sumy = sum(y)
sumy2 = sum(y**2)
sumN = sum(N_i)
sumN2 = sum(N_i**2)
S2 = ((sumN**2) * (sumy2-sumy)) / ((sumy**2) * sumN2)

sigmaa2 = (sum(zy**2) - sum(rowSums(z**2) * (y**2))) / (sumy**2)
sigma02 = log(S2) - sigmaa2
sigma2 = sigma02+sigmaa2
muhat = log(sumy / n) - log(sumN / n) - 0.5 * log(S2)
hhat = sigmaa2 / (sigmaa2 + sigma02)

T2 = 0
for(t in 1:ncol(z)){
  z_sub = z[,t] * y
  s_sub1 = sum(z_sub)
  s_sub2 = sum(z_sub**2)
  s_sub3 = sum(z_sub**3)
  s_sub4 = sum(z_sub**4)
  T2_sub = (s_sub1**4) - 6 * (s_sub1**2) * s_sub2 + 8 * s_sub1 * s_sub3 + 3 * (s_sub2**2) - 6 * s_sub4
  T2 = T2 + T2_sub
}

T1 = (sum(zy**2) - sum(rowSums(z**2) * (y**2))) / (n**2)
T2 = p / (n * (n-1) * (n-2) * (n-3)) * T2
psi_hat = T2 / (3 * (T1**2))
omega_hat = 1 / psi_hat

En1 = 1 / n * sum(N_i)
En2 = 1 / n * sum(N_i^2)
En3 = 1 / n * sum(N_i^3)
En4 = 1 / n * sum(N_i^4)
v02n1 = En2 / (En1^2)
v02n2 = En3 / (En1 * En2)
v02n3 = En4 / (En2^2)
b1 = 1 / n * sumy
b2 = (sum(zy**2) - sum(rowSums(z**2) * (y**2))) / (n^2)
b3 = 1 / n * (sumy2 - sumy)

v02 = 4 * (((sigmaa2 + 1)^2 + sigmaa2) * exp(sigma2) - 1) * v02n1 - 
  4 * ((2 * sigmaa2 + 1) * exp(2 * sigma2) - 1) * v02n2 + 
  (exp(4 * sigma2) - 1) * v02n3 +
  4 / b1 * (exp(sigma2) * (En1 * En3) / (En2^2) - (sigmaa2 + 1)) + 2 / b3
v02 = 1 / n * v02
va2 = sigmaa2^2 * (n / p * (3 / omega_hat - 1) + 4 * (exp(sigma2) * v02n1 + (b1 + b3) / b2))
va2 = 1 / n * va2

true_sigma02 = 1 - h
true_sigmaa2 = h * w
  
coverage02 = ifelse((sigma02 <= true_sigma02 + 1.96 * sqrt(v02)) & (sigma02 >= true_sigma02 - 1.96 * sqrt(v02)),1,0)
coveragea2 = ifelse((sigmaa2 <= true_sigmaa2 + 1.96 * sqrt(va2)) & (sigmaa2 >= true_sigmaa2 - 1.96 * sqrt(va2)),1,0)

end_time = Sys.time()
total_time = as.numeric(end_time - start_time, unit = "secs")
df = data.frame(sigma02 = sigma02, sigmaa2 = sigmaa2, v02 = v02, va2 = va2, 
                coverage02 = coverage02, coveragea2 = coveragea2, 
                heritability = hhat, mu = muhat, omega = omega_hat, time = total_time)

write.table(df, paste0(i, ".txt"), quote=F, col.names = T, row.names = F)