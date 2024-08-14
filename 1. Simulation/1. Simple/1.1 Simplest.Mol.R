#args = commandArgs(trailingOnly=TRUE)
#library(data.table)

#options(stringsAsFactors=F)
#output_folder = "/ysm-gpfs/pi/zhao/yz738/MolGWAS/Simplest/MOL"
#n = as.numeric(args[1])
#p = as.numeric(args[2])
#h = as.numeric(args[3])
#i = as.numeric(args[4])
n = 10000
p = 5000
w = 1
h = 0.4
i = 160
N = 10

#setwd(output_folder)
#newfolder = paste0("./n_", n, "_p_", p, "_h_", h, "_N_", N)
#dir.create(newfolder)
#setwd(newfolder)

set.seed(i)

z = matrix(rnorm(n * p), nrow=n)
alpha = rnorm(p, mean = 0, sd = sqrt(h / p))
epsilon = rnorm(n, mean = 0, sd = sqrt(1 - h))
eta = z%*%alpha + epsilon
y = rpois(n, exp(eta)*N)

start_time = Sys.time()
zy = t(z) %*% y
sumy = sum(y)
sigma12 = (sum(zy**2) - sum(rowSums(z**2) * (y ** 2))) / (sumy**2)
sigma02 = 2 * (log(sumy / n) - log(N)) - sigma12
hhat = sigma12 / (sigma12 + sigma02)

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
T2 = p /((n-1)*(n-2)*(n-3)*n) * T2
psi_hat = T2/(3*(T1**2))

end_time = Sys.time()
total_time = as.numeric(end_time - start_time, unit="secs")

df = data.frame(T1 = T1, T2 = T2, omega = 1/psi_hat, 
                sigma02=sigma02, sigma12=sigma12, heritability=hhat, time=total_time)
print(round(df,3))

#write.table(df, paste0(i, ".txt"), quote=F, col.names = T, row.names = F)