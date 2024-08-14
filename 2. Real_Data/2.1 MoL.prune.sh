# MoL estimation
module load miniconda
conda activate r_env

## Final data
setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/prune")

options(scipen = 999)

library(genio)
library(data.table)

# Calculate the variance results for each gender
MoL_each <- function(z,x,y,N_i){
    n = nrow(z)
    p = ncol(z)

    zy = t(z) %*% y
    xy = t(x) %*% y
    sumy = sum(y)
    sumy2 = sum(y**2)
    sumN = sum(N_i)
    sumN2 = sum(N_i**2)

    # sigma calculation
    sigmaa2 = (sum(zy**2) - sum(rowSums(z**2) * (y**2))) / (sumy**2)
    tao2 = (sum(xy**2) - sum(rowSums(x**2) * (y**2))) / (sumy**2)
    logS2 = log(sumy2-sumy) - 2*log(sumy) + 2 * log(sumN) - log(sumN2)
    sigma02 = logS2 - sigmaa2 - tao2
    sigma2 = sigmaa2 + sigma02
    muhat = log(sumy) - log(sumN) - 0.5 * logS2

    # b and v02 calculation
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

    # Omega_hat calculation
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

    # result summary
    total_df = data.frame(sigma02 = sigma02, sigmaa2 = sigmaa2, tao2 = tao2, v02 = v02, v02n1 = v02n1, v02n2 = v02n2, b1 = b1, b2 = b2, b3 = b3,T1 = T2, T2 = T2)
    return(total_df)
}

total_df_all = c()

# Data
thres=0.001
dist=100
    
for (trait in c("week_champagne","week_red_wine")){
geno <- read_plink(paste0(trait,"_prune_",thres,"_",dist))
cov <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/cov_",trait,".txt"))
pheno <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/pheno_",trait,".txt"))

geno_matrix = t(geno$X)
cov = cov[which(cov$FID %in% rownames(geno_matrix)),]
pheno = pheno[which(pheno$FID %in% rownames(geno_matrix)),]

# Geno matrix
col_m = apply(geno_matrix, 2, mean, na.rm = TRUE)
col_sd = apply(geno_matrix, 2, sd, na.rm = TRUE)

col_m_matrix = matrix(rep(col_m, nrow(geno_matrix)), ncol = length(col_m), byrow = TRUE)
col_sd_matrix = matrix(rep(col_sd, nrow(geno_matrix)), ncol = length(col_sd), byrow = TRUE)
sd_geno_matrix = (geno_matrix - col_m_matrix) / col_sd_matrix
sd_geno_matrix[which(is.na(sd_geno_matrix))] = 0

# Covaraites: continuous covariate: age + 10 PCs, discrete covaraites: sex
## Continous covariates should be standardized
continuous_cov = cov[,-c(1:3)]
continuous_cov = as.matrix(continuous_cov)
rownames(continuous_cov) = cov$FID

sd_continuous_cov = apply(continuous_cov,2,scale,scale = FALSE)
continuous_cov_eigen = eigen(cov(continuous_cov))
continuous_cov_eigenhalf = continuous_cov_eigen$vectors %*% diag(1 / sqrt(continuous_cov_eigen$values))
sd_continuous_cov =  sd_continuous_cov %*% continuous_cov_eigenhalf
rownames(sd_continuous_cov) = cov$FID

## We separate the data by male and female.
N = 1
discrete_cov = cov[,3]
maleID = which(discrete_cov$sex == 1) 
femaleID = which(discrete_cov$sex == 2)
n_male = length(maleID)
n_female = length(femaleID)
n_total = n_male + n_female

sd_geno_matrix_male = sd_geno_matrix[maleID,]
sd_geno_matrix_female = sd_geno_matrix[femaleID,]
sd_continuous_cov_male = sd_continuous_cov[maleID,]
sd_continuous_cov_female = sd_continuous_cov[femaleID,]
pheno_male = unlist(pheno[maleID,3])
pheno_female = unlist(pheno[femaleID,3])
N_male = rep(N,length(maleID))
N_female = rep(N,length(femaleID))

# Separate MoL calculation
male_result = MoL_each(z = sd_geno_matrix_male, x = sd_continuous_cov_male, y = pheno_male, N_i = N_male)
female_result = MoL_each(z = sd_geno_matrix_female, x = sd_continuous_cov_female, y = pheno_female, N_i = N_female)

# Combine two result together
sigma02 = (n_male * male_result$sigma02 + n_female * female_result$sigma02) / n_total
sigmaa2 = (n_male * male_result$sigmaa2 + n_female * female_result$sigmaa2) / n_total
sigma2 = sigmaa2 + sigma02
hhat = sigmaa2 / sigma2

gamma_male = n_male / ncol(sd_geno_matrix_male)
gamma_female =  n_female / ncol(sd_geno_matrix_female)
gamma_total = gamma_male + gamma_female

T1_square = (n_male * male_result$T1^2 + n_female * female_result$T1^2) / n_total
T2 = (n_male * male_result$T2 + n_female * female_result$T2) / n_total
omega_hat = 3 * T1_square / T2

v02 = (gamma_male * male_result$v02 + gamma_female * female_result$v02) / gamma_total
va2 = sigmaa2^2 * (gamma_total * (3 / omega_hat - 1) + 
                    4 / gamma_total * (gamma_male * (exp(sigma2 + male_result$tao2) * male_result$v02n1 + (male_result$b1 + male_result$b3) / male_result$b2) + 
                                       gamma_female * (exp(sigma2 + female_result$tao2) * female_result$v02n1 + (female_result$b1 + female_result$b3) / female_result$b2)))
v10 = 4 * exp(sigma2) * 1 / gamma_total *  (gamma_male * (sigmaa2 * exp(sigma2 + male_result$tao2) * male_result$v02n2 - (2 * sigmaa2^2 + male_result$tao2 + 1) * male_result$v02n1) + 
                                            gamma_female * (sigmaa2 * exp(sigma2 + female_result$tao2) * female_result$v02n2 - (2 * sigmaa2^2 + female_result$tao2 + 1) * female_result$v02n1))
v2 = (sigma2)^(-4) * (v02 * sigmaa2 ^ 2 - 2 * v10 * sigmaa2 * sigma02 + va2 * sigma02 ^ 2)

#if(omega_hat >=1){
#    omega_hat = 1
#} else if(omega_hat <= 0){
#    omega_hat = 0
#}

total_df = data.frame(trait = trait, threshold = thres, distance = dist, sample_size = nrow(sd_geno_matrix), n_SNPs = ncol(sd_geno_matrix), sigma02 = sigma02, sigmaa2 = sigmaa2, v_s02 = 1 / n_total * v02, v_sa2 = 1 / n_total * va2, heritability = hhat, v_h2 = 1 / n_total * v2, omega = omega_hat)
print(paste0(trait,"_",dist))
print(total_df)

total_df_all = rbind(total_df_all,total_df)
}

print(total_df_all)
write.table(total_df_all, file="/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/result/MoL_prune_table.tsv", append=F, sep="\t", quote=F, row.names=F, col.names=T)

## week_champagne, week_red_wine