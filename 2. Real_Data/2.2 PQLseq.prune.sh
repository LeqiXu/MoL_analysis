# 1. command generating
sizelist = c("s","l","ll")

dist=100
thres=0.001

for (trait in c("week_champagne","week_red_wine")){

for(i in 1:length(sizelist)){
joblist = c()
joblist = c(joblist, paste0("module load miniconda;conda activate r_env;Rscript --vanilla PQLseq_prune_",dist,"_",thres,"_",trait,"_",sizelist[i],".R ", 1:100))
write.table(joblist, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/code/PQLseq_prune_",dist,"_",thres,"_",trait,"_",sizelist[i],".sh"), quote=F, col.names = F, row.names = F)
}
}

#2. dsq generating
## choice 1
module load dSQ
for trait in week_champagne week_red_wine
do
dsq --job-file PQLseq_prune_${dist}_${thres}_${trait}_${size}.sh --job-name dsq-dsq-PQLseq_prune_${dist}_${thres}_${trait}_${size} --requeue --partition=bigmem -n 1 --cpus-per-task=1 --mem=100g -t 1-00:00:00 --mail-type=ALL --status-dir /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/code/joblist --output /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/code/joblist/dsq-PQLseq_prune_${dist}_${thres}_${trait}_${size}.out 
done

## choice 2
#!/bin/bash
#SBATCH --output /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/code/joblist/dsq-PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE.out
#SBATCH --array 0-99
#SBATCH --job-name dsq-dsq-PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE
#SBATCH --requeue --partition=bigmem -n 1 --cpus-per-task=1 --mem=100g -t 1-00:00:00 --mail-type=ALL

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/code/PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE.sh --status-dir /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/code/joblist

vim PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE.pbs

size=l
dist=100
thres=0.001

for trait in week_champagne week_red_wine
do
cp PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE.R PQLseq_prune_${dist}_${thres}_${trait}_${size}.R
sed -i "s/THESIZE/${size}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.R
sed -i "s/THETRAIT/${trait}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.R
sed -i "s/THETHRES/${thres}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.R
sed -i "s/THEDIST/${dist}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.R
done

for trait in week_champagne week_red_wine
do
cp PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE.pbs PQLseq_prune_${dist}_${thres}_${trait}_${size}.pbs
sed -i "s/THESIZE/${size}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.pbs; sed -i "s/THETRAIT/${trait}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.pbs
sed -i "s/THETHRES/${thres}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.pbs; sed -i "s/THEDIST/${dist}/g" PQLseq_prune_${dist}_${thres}_${trait}_${size}.pbs
done

for trait in week_champagne week_red_wine
do 
sbatch PQLseq_prune_${dist}_${thres}_${trait}_${size}.pbs
done

## Final data
setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/prune")

library(genio)
library(data.table)
library(PQLseq)

options(stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
j = as.numeric(args[1])

geno <- read_plink("THETRAIT_prune_THETHRES_THEDIST")
cov <- fread("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/cov_THETRAIT.txt")
pheno <- fread("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/pheno_THETRAIT.txt")
id <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/keep_THETRAIT_THESIZE_",j,".txt"))

cov = cov[which(cov$FID %in% id$FID),]
pheno = pheno[which(pheno$FID %in% id$FID),]

# Geno matrix
geno_matrix = t(geno$X)
geno_matrix = geno_matrix[which(rownames(geno_matrix) %in% id$FID),]
cov = cov[which(cov$FID %in% rownames(geno_matrix)),]
pheno = pheno[which(pheno$FID %in% rownames(geno_matrix)),]
n = nrow(pheno)
print(n)

col_m = apply(geno_matrix, 2, mean, na.rm = TRUE)
col_sd = apply(geno_matrix, 2, sd, na.rm = TRUE)

col_m_matrix = matrix(rep(col_m, n), ncol = length(col_m), byrow = TRUE)
col_sd_matrix = matrix(rep(col_sd, n), ncol = length(col_sd), byrow = TRUE)
sd_geno_matrix = (geno_matrix - col_m_matrix) / col_sd_matrix
sd_geno_matrix[which(is.na(sd_geno_matrix))] = 0

relatednessmatrix = sd_geno_matrix %*% t(sd_geno_matrix) / sum(diag(sd_geno_matrix %*% t(sd_geno_matrix))) * n

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

# Phneotype
y = unlist(pheno[,3])
dfy = data.frame(t(y))

# Predictor
predictor = rnorm(n)

N = 7
N_i = rep(N, n)
totalcount = data.frame(t(N_i))

## We separate the data by male and female.
discrete_cov = as.factor(cov$sex)
names(discrete_cov) = cov$FID
covariate = cbind(discrete_cov,sd_continuous_cov)

# PQLseq method
model_DNA = pqlseq(RawCountDataSet = dfy, Covariates = covariate,
                   Phenotypes = predictor,
                   RelatednessMatrix = relatednessmatrix,
                   LibSize = totalcount)

write.table(model_DNA, paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/result/PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE_",j, "_result.txt"), quote=F, col.names = T, row.names = F)

vim PQLseq_prune_THEDIST_THETHRES_THETRAIT_THESIZE.R