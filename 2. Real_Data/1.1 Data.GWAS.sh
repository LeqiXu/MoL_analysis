## select age, sex and 10 PCs as covariates
cp /gpfs/ysm/pi/zhao-data/yy496/ukbb_pheno/MR_time/adjust/adjust.csv /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/adjust.csv

#1568   Average weekly red wine intake  Alcohol
#1578   Average weekly champagne plus white wine intake Alcohol

## 1. Phenotype preparation
## See the distribution of 8 phenotypes
setwd("/Users/xuleqi/Desktop/Research/Project/3.1 Count_Data_GWAS_Method：MoL/Data")
library(data.table)
library(tidyr)
library(ggplot2)

cov_data <- fread("covariates.csv")
cov_data = cov_data[which(cov_data$kin == 0 & cov_data$eur == 1),]

count_data = fread("count_data.csv")
count_data = count_data[which(count_data$eid %in% cov_data$eid),]
pheno_data = count_data[,-c(1:3)]
colnames(pheno_data) = c("week_red_wine","week_champagne")
colSums(!is.na(pheno_data) & pheno_data > 0)
pheno_data_gather = gather(pheno_data)
pheno_data_gather = pheno_data_gather[which(pheno_data_gather$value > 0),]
ggplot(pheno_data_gather,aes(value)) + geom_histogram() + facet_wrap(~key, scales = 'free', ncol = 5) + 
        theme(axis.text = element_text(size = 16),strip.text = element_text(size = 16), axis.title = element_text(size = 18))

## Extract hm3 SNPs
hm3_origianl <- read.table("/gpfs/gibbs/pi/zhao/lx94/trans_PRS/data/eur_w_ld_chr/w_hm3.snplist", header=T, stringsAsFactors=F, fill = TRUE)
hm3_SNP = data.frame(SNP = hm3_origianl$SNP)
write.table(hm3_SNP, file="/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/w_hm3.snplist", append=F, sep="\t", quote=F, row.names=F, col.names=T)

## data cleaning for each trait
# covariates.csv: We restrict our analysis to unrelated European population (kin==0,eur==1), 
# when conducting GWAS analysis, sex age_recruit and first 10 PCs should be adjusted.
library(data.table)

setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data")
cov_data <- fread("covariates.csv")
cov_data = cov_data[which(cov_data$kin == 0 & cov_data$eur == 1), c("eid","age_recruit","sex",paste0("PC",1:10))]

count_data = fread("count_data.csv")
colnames(count_data) = c("eid","age_recruit","date_recruit","week_red_wine","week_champagne")

## week_red_wine: 133610
week_red_wine_data <- count_data[,c("eid","age_recruit","week_red_wine")]
week_red_wine_data = na.omit(week_red_wine_data)
week_red_wine_data = week_red_wine_data[which(week_red_wine_data$week_red_wine > 0 & week_red_wine_data$eid %in% cov_data$eid),] 

keep_week_red_wine = data.frame(FID = week_red_wine_data$eid, IID = week_red_wine_data$eid)
pheno_week_red_wine = data.frame(FID = week_red_wine_data$eid, IID = week_red_wine_data$eid, week_red_wine = week_red_wine_data$week_red_wine)
cov_week_red_wine = cov_data[which(cov_data$eid %in% week_red_wine_data$eid)]
cov_week_red_wine$FID = cov_week_red_wine$eid
cov_week_red_wine$IID = cov_week_red_wine$eid
cov_week_red_wine$sex = ifelse(cov_week_red_wine$sex == "Male",1,2)
cov_week_red_wine = cov_week_red_wine[,c("FID","IID","sex","age_recruit",paste0("PC",1:10))]

fwrite(keep_week_red_wine,"keep_week_red_wine.txt", sep = "\t")
fwrite(pheno_week_red_wine,"pheno_week_red_wine.txt", sep = "\t")
fwrite(cov_week_red_wine,"cov_week_red_wine.txt", sep = "\t")

## week_champagne: 111351
week_champagne_data <- count_data[,c("eid","age_recruit","week_champagne")]
week_champagne_data = na.omit(week_champagne_data)
week_champagne_data = week_champagne_data[which(week_champagne_data$week_champagne > 0 & week_champagne_data$eid %in% cov_data$eid),] 

keep_week_champagne = data.frame(FID = week_champagne_data$eid, IID = week_champagne_data$eid)
pheno_week_champagne = data.frame(FID = week_champagne_data$eid, IID = week_champagne_data$eid, week_champagne = week_champagne_data$week_champagne)
cov_week_champagne = cov_data[which(cov_data$eid %in% week_champagne_data$eid)]
cov_week_champagne$FID = cov_week_champagne$eid
cov_week_champagne$IID = cov_week_champagne$eid
cov_week_champagne$sex = ifelse(cov_week_champagne$sex == "Male",1,2)
cov_week_champagne = cov_week_champagne[,c("FID","IID","sex","age_recruit",paste0("PC",1:10))]

fwrite(keep_week_champagne,"keep_week_champagne.txt", sep = "\t")
fwrite(pheno_week_champagne,"pheno_week_champagne.txt", sep = "\t")
fwrite(cov_week_champagne,"cov_week_champagne.txt", sep = "\t")

## 2. GWAS study: week_red_wine,week_champagne
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=120G
#SBATCH --cpus-per-task=30
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=gwas_week_red_wine
#SBATCH --output=out_gwas_week_red_wine.txt

## load modules in mccleary ##
module load PLINK/2_avx2_20221024

## run plink2 in exposure ##
plink2 --bfile /gpfs/ysm/pi/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal \
--double-id \
--extract /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/w_hm3.snplist \
--keep /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/keep_week_red_wine.txt \
--pheno /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/pheno_week_red_wine.txt \
--covar /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/cov_week_red_wine.txt \
--glm hide-covar \
--out /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm/week_red_wine

#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=120G
#SBATCH --cpus-per-task=30
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=gwas_week_champagne
#SBATCH --output=out_gwas_week_champagne.txt

## load modules in mccleary ##
module load PLINK/2_avx2_20221024

## run plink2 in exposure ##
plink2 --bfile /gpfs/ysm/pi/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal \
--double-id \
--extract /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/w_hm3.snplist \
--keep /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/keep_week_champagne.txt \
--pheno /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/pheno_week_champagne.txt \
--covar /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/cov_week_champagne.txt \
--glm hide-covar \
--out /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm/week_champagne

## 3. Organize GWAS
## select GWAS column
setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm")
library(data.table)

for (trait in c("week_red_wine","week_champagne")){
glm_data <- fread(paste0(trait,".",trait,".glm.linear"))
glm_data = na.omit(glm_data)
sumstat = glm_data[,c(3,6,4,8,12,9)]
colnames(sumstat) = c("SNP","A1","A2","N","P","BETA")
fwrite(sumstat, paste0(trait,"_glm.txt"),sep = "\t")
}

## clean the gwas result
module load miniconda
conda activate ldsc
cd /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data

# merge sumstat
for trait in week_red_wine week_champagne
do
python /gpfs/gibbs/pi/zhao/lx94/software/ldsc/munge_sumstats.py \
--sumstats glm/${trait}_glm.txt \
--out glm/${trait} \
--merge-alleles /gpfs/gibbs/pi/zhao/lx94/trans_PRS/data/eur_w_ld_chr/w_hm3.snplist
done

## 4. Obtain clean GWAS
# Previous smoke pos
setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm")
library(data.table)
library(readr)

for (trait in c("week_red_wine","week_champagne")){
glm_data <- fread(paste0(trait,"_glm.txt"))
clean_snp = read_tsv(paste0(trait,".sumstats.gz"))
clean_snp = na.omit(clean_snp)
glm_data = glm_data[which(glm_data$SNP %in% clean_snp$SNP),c("SNP","N","P","BETA")]
clean_snp = clean_snp[,c("SNP","A1","A2")]

final_data = merge(glm_data,clean_snp,by=c("SNP"))
final_data = na.omit(final_data)
fwrite(final_data, paste0("clean_",trait,"_glm.txt"),sep = "\t")
}