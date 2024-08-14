# 1. Genotype pruning
# extract P less than 0.001 SNPs
library(data.table)

for(trait in c("week_red_wine","week_champagne")){
clean_glm = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm/clean_",trait,"_glm.txt"))
SNP_0.001 = data.table(SNP=clean_glm$SNP[which(clean_glm$P <=0.001)])
write.table(SNP_0.001, file=paste0("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm/",trait,"_SNP_0.001.txt"), append=F, sep="\t", quote=F, row.names=F, col.names=T)
}

#!/bin/bash
#SBATCH --requeue
#SBATCH --partition=bigmem,scavenge
#SBATCH --mem=150G
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=prune_THETRAIT_THETHRES_THEDIST
#SBATCH --output=out_prune_THETRAIT_THETHRES_THEDIST.txt

# mccleary
module load PLINK/1.9b_6.21-x86_64 
# grace
#module load PLINK/1.90-beta6.9

## prune for previous smoke pos##
trait=THETRAIT
thres=THETHRES
dist=THEDIST
plink --bfile /gpfs/ysm/pi/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal \
--keep /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/keep_${trait}.txt \
--extract /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/glm/${trait}_SNP_${thres}.txt \
--indep ${dist} 5 1 \
--out /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/prune/${trait}_prune_${thres}_${dist}

vim prune_THETRAIT_THETHRES_THEDIST.sh

thres=0.001
dist=100

for trait in week_red_wine week_champagne
do
cp prune_THETRAIT_THETHRES_THEDIST.sh prune_${trait}_${thres}_${dist}.sh
sed -i "s/THETRAIT/${trait}/g" prune_${trait}_${thres}_${dist}.sh
sed -i "s/THETHRES/${thres}/g" prune_${trait}_${thres}_${dist}.sh
sed -i "s/THEDIST/${dist}/g" prune_${trait}_${thres}_${dist}.sh
sbatch prune_${trait}_${thres}_${dist}.sh
done

## Organize SNPs
module load miniconda
conda activate r_env

setwd("/gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/prune")

thres=0.001
dist=100

for (trait in c("week_red_wine","week_champagne")){
prune = read.table(paste0(trait,"_prune_",thres,"_",dist,".prune.in"), header=F, stringsAsFactors=F)
prune_SNP = data.frame(SNP = prune$V1)
write.table(prune_SNP, file=paste0(trait,"_prune_",thres,"_",dist,".snplist"), append=F, sep="\t", quote=F, row.names=F, col.names=T)
}


## Extract prune SNPs
## load modules in mccleary ##
module load PLINK/2_avx2_20221024

## run plink2 for previous smoke pos ##
thres=0.001
dist=100

for trait in week_red_wine week_champagne
do
plink2 --bfile /gpfs/ysm/pi/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal \
--double-id \
--keep /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/ukb_data/keep_${trait}.txt \
--extract /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/prune/${trait}_prune_${thres}_${dist}.snplist \
--make-bed \
--out /gpfs/gibbs/pi/zhao/lx94/MoLGWAS/Real_data/result_data/prune/${trait}_prune_${thres}_${dist}
done