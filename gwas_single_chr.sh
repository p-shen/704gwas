module load gcc/6.2.0
module load R/3.4.1

Rscript gwas.r -a "WTCCC/58C/Affx_gt_58C_Chiamo_$2.tped.gz" -b "WTCCC/NBS/Affx_gt_NBS_Chiamo_$2.tped.gz" -d "WTCCC/$1/Affx_gt_$1_Chiamo_$2.tped.gz" -n $2 -s "WTCCC/58C/snps/snps_$2" -g $1
