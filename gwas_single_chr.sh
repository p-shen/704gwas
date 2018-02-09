module load gcc/6.2.0
module load R/3.4.1

Rscript gwas.r -a "WTCCC/58C/Affx_gt_58C_Chiamo_$1.tped.gz" -b "WTCCC/NBS/Affx_gt_NBS_Chiamo_$1.tped.gz" -d "WTCCC/T2D/Affx_gt_T2D_Chiamo_$1.tped.gz" -n $1 -s "WTCCC/58C/snps/snps_$1"
